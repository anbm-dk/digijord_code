# Function to optimize xgboost

bounds <- list(
  eta = c(0.1, 1),
  max_depth = c(1L, 60L),
  min_child_weight_sqrt = c(1, sqrt(64)),
  gamma_sqrt = c(0, sqrt(40)),
  colsample_bytree = c(0.1, 1),
  subsample = c(0.1, 1),
  colsample_bylevel = c(0.1, 1),
  ogcs_index = c(1L, 7L),
  total_imp = c(0.5, 1),
  
)

optimize_xgboost <- function(
    target = NULL,  # character vector (length 1), target variable.
    cov_names = NULL,  # Character vector, covariate names,
    data = NULL, # data frame, input data
    bounds_bayes = NULL, # named list with bounds for bayesian opt.
    bounds_pred = NULL, # numeric, length 2, bounds for predicted values
    cores = 19, # number cores for parallelization
    trgrid = NULL, # data frame with tuning parameters to be tested in basic model
    folds = NULL, # list with indices, folds for cross validation
    sumfun = NULL, # summary function for accuracy assessment
    metric = NULL, # character, length 1, name of evaluation metric
    max_metric = NULL, # logical, should the evaluation metric be maximized
    weights = NULL, # numeric, weights for model training and evaluation
    trees_per_round = NULL, # numeric, length 1, number of trees that xgboost should train in each round
    objective = NULL, # character, length 1, objective function for xgboost
    colsample_bylevel_basic = 0.75, # numeric, colsample_bylevel for basic model
) {
  require(ParBayesianOptimization)
  require(caret)
  require(xgboost)
  require(magrittr)
  require(dplyr)
  # Identify OGCs and make a list with the numbers of OGCs to be tested in the
  # models
  ogcs_names <- cov_names %>%
    grep('ogc_pi', ., value = TRUE)
  ogcs_names_list <- list(ogcs_names)
  n_ogcs_v <- numeric()
  m <- 1
  n_ogcs <- length(ogcs_names_list[[m]])
  n_ogcs_v[m] <- n_ogcs
  while (n_ogcs > 2) {
    m <- m + 1
    ogcs_names_list[[m]] <- ogcs_names_list[[m - 1]][c(TRUE, FALSE)]
    n_ogcs <- length(ogcs_names_list[[m]])
    n_ogcs_v[m] <- n_ogcs
  }
  ogcs_names_list[[length(ogcs_names_list) + 1]] <- character()
  n_ogcs_v %<>% c(., 0)
  
  # Identify covariates that are not OGCs
  cov_p_i <- cov_names %>%
    grep('ogc_pi', ., value = TRUE, invert = TRUE) %>%
    paste0(collapse = " + ")
  
  # Basic model
  formula_basic <- paste0(target, " ~ ", cov_p_i) %>%
    as.formula()
  
  showConnections()
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, { require(boot) } )
  clusterExport(cl, c("get_RMSEw", "get_R2w"))
  
  set.seed(1)
  
  basic_model <- caret::train(
    form = formula_basic,
    data = data,
    method = "xgbTree",
    na.action = na.pass,
    tuneGrid = trgrid,
    trControl = trainControl(
      index = folds,
      savePredictions = "final",
      predictionBounds = bounds_pred,
      summaryFunction = sumfun,
      allowParallel = TRUE
    ),
    metric = metric,
    maximize = max_metric,
    weights = weights,
    num_parallel_tree = trees_per_round,
    objective = objective,
    colsample_bylevel = colsample_bylevel_basic,
    nthread = 1
  )
  
  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)
  
  tr_step <- 1
  tr_summaries <- list()
  tr_summaries[[tr_step]] <- basic_model$results
  tr_step %<>% `+`(1)
  
  # Scaled cumulative covariate importance
  cov_ranked <- basic_model %>%
    varImp(basic_model) %>%
    .$importance %>%
    rownames_to_column() %>%
    arrange(desc(Overall)) %>%
    mutate(
      scaled = Overall/sum(Overall),
      cumul = cumsum(scaled)
    )
  
  # Bayesian optimization
  scoringFunction <- function(
    eta,  # OK
    max_depth,  # OK
    min_child_weight_sqrt,  # OK
    gamma_sqrt,  # OK
    colsample_bytree,  # OK
    subsample,  # OK
    colsample_bylevel,
    ogcs_index,  # OK
    total_imp  # OK
  ) {
    # Drop unimportant covariates
    cov_i_filtered <- cov_ranked %>%
      filter(cumul < total_imp) %>%  #!
      .$rowname
    
    # Make sure SOM removal is a covariate
    if ((i %in% 1:4) & !("SOM_removed" %in% cov_i_filtered)) {
      cov_i_filtered %<>% c(., "SOM_removed")
    }
    
    # Add OGCs
    cov_i_filtered %<>% c(., ogcs_names_list[[ogcs_index]])  # !
    
    # Make formula
    cov_formula <- cov_i_filtered %>% paste0(collapse = " + ")
    
    formula_i <- paste0(frac, " ~ ", cov_formula) %>%
      as.formula()
    
    my_gamma <- gamma_sqrt^2
    my_min_child_weight <- min_child_weight_sqrt^2
    
    showConnections()
    
    set.seed(1)
    
    model_out <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = tgrid$nrounds,
        eta = eta,  # !
        max_depth = max_depth,  # !
        min_child_weight = my_min_child_weight, # !
        gamma = my_gamma, # !
        colsample_bytree = colsample_bytree, # !
        subsample = subsample # !
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = bounds_pred,
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metrics_i,
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives_i,
      colsample_bylevel = colsample_bylevel,
      nthread = 1
    )
    
    min_RMSEw <- model_out$results %>%
      select(any_of(metrics_i)) %>%
      min()
    
    return(
      list(
        Score = 0 - min_RMSEw,
        n_ogcs = length(ogcs_names_list[[ogcs_index]]),
        gamma = my_gamma,
        min_child_weight = my_min_child_weight,
        n_cov = length(cov_i_filtered)
      )
    )
  }
  
  
}

showConnections()

cl <- makeCluster(19)
registerDoParallel(cl)
clusterEvalQ(
  cl,
  {
    library(caret)
    library(xgboost)
    library(magrittr)
    library(dplyr)
    library(tools)
    library(boot)
  }
)

bounds_lower_i <- bounds_lower[i]
bounds_upper_i <- bounds_upper[i]
metrics_i <- metrics[i]
objectives_i <- objectives[i]

clusterExport(
  cl,
  c("i",
    "frac",
    "bounds_lower_i",
    "bounds_upper_i",
    "cov_i_ranked",
    "folds_i",
    "get_RMSEw",
    "get_R2w",
    "metrics_i",
    "objectives_i",
    "ogcs_names_list",
    "sumfun",
    "tgrid",
    "trdat",
    "trees_per_round"
  )
)

set.seed(321)

models_scoreresults[[i]] <- bayesOpt(
  FUN = scoringFunction,
  bounds = bounds_bayes,
  initPoints = 19,
  iters.n = 190,
  iters.k = 19,
  acq =  "ucb",
  gsPoints = 190,
  parallel = TRUE,
  verbose = 1,
  acqThresh = 0.95
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)

print(
  models_scoreresults[[i]]$scoreSummary
)

models_bestscores[[i]] <- models_scoreresults[[i]]$scoreSummary %>%
  filter(Score == max(Score, na.rm = TRUE))

print(
  models_bestscores[[i]]
)

best_pars <- getBestPars(models_scoreresults[[i]])

print("Training final model")

# Drop unimportant covariates
cov_i_filtered <- cov_i_ranked %>%
  filter(cumul < best_pars$total_imp) %>%  #!
  .$rowname

# Make sure SOM removal is a covariate
if (i %in% 1:4 & !"SOM_removed" %in% cov_i_filtered) {
  cov_i_filtered %<>% c(., "SOM_removed")
}

total_imp_models[i] <- best_pars$total_imp

# Add OGCs
cov_i_filtered %<>% c(., ogcs_names_list[[best_pars$ogcs_index]])  # !
n_ogcs_models[i] <- n_ogcs_v[best_pars$ogcs_index]

# Make formula
cov_formula <- cov_i_filtered %>% paste0(collapse = " + ")

formula_i <- paste0(frac, " ~ ", cov_formula) %>%
  as.formula()

# Lower eta

eta_test_final <- best_pars$eta %>%
  log() %>%
  seq(., . + log(0.01), length.out = 9) %>%
  exp() %>%
  round(3)

showConnections()
cl <- makePSOCKcluster(19)
registerDoParallel(cl)

clusterEvalQ(
  cl,
  {
    library(boot)
  }
)

clusterExport(
  cl,
  c(
    "get_RMSEw",
    "get_R2w"
  )
)

set.seed(1)

model_final <- caret::train(
  form = formula_i,
  data = trdat,
  method = "xgbTree",
  na.action = na.pass,
  tuneGrid = expand.grid(
    nrounds = tgrid$nrounds*10,
    eta = eta_test_final, # NB
    max_depth = best_pars$max_depth,
    min_child_weight = best_pars$min_child_weight_sqrt^2,
    gamma = best_pars$gamma_sqrt^2,
    colsample_bytree = best_pars$colsample_bytree,
    subsample = best_pars$subsample
  ),
  trControl = trainControl(
    index = folds_i,
    savePredictions = "final",
    predictionBounds = c(bounds_lower[i], bounds_upper[i]),
    summaryFunction = sumfun,
    allowParallel = TRUE
  ),
  metric = metrics[i],
  maximize = FALSE,
  weights = trdat$w,
  num_parallel_tree = trees_per_round,
  objective = objective,
  colsample_bylevel = best_pars$colsample_bylevel,
  nthread = 1
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)

# END