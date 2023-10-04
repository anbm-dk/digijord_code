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
    cov_keep = NULL, # Character vector, covariates that should always be present
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
  # Identify covariates that are not OGCs and make formula
  formula_basic <- cov_names %>%
    grep('ogc_pi', ., value = TRUE, invert = TRUE) %>%
    paste0(collapse = " + ") %>%
    paste0(target, " ~ ", .) %>%
    as.formula()
  # Basic model
  showConnections()
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  set.seed(321)
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
  showConnections()
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
  # Scoring function for Bayesian optimization
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
    # Drop unimportant covariates and add OGCs
    cov_i_filtered <- cov_ranked %>%
      filter(cumul < total_imp) %>%  #!
      .$rowname %>%
      c(., ogcs_names_list[[ogcs_index]])  # !
    # Make sure SOM removal is a covariate
    if (!is.null(cov_keep)) {
      cov_i_filtered %<>%
        c(., cov_keep) %>%
        unique()
    }
    # Make formula
    formula_i <- cov_i_filtered %>%
      paste0(collapse = " + ") %>%
      paste0(target, " ~ ", .) %>%
      as.formula()
    my_gamma <- gamma_sqrt^2
    my_min_child_weight <- min_child_weight_sqrt^2
    # Train model
    set.seed(321)
    model_out <- caret::train(
      form = formula_i,
      data = data,
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
        index = folds,
        savePredictions = "final",
        predictionBounds = bounds_pred,
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metric,
      maximize = max_metric,
      weights = weights,
      num_parallel_tree = trees_per_round,
      objective = objective,
      colsample_bylevel = colsample_bylevel,
      nthread = 1
    )
    if (max_metric) {
      out_score <- model_out$results %>%
        select(any_of(metric)) %>%
        max()
    } else {
      out_score <- model_out$results %>%
        select(any_of(metric)) %>%
        min() %>%
        "*"(-1)
    }
    return(
      list(
        Score = out_score,
        n_ogcs = length(ogcs_names_list[[ogcs_index]]),
        gamma = my_gamma,
        min_child_weight = my_min_child_weight,
        n_cov = length(cov_i_filtered)
      )
    )
  }
  # Bayesian optimization
  showConnections()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterEvalQ(
    cl,
    {
      require(caret)
      require(xgboost)
      require(magrittr)
      require(dplyr)
      require(tools)
      require(boot)
    }
  )
  clusterExport(
    cl,
    c(
      "target",
      "bounds_pred",
      "cov_i_ranked",
      "folds",
      "sumfun",
      "metric",
      "objective",
      "ogcs_names_list",
      "tgrid",
      "data",
      "trees_per_round"
    ),
    envir = environment()
  )
  set.seed(321)
  scoreresults <- bayesOpt(
    FUN = scoringFunction,
    bounds = bounds_bayes,
    initPoints = cores,
    iters.n = cores*10,
    iters.k = cores,
    acq = "ucb",
    gsPoints = cores*10,
    parallel = TRUE,
    verbose = 0,
    acqThresh = 0.95
  )
  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)
  showConnections()
  bestscores <- scoreresults$scoreSummary %>%
    filter(Score == max(Score, na.rm = TRUE))
  best_pars <- getBestPars(scoreresults)
}


# Final model
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
  ),
  envir = environment()
)

set.seed(321)

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