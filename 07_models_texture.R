# 07: Train texture models

# 1: Start up

# library(Cubist)
library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)
library(xgboost)

library(doParallel)
library(spatstat) # weights

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

source("f_predict_passna.R")

# To do:
# Pdp with depth
# Profile examples

# Test 1 - 8: Cubist
# Test 8: New covariates (chelsa, river valley bottoms, hillyness)
# Test 9: xgboost
# Test 10: Predicted row and column (poor accuracy)
# Test 11: Fever data, more
# Test 12: Depth boundaries as covariates, stepwise xgb optimization
# Test 13: Bayesian optimization
testn <- 13
mycrs <- "EPSG:25832"

# Results folder

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") %T>%
  dir.create()

# 2: Load observations

dir_obs_proc <- dir_dat %>%
  paste0(., "/observations/processed/")

dsc <- dir_obs_proc %>%
  paste0(., "dsc.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

SEGES <- dir_obs_proc %>%
  paste0(., "SEGES.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

SINKS <- dir_obs_proc %>%
  paste0(., "SINKS.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

profiles_texture <- dir_obs_proc %>%
  paste0(., "profiles_texture.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

forest_samples <- dir_obs_proc %>%
  paste0(., "forest_samples.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

# 3: Load folds

dir_folds <- dir_dat %>%
  paste0(., "/folds/")

dsc_folds <- dir_folds %>%
  paste0(., "dsc_folds.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )

SEGES_folds <- dir_folds %>%
  paste0(., "SEGES_folds.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )

SINKS_folds <- dir_folds %>%
  paste0(., "SINKS_folds.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )

profiles_folds <- dir_folds %>%
  paste0(., "profiles_folds.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  right_join(values(profiles_texture)) %>%
  select(lyr.1)

forest_folds <- dir_folds %>%
  paste0(., "forest_folds.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )

# 4: Load covariate data

dir_cov <- dir_dat %>% paste0(., "/covariates")

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230712.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

cov_files <- dir_cov %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_names %>%
  write.table(
    paste0("cov_names_", Sys.Date(), ".csv")
  )

cov_names[!cov_names %in% cov_cats$name]


# 5: Load extracted covariates

dir_extr <- dir_dat %>%
  paste0(., "/extracts/")

usebuffer <- FALSE

if (usebuffer) {
  dsc_extr <- dir_extr %>%
    paste0(., "/buffer_dsc_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )

  SEGES_extr <- dir_extr %>%
    paste0(., "/buffer_SEGES_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )
} else {
  dsc_extr <- dir_extr %>%
    paste0(., "/dsc_extr.rds") %>%
    readRDS()

  SEGES_extr <- dir_extr %>%
    paste0(., "/SEGES_extr.rds") %>%
    readRDS()
}

SINKS_extr <- dir_extr %>%
  paste0(., "/SINKS_extr.rds") %>%
  readRDS()

profiles_extr <- dir_extr %>%
  paste0(., "profiles_extr.rds") %>%
  readRDS() %>%
  right_join(values(profiles_texture)) %>%
  select(any_of(cov_names))

SINKS_extr <- dir_extr %>%
  paste0(., "/SINKS_extr.rds") %>%
  readRDS()

forests_extr <- dir_extr %>%
  paste0(., "/forests_extr.rds") %>%
  readRDS()

# 6: Merge data and transform the target variables

obs_data <- list(dsc, SEGES, SINKS, profiles_texture, forest_samples) %>%
  vect() %>%
  values() %>%
  mutate(
    logSOC = log(SOC),
    logCaCO3 = log(CaCO3),
    year = date %>%
      as.character() %>%
      substr(start = 1, stop = 4) %>%
      as.numeric()
  )

fractions <- c("clay", "silt", "fine_sand", "coarse_sand", "logSOC", "logCaCO3")

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)

bounds_lower <- rep(0, 6)
bounds_upper <- rep(100, 6)


# 7: Make training data

folds <- bind_rows(
  dsc_folds,
  SEGES_folds,
  SINKS_folds,
  profiles_folds,
  forest_folds
)

names(folds) <- "fold"

extr <- bind_rows(
  dsc_extr,
  SEGES_extr,
  SINKS_extr,
  profiles_extr,
  forests_extr
)

obs <- cbind(obs_data, extr, folds) %>%
  filter(!is.na(UTMX) & !is.na(UTMY))

obs %<>%
  rownames_to_column() %>%
  mutate(ID_new = rowname, .before = everything()) %>%
  select(-rowname)

write.table(
  obs,
  paste0(dir_results, "observations_texture.csv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = ";"
)

obs_top <- obs %>%
  filter(
    upper < 25,
    is.finite(fold)
  )

obs_prf <- obs %>%
  filter(
    db == "Profile database",
    is.finite(fold)
  )

# Make new ID

obs_top_v <- obs_top %>% vect(geom = c("UTMX", "UTMY"))

library(viridisLite)

tiff(
  paste0(dir_results, "/obs_map_test", testn, ".tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

plot(
  obs_top_v, "clay",
  breaks = 5, breakby = "cases", col = cividis(5),
  cex = 0.2
)

try(dev.off())
try(dev.off())

plot(
  obs_top_v, "clay",
  breaks = 5, breakby = "cases", col = cividis(5),
  cex = 0.4
)

# 8: Set up models

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  dplyr::select(., name) %>%
  unlist() %>%
  unname()


# Template for custom eval
# evalerror <- function(preds, dtrain) {
#   labels <- getinfo(dtrain, "label")
#   err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
#   return(list(metric = "error", value = err))
# }

# Weighted RMSE
get_RMSEw <- function(d, w) {
  sqe <- w * (d[, 1] - d[, 2])^2
  msqe <- sum(sqe) / sum(w)
  out <- sqrt(msqe)
  return(out)
}

# Weighted R^2
get_R2w <- function(d, w) {
  require(boot)
  out <- boot::corr(d[, 1:2], w)^2
  return(out)
}

# Weighted summary function
WeightedSummary <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  out <- numeric()
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw", "R2w")
  return(out)
}

# Weighted summary function with log transformation
WeightedSummary_log <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  out <- numeric()
  data[, 1:2] <- log(data[, 1:2])
  data <- data[is.finite(rowSums(data)), ]
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw_log", "R2w_log")
  return(out)
}

# Weighted summary function with square root transformation
WeightedSummary_sqrt <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  out <- numeric()
  data[, 1:2] <- sqrt(data[, 1:2])
  data <- data[is.finite(rowSums(data)), ]
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw_sqrt", "R2w_sqrt")
  return(out)
}

metrics <- rep("RMSEw", length(fractions))
metrics[fractions == "SOC"] <- "RMSEw_log"
metrics[fractions == "CaCO3"] <- "RMSEw_sqrt"

# Function to calculate point densities

qnorm(seq(0.55, 0.95, 0.1), 0, 1)

get_dens <- function(datxy, sig) {
  dens_out <- ppp(
    datxy$UTMX,
    datxy$UTMY,
    c(441000, 894000),
    c(6049000, 6403000)
  ) %>%
    density(
      sigma = sig,
      at = "points",
      leaveoneout = FALSE
    )

  attributes(dens_out) <- NULL

  return(dens_out)
}

# Tuning grid

tgrid <- expand.grid(
  nrounds = 10,
  eta = 0.3,
  max_depth = 6,
  min_child_weight = 1,
  gamma = 0,
  colsample_bytree = 0.75,
  subsample = 0.75
)

eta_test <- seq(0.1, 1, 0.1)
max_depth_test <- seq(1, 30, 3)
min_child_weight_test <- c(1, 2, 4, 8, 16, 32, 64)
gamma_test <- seq(0, 0.6, 0.1)
colsample_bytree_test <- seq(0.1, 1, 0.1)
subsample_test <- seq(0.1, 1, 0.1)

objectives <- c(rep("reg:squarederror", 4), rep("reg:tweedie", 2))

trees_per_round <- 10

# Identify OGCs and make a list with the numbers of OGCs to be tested in the
# models

ogcs_names <- extr %>%
  names() %>%
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

# Bayesian optimization

library(ParBayesianOptimization)

bounds <- list(
  eta = c(0.1, 1),
  max_depth = c(1L, 50L),
  min_child_weight_sqrt = c(1, sqrt(64)),
  gamma_sqrt = c(0, sqrt(30)),
  colsample_bytree = c(0.1, 1),
  subsample = c(0.1, 1),
  colsample_bylevel = c(0.1, 1),
  ogcs_index = c(1L, 7L),
  total_imp = c(0.5, 1)
)

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
  cov_i_filtered <- cov_i_ranked %>%
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
      predictionBounds = c(bounds_lower_i, bounds_upper_i),
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

xgb_opt_stepwise <- FALSE

# Small random sample for testing
# Remember to include full dataset in the final model
n <- 1000

use_all_points <- TRUE
# use_all_points <- FALSE

extra_tuning_xgb <- TRUE
# extra_tuning_xgb <- FALSE

# 9: Train models

n_ogcs_models <- numeric()
total_imp_models <- numeric()

weights_objects <- list()

models_tr_summaries <- list()

models_scoreresults <- list()
models_bestscores <- list()

# Covariate selection:
# Step 1: Decide the optimal number of OGCs
# Step 2: Drop unimportant covariates

# xgb optimization:
# Step 1: Adjust learning rate
# Step 2: Fit max_depth and min_child_weight
# Step 3: Tune gamma
# Step 4: Adjust subsampling
# Step 5: Increase nrounds, readjust learning rate

# Test OGC after model tuning?

models_predictions <- matrix(
  numeric(),
  nrow = nrow(obs),
  ncol = length(fractions)
  )

models_weights <- matrix(
  numeric(),
  nrow = nrow(obs),
  ncol = length(fractions)
)

models_indices <- matrix(
  numeric(),
  nrow = nrow(obs),
  ncol = length(fractions)
)

colnames(models_predictions) <- fractions
colnames(models_weights) <- fractions

models <- list()

for (i in 1:length(fractions))
{
  frac <- fractions[i]

  print(frac)
  
  models_tr_summaries[[i]] <- list()
  tr_step <- 1

  if (metrics[i] == "RMSEw_log") {
    sumfun <- WeightedSummary_log
  } else {
    if (metrics[i] == "RMSEw_sqrt") {
      sumfun <- WeightedSummary_sqrt
    } else {
      sumfun <- WeightedSummary
    }
  }
  
  trdat <- obs %>%
    filter(
      is.finite(.data[[frac]]),
      !is.na(UTMX),
      !is.na(UTMY),
      lower > 0,
      upper < 200,
      !is.na(dhm2015_terraen_10m)
    )

  # Weighting by depth intervals
  print("Calculating weights")
  w_interval <- 10
  w_increment <- 1
  w_startdepth <- 0
  w_maxdepth <- 200
  w_depths <- seq(w_startdepth, w_maxdepth, w_increment)
  w_iterations <- length(w_depths)

  w_mat <- matrix(numeric(), nrow = nrow(trdat), ncol = w_iterations)
  cm_mat <- matrix(numeric(), nrow = nrow(trdat), ncol = w_iterations)

  for (j in 1:w_iterations)
  {
    upper_j <- w_depths[j] - w_interval
    lower_j <- upper_j + w_interval * 2

    trdat_ind <- trdat$lower > upper_j & trdat$upper < lower_j
    trdat_ind[is.na(trdat_ind)] <- FALSE

    trdat_j <- trdat[trdat_ind, ]
    
    trdat_j %<>%
      mutate(
        thickness = lower - upper,
        upper_int = case_when(
          upper > upper_j ~ upper,
          .default = upper_j
        ),
        lower_int = case_when(
          lower < lower_j ~ lower,
          .default = lower_j
        ),
        cm_int = case_when(
          thickness == 0 ~ 1,
          .default = lower_int - upper_int
        )
      )
    
    cm_mat[trdat_ind, j] <- trdat_j$cm_int

    # # Sigma equal to the radius of a circle with an equal area per sample
    # sigma_j <- sqrt(43107 / (nrow(trdat_j) * pi)) * 1000

    # Use the expected mean density as a baseline
    # Do this calculation for the entire area, as a specific calculation for
    # wetlands will give too much weight to these areas.
    mean_dens_j <- nrow(trdat_j) / (43107 * 10^6)

    # For SOC:
    # Separate densities for wetlands and uplands
    if (frac == "SOC") {
      areas <- c(39807, 3299)
      
      w_j <- numeric(nrow(trdat_j))

      for (k in 0:1) {
        trdat_j_wl_ind <- trdat_j$cwl_10m_crisp == k
        
        trdat_j_wl_ind %<>% { ifelse(is.na(.), FALSE, .) }

        trdat_jk <- trdat_j[trdat_j_wl_ind, ]
        
        # Sigma equal to the radius of a circle with an equal area per sample
        sigma_jk <- sqrt(areas[k + 1] / (nrow(trdat_jk) * pi)) * 1000

        dens_jk <- get_dens(trdat_jk, sigma_jk)
        
        w_j[trdat_j_wl_ind] <- mean_dens_j / dens_jk
      }
    } else {
      # Sigma equal to the radius of a circle with an equal area per sample
      sigma_j <- sqrt(43107 / (nrow(trdat_j) * pi)) * 1000
      
      dens_j <- get_dens(trdat_j, sigma_j)
      
      w_j <- mean_dens_j / dens_j
    }

    w_j[w_j > 1] <- 1
    w_mat[trdat_ind, j] <- w_j
  }

  cm_mat[is.na(cm_mat)] <- 0
  
  cm_sums <- apply(cm_mat, 1, sum)
  
  w_cm_mat <- w_mat*cm_mat
  
  w_cm_sums <- apply(
    w_cm_mat,
    1,
    function(x) {
      out <- sum(x, na.rm = TRUE)
      return(out)
    }
  )
  
  w_depth <- w_cm_sums / cm_sums
  w_depth[!is.finite(w_depth)] <- 0
  
  # Using the year as a covariate causes the model to identify the SINKS points
  # based only on the sampling year, as the campaign took place over only two
  # years, which were also underrepresented in the remaining training data.
  if (frac == "SOC") {
    w_year <- 0.99^(max(trdat$year, na.rm = TRUE) - trdat$year)
    w_year %<>% { ifelse(is.na(.), 0, .) }
    w_depth %<>% `*`(w_year)
  }
  
  weights_objects[[i]] <- list()
  weights_objects[[i]]$cm_mat <- cm_mat
  weights_objects[[i]]$cm_sums <- cm_sums
  weights_objects[[i]]$w_cm_mat <- w_cm_mat
  weights_objects[[i]]$w_cm_sums <- w_cm_sums
  weights_objects[[i]]$w_depth <- w_depth

  trdat$w <- w_depth
  
  trdat_w_indices <- which(obs$ID_new %in% trdat$ID_new)
  
  weights_objects[[i]]$indices <- trdat_w_indices
  
  models_weights[trdat_w_indices, i] <- w_depth
  
  # Three folds (placeholder)
  trdat %<>% mutate(
    fold = ceiling(fold / 3)
  )
  
  trdat %<>% filter(fold < 4)
  
  if (!use_all_points) {
    set.seed(1)
    trdat %<>% sample_n(n)
  }
  
  trdat_indices <- which(obs$ID_new %in% trdat$ID_new)
  
  models_indices[, i] <- obs$ID_new %in% trdat$ID_new
  
  holdout_i <- obs[-trdat_indices, ]
  holdout_indices <- which(obs$ID_new %in% holdout_i$ID_new)

  # List of folds

  folds_i <- lapply(
    unique(trdat$fold),
    function(x) {
      out <- trdat %>%
        mutate(
          is_j = fold != x,
          rnum = row_number(),
          ind_j = is_j * rnum
        ) %>%
        filter(ind_j != 0) %>%
        dplyr::select(., ind_j) %>%
        unlist() %>%
        unname()
    }
  )

  showConnections()

  # Add depth boundaries and methods as covariates
  cov_c_i <- cov_selected %>%
    c("upper", "lower")
  if (i %in% 1:4) {
    cov_c_i <- cov_selected %>%
      c("upper", "lower") %>%
      c(., "SOM_removed")
  } 
  
  # Identify covariates that are not OGCs
  
  covs_not_ogc <- grep('ogc_pi', cov_c_i, value = TRUE, invert = TRUE)
  
  cov_p_i <- covs_not_ogc %>% paste0(collapse = " + ")
  
  formula_i <- paste0(frac, " ~ ", cov_p_i) %>%
    as.formula()
  
  # xgboost optimization
  # 1: Fit learning rate (eta)
  print("Step 1: Fit learning rate (eta)")
  
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

  models[[i]] <- caret::train(
    form = formula_i,
    data = trdat,
    method = "xgbTree",
    na.action = na.pass,
    tuneGrid = expand.grid(
      nrounds = tgrid$nrounds,
      eta = eta_test, # NB
      max_depth = tgrid$max_depth,
      min_child_weight = tgrid$min_child_weight,
      gamma = tgrid$gamma,
      colsample_bytree = tgrid$colsample_bytree,
      subsample = tgrid$subsample
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
    objective = objectives[i],
    colsample_bylevel = 0.75,
    nthread = 1
  )
  
  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)
  
  models_tr_summaries[[i]][[tr_step]] <- models[[i]]$results
  print(models_tr_summaries[[i]][[tr_step]])
  tr_step %<>% `+`(1)
  
  # Scaled cumulative covariate importance
  cov_i_ranked <- varImp(models[[i]])$importance %>%
    rownames_to_column() %>%
    mutate(
      scaled = Overall/sum(Overall),
      cumul = cumsum(scaled)
    )
  
  if (extra_tuning_xgb & xgb_opt_stepwise) {
    # CS Step 1: Drop unimportant covariates
    cov_c_i <- cov_i_ranked %>%
      filter(cumul < 0.99) %>%
      .$rowname
    
    cov_p_i <- cov_c_i %>% paste0(collapse = " + ")
    
    formula_i <- paste0(frac, " ~ ", cov_p_i) %>%
      as.formula()
    
    # xgb opt Step 2: Fit max_depth and min_child_weight
    print("Step 2: Fit max_depth and min_child_weight")

    set.seed(1)

    model2 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = max_depth_test, # NB
        min_child_weight = min_child_weight_test, # NB
        gamma = models[[i]]$bestTune$gamma,
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metrics[i],
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )
    
    models[[i]] <- model2
    
    models_tr_summaries[[i]][[tr_step]] <- models[[i]]$results
    print(models_tr_summaries[[i]][[tr_step]])
    tr_step %<>% `+`(1)

    # xgb opt Step 3: Tune gamma
    print("Step 3: Tune gamma")

    set.seed(1)

    model3 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = models[[i]]$bestTune$max_depth,
        min_child_weight = models[[i]]$bestTune$min_child_weight,
        gamma = gamma_test, # NB
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metrics[i],
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )

    models[[i]] <- model3
    
    models_tr_summaries[[i]][[tr_step]] <- models[[i]]$results
    print(models_tr_summaries[[i]][[tr_step]])
    tr_step %<>% `+`(1)
    
    # xgb opt Step 4: Adjust subsampling
    print("Step 4: Adjust subsampling")
    
    set.seed(1)
    
    model4 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = models[[i]]$bestTune$max_depth,
        min_child_weight = models[[i]]$bestTune$min_child_weight,
        gamma = models[[i]]$bestTune$gamma,
        colsample_bytree = colsample_bytree_test,
        subsample = subsample_test
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metrics[i],
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )
    
    models[[i]] <- model4
    
    models_tr_summaries[[i]][[tr_step]] <- models[[i]]$results
    print(models_tr_summaries[[i]][[tr_step]])
    tr_step %<>% `+`(1)
    
    # CS Step 2: Decide the optimal number of OGCs
    cov_c_i <- varImp(models[[i]])$importance %>%
      rownames_to_column() %>%
      .$rowname
    
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
    
    ogcs_names_list %<>% lapply(., function(x) {c(cov_c_i, x)})
    ogcs_names_list[[length(ogcs_names_list) + 1]] <- cov_c_i
    n_ogcs_v %<>% c(., 0)
    
    print("Testing OGCs")
    
    models_ogc_test <- ogcs_names_list %>%
      lapply(
        .,
        function(x) {
          cov_p_i <- x %>% paste0(collapse = " + ")
          
          formula_i <- paste0(frac, " ~ ", cov_p_i) %>%
            as.formula()
          
          set.seed(1)
          
          out <- caret::train(
            form = formula_i,
            data = trdat,
            method = "xgbTree",
            na.action = na.pass,
            tuneGrid = expand.grid(
              nrounds = models[[i]]$bestTune$nrounds,
              eta = models[[i]]$bestTune$eta,
              max_depth = models[[i]]$bestTune$max_depth,
              min_child_weight = models[[i]]$bestTune$min_child_weight,
              gamma = models[[i]]$bestTune$gamma,
              colsample_bytree = models[[i]]$bestTune$colsample_bytree,
              subsample = models[[i]]$bestTune$subsample
            ),
            trControl = trainControl(
              index = folds_i,
              savePredictions = "final",
              predictionBounds = c(bounds_lower[i], bounds_upper[i]),
              summaryFunction = sumfun,
              allowParallel = FALSE
            ),
            metric = metrics[i],
            maximize = FALSE,
            weights = trdat$w,
            num_parallel_tree = trees_per_round,
            objective = objectives[i]
          )
          return(out)
        }
      )
    
    ogc_results <- models_ogc_test %>% lapply(
      ., function(x) x$results %>% select(any_of(metrics[i])) %>% min()
    ) %>%
      unlist()
    
    which_ogc_ind <- which.min(ogc_results)
    
    ogc_df <- data.frame(
      fraction = frac,
      n_ogcs = n_ogcs_v,
      acc = ogc_results,
      metric = metrics[i]
    )
    
    write.table(
      ogc_df,
      paste0(dir_results, "ogc_acc_", frac, ".csv"),
      row.names = FALSE,
      col.names = TRUE,
      sep = ";"
    )
    
    models_tr_summaries[[i]][[tr_step]] <- ogc_df
    print(models_tr_summaries[[i]][[tr_step]])
    tr_step %<>% `+`(1)
    
    n_ogcs_models[i] <- n_ogcs_v[which_ogc_ind]
    
    models[[i]] <- models_ogc_test[[which_ogc_ind]]
    
    cov_c_i <- varImp(models_ogc_test[[which_ogc_ind]])$importance %>%
      rownames_to_column() %>%
      .$rowname

    cov_p_i <- cov_c_i %>% paste0(collapse = " + ")

    formula_i <- paste0(frac, " ~ ", cov_p_i) %>%
      as.formula()
    
    # xgb opt Step 5: Increase nrounds, readjust learning rate
    print("Step 5")
    
    eta_test_final <- model4$bestTune$eta %>%
      log() %>%
      seq(., . + log(0.01), length.out = 9) %>%
      exp() %>%
      round(3)
    
    set.seed(1)
    
    model5 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds*10,
        eta = eta_test_final, # NB
        max_depth = models[[i]]$bestTune$max_depth,
        min_child_weight = models[[i]]$bestTune$min_child_weight,
        gamma = models[[i]]$bestTune$gamma,
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = sumfun,
        allowParallel = FALSE
      ),
      metric = metrics[i],
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )
    
    models[[i]] <- model5
    
    models_tr_summaries[[i]][[tr_step]] <- models[[i]]$results
    print(models_tr_summaries[[i]][[tr_step]])
    tr_step %<>% `+`(1)
  }
  
  # Bayes optimization
  if (extra_tuning_xgb & !xgb_opt_stepwise) {
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
      bounds = bounds,
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
      objective = objectives[i],
      colsample_bylevel = best_pars$colsample_bylevel,
      nthread = 1
    )
    
    stopCluster(cl)
    foreach::registerDoSEQ()
    rm(cl)
    
    models[[i]] <- model_final
  }
  
  print(models[[i]])
  
  models_predictions[trdat_indices, i] <- models[[i]]$pred %>%
    arrange(rowIndex) %>%
    distinct(rowIndex, .keep_all = TRUE) %>%
    dplyr::select(., pred) %>%
    unlist() %>%
    unname()
  
  if (i %in% 1:4) {
    n_const_i <- 3
  } else {
    n_const_i <- 2
  }
  
  models_predictions[holdout_indices, i] <- predict_passna(
    models[[i]],
    holdout_i,
    n_const = n_const_i
    )
  
  saveRDS(
    models[[i]],
    paste0(dir_results, "/model_", frac, ".rds")
  )
}

# End of model training

names(weights_objects) <- fractions

saveRDS(
  weights_objects,
  paste0(dir_results, "/weights_objects.rds")
)

saveRDS(
  models_tr_summaries,
  paste0(dir_results, "/models_tr_summaries.rds")
)

saveRDS(
  models_scoreresults,
  paste0(dir_results, "/models_scoreresults.rds")
)

saveRDS(
  models_bestscores,
  paste0(dir_results, "/models_bestscores.rds")
)

write.table(
  models_weights,
  file = paste0(dir_results, "/models_weights.csv"),
  sep = ";",
  row.names = FALSE
)

write.table(
  models_indices,
  file = paste0(dir_results, "/models_indices.csv"),
  sep = ";",
  row.names = FALSE
)

models_loaded <- lapply(
  1:6,
  function(x) {
    out <- fractions[x] %>%
      paste0(dir_results, "/model_", ., ".rds") %>%
      readRDS()
    return(out)
  }
)

models <- models_loaded

names(models) <- fractions

# Model summary

models_sum <- lapply(models, function(x) {
  results <- x$results
  if ("RMSEw" %in% names(results)) {
    out <- results %>% filter(RMSEw == min(RMSEw))
  } else {
    if ("RMSEw_sqrt" %in% names(results)) {
      out <- results %>% filter(RMSEw_sqrt == min(RMSEw_sqrt))
    } else {
      out <- results %>%
        filter(RMSEw_log == min(RMSEw_log))
    }
  }
  return(out)
}) %>%
  bind_rows() %>%
  mutate(
    Fraction = fractions,
    .before = 1
  ) %>% mutate(
    ogcs = n_ogcs_models,
    total_imp = total_imp_models,
    .after = 1
  ) %T>%
  write.table(
    file = paste0(dir_results, "/models_sum.csv"),
    sep = ";",
    row.names = FALSE
  )

models_sum

# Covariate importance

imp_all <- models %>%
  seq_along() %>%
  lapply(
    function(x) {
      out <- varImp(models[[x]])$importance %>%
        as.data.frame() %>%
        rownames_to_column(var = "covariate") %>%
        mutate(fraction = names(models)[x])
      return(out)
    }
  ) %>%
  bind_rows() %>%
  pivot_wider(
    id_cols = covariate,
    names_from = fraction,
    values_from = Overall
  ) %>%
  replace(is.na(.), 0) %>%
  rowwise() %>%
  mutate(mean_imp = mean(c_across(-covariate))) %>%
  arrange(-mean_imp) %T>%
  write.table(
    file = paste0(dir_results, "/var_imp.csv"),
    sep = ";",
    row.names = FALSE
  )

imp_all

# Inspect models

get_acc <- function(x2, i2) {
  df <- x2$pred %>%
    arrange(rowIndex) %>%
    distinct(rowIndex, .keep_all = TRUE) %>%
    dplyr::select(., c(pred, obs, weights)) %>%
    mutate(
      pred = ifelse(pred < 0, 0, pred)
    )

  # if (i2 > 4) df %<>% exp

  df %<>% bind_cols(x2$trainingData)

  r2_all <- df %$% get_R2w(cbind(pred, obs), weights)

  r2_bare <- df %>%
    filter(!is.na(s2_geomedian_b2)) %$%
    get_R2w(cbind(pred, obs), weights)

  r2_covered <- df %>%
    filter(is.na(s2_geomedian_b2)) %$%
    get_R2w(cbind(pred, obs), weights)

  rmse_all <- df %$% get_RMSEw(cbind(pred, obs), weights)

  rmse_bare <- df %>%
    filter(!is.na(s2_geomedian_b2)) %$%
    get_RMSEw(cbind(pred, obs), weights)

  rmse_covered <- df %>%
    filter(is.na(s2_geomedian_b2)) %$%
    get_RMSEw(cbind(pred, obs), weights)

  out <- data.frame(
    r2_all,
    r2_bare,
    r2_covered,
    rmse_all,
    rmse_bare,
    rmse_covered
  )

  return(out)
}

acc_all <- foreach(i = 1:6, .combine = rbind) %do%
  get_acc(models[[i]], i)

acc_all %<>% mutate(fraction = fraction_names, .before = 1)

write.table(
  acc_all,
  paste0(dir_results, "/acc_all_test", testn, ".csv"),
  sep = ";",
  row.names = FALSE
)

# Standardize and transform predictions

models_predictions

models_predictions %<>%
  apply(
    1, function(x) {
      sum_min <- sum(x[1:4])
      
      if (!is.na(sum_min)) {
        x[1:4] %<>% `/`(sum_min) %>% `*`(100)
      }
      
      x %<>% matrix(ncol = length(x)) %>% as.data.frame()
      
      return(x)
    },
    simplify = FALSE
  ) %>%
  bind_rows()

colnames(models_predictions) <- fractions

write.table(
  models_predictions,
  file = paste0(dir_results, "/models_predictions_raw.csv"),
  sep = ";",
  row.names = FALSE
)

source("f_classify_soil_JB.R")

JB_predicted <- classify_soil_JB(
  clay = models_predictions[, 1],
  silt = models_predictions[, 2],
  sand_f = models_predictions[, 3],
  SOM = models_predictions[, 5] / 0.568,
  CaCO3 = models_predictions[, 6]
) %>%
  factor(
    levels = 1:12,
    labels = paste0("JB", 1:12)
    )

JB_observed <- classify_soil_JB(
  clay = obs$clay,
  silt = obs$silt,
  sand_f = obs$fine_sand,
  SOM = obs$SOC / 0.568,
  CaCO3 = obs$CaCO3
) %>%
  factor(
    levels = 1:12,
    labels = paste0("JB", 1:12)
  )

mean_weights <- apply(models_weights, 1, function(x) { mean(x, na.rm = TRUE) })

JB_df <- data.frame(
  ID_new = obs$ID_new,
  observed = JB_observed,
  predicted = JB_predicted,
  upper = obs$upper,
  lower = obs$lower,
  imputed = obs$imputed,
  w = NA,
  fold = obs$fold
) %>%
  filter(imputed == FALSE) %>%
  select(-imputed)

JB_df %>%
  filter(fold != 10) %>%
  select(c(predicted, observed)) %>%
  na.omit() %>%
  table() %>%
  confusionMatrix()


# Analyse predictions and observations

getpred <- function(
    x2, i2, axmaxmin = TRUE, drop_imputed = TRUE, drop_inf = TRUE
) {
  
  df <- obs %>% select(any_of(fractions[i2]))
  colnames(df) <- "obs"
  df$pred = models_predictions[, i2]
  
  if (i2 == 5) df %<>% log()
  if (i2 == 6) df %<>% sqrt()
  
  depths_x <- obs %>%
    select(upper, lower, fold, ID_new, imputed)
  df %<>% bind_cols(depths_x)
  df %<>% mutate(
    fraction = fractions[i2],
    w = models_weights[, i2]
  ) %>%
    filter(
      !is.na(obs),
      !is.na(w)
      )
  
  if (drop_inf) {
    df %<>% filter(is.finite(obs))
  }
  
  if (drop_imputed) {
    df %<>% filter(imputed == FALSE)
  }
  
  if (axmaxmin) {
    if (i2 == 5) {
      df %<>%
        mutate(
          ax_min = quantile(obs, 0.01)
        )
    } else {
      df %<>% mutate(
        ax_min = 0
      )
    }
    df %<>% mutate(
      ax_max = quantile(obs, 0.99)
    )
    df %<>%
      filter(
        obs < ax_max,
        pred < ax_max,
        obs >= ax_min,
        pred >= ax_min
      )
  }
  
  return(df)
}

allpred <- foreach(i = 1:6, .combine = rbind) %do%
  getpred(models[[i]], i)

allpred$fraction %<>% factor(levels = fractions)

levels(allpred$fraction) <- c(
  expression("Clay"~"(%)"),
  expression("Silt"~"(%)"),
  expression("Fine"~"sand"~"(%)"),
  expression("Coarse"~"sand"~"(%)"),
  expression("log[SOC"~"(%)]"),
  expression(sqrt("CaCO"[3]~"(%)"))
)

my_breaks_dens <- c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1)

# mycolors1 <- rgb(
#   seq(1, 0, length.out = 10),
#   seq(1, 1, length.out = 10),
#   seq(1, 1, length.out = 10)
#   )
mycolors2 <- rgb(
  seq(0, 0, length.out = 20),
  seq(1, 0, length.out = 20),
  seq(1, 1, length.out = 20)
)
mycolors3 <- rgb(
  seq(0, 0, length.out = 21),
  seq(0, 0, length.out = 21),
  seq(1, 0, length.out = 21)
)

mycolorgradient <- c(
  # mycolors1,
  mycolors2[-1],
  mycolors3[-1])

breaks <- c(0, 30, 60, 100, 200)

tr_partitions <- c("CV", "Holdout")

for (k in 1:2) {
  if (k == 1) {
    allpred_k <- allpred %>%
      filter(fold < 10)
  } else {
    allpred_k <- allpred %>%
      filter(fold == 10)
  }
  
  for (i in 1:(length(breaks) - 1)) {
    allpred_i <- allpred_k %>%
      filter(
        upper < breaks[i + 1],
        lower > breaks[i]
      )
    
    if (nrow(allpred_i) > 0) {
      myplot <- allpred_i %>%
        ggplot(aes(x = obs, y = pred, weight = w)) +
        ggtitle(
          paste0(
            "Accuracy at ", breaks[i], " - ", breaks[i + 1], " cm depth, ",
            tr_partitions[k]
          )
        ) + 
        geom_hex(
          aes(
            alpha = after_stat(ndensity),
            fill = after_stat(ndensity)
          ),
          bins = 30
        ) +
        scale_fill_gradientn(
          colours = mycolorgradient,
          aesthetics = "fill",
          breaks = my_breaks_dens,
          limits = c(0, 1),
          na.value = 0,
          trans = "sqrt"
        ) +
        guides(
          fill = "legend"
          , alpha = "legend"
        ) +
        scale_alpha_continuous(
          limits= c(0, 1),
          range = c(0.05, 1),
          breaks = my_breaks_dens,
          na.value = 0,
          trans = "sqrt"
        ) +
        facet_wrap(~fraction, nrow = 2, scales = "free", labeller = label_parsed) +
        theme_bw() +
        theme(aspect.ratio = 1) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        geom_abline(col = "red") +
        geom_blank(aes(y = ax_max)) +
        geom_blank(aes(x = ax_max)) +
        geom_blank(aes(y = ax_min)) +
        geom_blank(aes(x = ax_min)) +
        xlab("Observation") +
        ylab("Prediction")
      
      tiff(
        paste0(dir_results, "/", tr_partitions[k], "_accuracy_", LETTERS[i],
               "_test", testn, ".tiff"),
        width = 15,
        height = 10,
        units = "cm",
        res = 300
      )
      
      print(myplot)
      
      try(dev.off())
    }
  }
}

tiff(
  paste0(dir_results, "/accuracy_all_test", testn, ".tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

allpred %>%
  ggplot(aes(x = obs, y = pred, weight = w)) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
        ),
    bins = 30
    ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits= c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
    ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  facet_wrap(~fraction, nrow = 2, scales = "free", labeller = label_parsed) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(col = "red") +
  geom_blank(aes(y = ax_max)) +
  geom_blank(aes(x = ax_max)) +
  geom_blank(aes(y = ax_min)) +
  geom_blank(aes(x = ax_min)) +
  xlab("Observation") +
  ylab("Prediction")

try(dev.off())


# Accuracy vs depth

d_out <- list()

# Also calculate statistics by depth for the observations (to do)

allpred_depth <- foreach(i = 1:6) %do%
  getpred(
    models[[i]], i, axmaxmin = FALSE, drop_imputed = FALSE, drop_inf = FALSE
  )

for (i in 1:length(models)) {
  depth_weights <- weights_objects[[i]]$w_cm_mat / weights_objects[[i]]$w_cm_sums
  
  mdata <- allpred_depth[[i]]
  
  d_out[[i]] <- apply(
    depth_weights, 2, 
    function(x) {
      ddat <- mdata %>%
        mutate(
          w_div = x
        ) %>%
        filter(
          is.finite(w_div),
          is.finite(obs),
          is.finite(pred)
        )
      
      if (i == 5) {
        ddat %<>% filter(imputed == FALSE)
      }
      
      RMSEws <- numeric()
      R2ws <- numeric()
      weights_ks <- numeric()
      
      for (k in 1:2) {
        if (k == 1) {
          ddat_k <- ddat %>%
            filter(fold < 10)
        } else {
          ddat_k <- ddat %>%
            filter(fold == 10)
        }
        
        if (nrow(ddat_k) > 2) {
          RMSEws[k] <- get_RMSEw(
            select(ddat_k, c(pred, obs)),
            ddat_k$w_div
          )
          R2ws[k] <- get_R2w(
            select(ddat_k, c(pred, obs)),
            ddat_k$w_div
          )
          weights_ks[k] <- sum(ddat_k$w_div)
        } else {
          
          RMSEws[k] <- NA
          R2ws[k] <- NA
          weights_ks[k] <- NA
        }
      }
      
      out <- data.frame(
        RMSEw = RMSEws,
        R2w = R2ws,
        Weights = weights_ks,
        Dataset = tr_partitions
      )
      
      return(out)
    }
  ) %>% bind_rows()
  
  d_out[[i]]$Fraction <- fractions[i]
  d_out[[i]]$Depth = rep(w_depths, each = 2)
}

d_out %<>%
  bind_rows() %>%
  mutate(
    Fraction = factor(
      Fraction,
      levels = fractions,
      labels = c(
        expression("Clay"~"(%)"),
        expression("Silt"~"(%)"),
        expression("Fine"~"sand"~"(%)"),
        expression("Coarse"~"sand"~"(%)"),
        expression("log[SOC"~"(%)]"),
        expression(sqrt("CaCO"[3]~"(%)"))
      )
    ),
    Dataset = as.factor(Dataset)
  )

myBreaks <- function(x) {
  breaks <- c(
    0,
    (max(x)*1/2) %>% log2() %>% round(., digits = 0) %>% (2)^.
    )
  names(breaks) <- attr(breaks,"labels")
  return(breaks)
}

tiff(
  paste0(dir_results, "/depth_RMSEw_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

d_out %>%
  ggplot(aes(x = RMSEw, y = Depth, color = Dataset)) +
  facet_wrap(~Fraction, nrow = 1, scales = "free_x", labeller = label_parsed) +
  geom_path() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(
    guide = guide_axis(check.overlap = TRUE),
    breaks = myBreaks,
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.2), 
                       add = c(0, 0))
    ) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 6)) +
  labs(
    x = bquote(RMSE[w]),
    y = "Depth (cm)"
  )


try(dev.off())
try(dev.off())

tiff(
  paste0(dir_results, "/depth_R2w_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

d_out %>%
  ggplot(aes(x = R2w, y = Depth, color = Dataset)) +
  facet_wrap(~Fraction, nrow = 1, labeller = label_parsed) +
  geom_path() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(
    guide = guide_axis(check.overlap = TRUE),
    breaks = c(0, 0.3, 0.6),
    expand = expansion(mult = c(0, 0), 
                 add = c(0.01, 0.02))
    ) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 6)) +
  labs(
    x = bquote({R^2}[w]),
    y = "Depth (cm)"
  )

try(dev.off())
try(dev.off())

tiff(
  paste0(dir_results, "/depth_weights_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

d_out %>%
  ggplot(aes(x = Weights, y = Depth, color = Dataset)) +
  facet_wrap(~Fraction, nrow = 1, labeller = label_parsed) +
  geom_path() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(
    guide = guide_axis(check.overlap = TRUE),
    breaks = c(0, 2000),
    expand = expansion(mult = c(0.01, 0.1), 
                       add = c(0.01, 0.1))
    ) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 6)) +
  labs(
    x = "Total weights",
    y = "Depth (cm)"
  )

try(dev.off())
try(dev.off())


# Covariate importance

library(colorRamps)
library(rcartocolor) # for colorblind palette

mycolors <- carto_pal(12, "Safe") %>% sort()

library(TSP)

l <- list()

ntop <- 20

for (i in 1:length(models))
{
  l[[i]] <- varImp(models[[i]])$importance %>%
    as_tibble(rownames = "covariate") %>%
    drop_na() %>%
    arrange(-Overall) %>%
    slice_head(n = ntop) %>%
    mutate(target = fractions[i]) %>%
    mutate(rank = 1:ntop)
}

l %<>% bind_rows() %>%
  mutate(
    target = factor(
      target,
      levels = fractions
    )
  )

l_cat <- cov_cats %>%
  mutate(
    covariate = name,
    category = case_when(
      category == "basic" ~ scorpan,
      category == "WATEM" ~ "OR",
      category == "sentinel_composite" ~ "S2 time series",
      category == "bare_soil" ~ "Bare soil",
      .default = "Other"
    )
  )

l %<>%
  left_join(l_cat)

l %<>%
  ungroup() %>%
  arrange(target, Overall) %>%
  mutate(order = row_number())

l %<>% mutate(
  category = case_when(
    covariate == "upper" ~ "Depth",
    covariate == "lower" ~ "Depth",
    covariate == "year" ~ "Time",
    category == "N" ~ "Spatial position",
    category == "R" ~ "Topography",
    category == "C" ~ "Climate",
    category == "C " ~ "Climate",
    category == "P" ~ "Parent materials",
    category == "S" ~ "Soil",
    category == "SO" ~ "Soil and organisms",
    category == "CR" ~ "Climate and topography",
    category == "OR" ~ "Organisms and topography",
    category == "O" ~ "Organisms",
    category == "RP" ~ "Topography and parent materials",
    .default = category
  )
)

l$category %<>% as.factor()

catcolors <- l$category %>%
  levels() %>%
  length() %>%
  carto_pal(., "Safe")
names(catcolors) <- levels(l$category)
colScale <- scale_fill_manual(name = "category", values = catcolors)

# Plot covariate importance

tiff(
  paste0(dir_results, "/importance_test", testn, ".tiff"),
  width = 40,
  height = 20,
  units = "cm",
  res = 300
)

l %>%
  ggplot(aes(x = order, y = Overall, bg = category)) +
  geom_col() +
  facet_wrap(
    ~target,
    ncol = 3,
    scales = "free"
  ) +
  # xlim(1, ntop) +
  ylim(0, NA) +
  coord_flip() +
  scale_x_continuous(
    breaks = l$order,
    labels = l$covariate,
    expand = c(0, 0)
  ) +
  colScale

try(dev.off())

# 10: Make maps for the test area

dir_cov_10km <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

predfolder <- dir_results %>%
  paste0(., "/predictions_testarea/") %T>%
  dir.create()

source("f_predict_passna.R")

# Make the maps

breaks <- c(0, 30, 60, 100, 200)

uppers <- breaks %>% rev() %>% .[-1] %>% rev()
lowers <- breaks %>% .[-1]

map_spec <- expand_grid(
  fraction_i = 1:6,
  interval = 1:4
)

showConnections()

numCores <- 19

cl <- makeCluster(numCores)

clusterEvalQ(
  cl,
  {
    library(terra)
    library(caret)
    library(xgboost)
    library(magrittr)
    library(dplyr)
    library(tools)
  }
)

clusterExport(
  cl,
  c(
    "uppers",
    "lowers",
    "map_spec",
    "predfolder",
    "dir_cov_10km",
    "models",
    "cov_selected",
    "predict_passna",
    "dir_dat",
    "fractions"
  )
)

parSapplyLB(
  cl,
  1:nrow(map_spec),
  function(x) {
    tmpfolder <- paste0(dir_dat, "/Temp/")

    terraOptions(memfrac = 0.02, tempdir = tmpfolder)

    cov_10km <- dir_cov_10km %>%
      list.files(full.names = TRUE) %>%
      rast() %>%
      subset(cov_selected)

    outname <- predfolder %>%
      paste0(
        ., "/", fractions[map_spec$fraction_i[x]],
        "_depth", map_spec$interval[x],
        ".tif"
      )

    predict(
      cov_10km,
      models[[map_spec$fraction_i[x]]],
      fun = predict_passna,
      na.rm = FALSE,
      const = data.frame(
        SOM_removed = 1,
        # year = 2010,
        upper = uppers[map_spec$interval[x]],
        lower = lowers[map_spec$interval[x]]
      ),
      n_const = 3,
      n_digits = 1,
      filename = outname,
      overwrite = TRUE
    )

    return(NA)
  }
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)

maps_10_km <- list()

for (i in 1:length(fractions)) {
  maps_10_km[[i]] <- c(1:4) %>%
    paste0(
      predfolder, "/", fractions[i],
      "_depth", .,
      ".tif"
    ) %>%
    rast()
  names(maps_10_km[[i]]) <- paste0(
    fraction_names[i], ", ", uppers, " - ", lowers, " cm"
  )
}

# Standardize mineral sum to 100

maps_10_km_mineral_fin <- lapply(
  1:4, function(i) {
    mineral_raw <- c(
      maps_10_km[[1]][[i]],
      maps_10_km[[2]][[i]],
      maps_10_km[[3]][[i]],
      maps_10_km[[4]][[i]]
    )
    
    mineral_sum_r <- mineral_raw %>% sum() 
    
    mineral_final <- mineral_raw*100 / mineral_sum_r
    
    mineral_final %<>% round(., digits = 1)
    
    return(mineral_final)
  }
)

maps_10_km_mineral_fin_frac <- lapply(
  1:4, function(x) {
    out <- c(
      maps_10_km_mineral_fin[[1]][[x]],
      maps_10_km_mineral_fin[[2]][[x]],
      maps_10_km_mineral_fin[[3]][[x]],
      maps_10_km_mineral_fin[[4]][[x]]
    )
  }
)

for (i in 1:length(maps_10_km_mineral_fin_frac)) {
  maps_10_km[[i]] <- maps_10_km_mineral_fin_frac[[i]]
}

# Figures for 10 km maps

library(viridisLite)
library(tidyterra)

try(dev.off())
lapply(1:6, function(x) {
  fname <- paste0(dir_results, "/", fractions[x], "_10km_test", testn, ".tiff")

  myplot <- autoplot(maps_10_km[[x]]) +
    scale_fill_gradientn(colours = viridis(100), na.value = NA)

  tiff(
    fname,
    width = 16,
    height = 14,
    units = "cm",
    res = 300
  )

  print(myplot)

  try(dev.off())
  try(dev.off())
})

# Map and plot soil class (JBNR)

source("f_classify_soil_JB.R")

maps_10km_jb <- lapply(
  1:length(uppers),
  function(x) {
    maps_10_km_s2 <- c(
      maps_10_km[[1]][[x]],
      maps_10_km[[2]][[x]],
      maps_10_km[[3]][[x]],
      maps_10_km[[5]][[x]] / 0.568,
      maps_10_km[[6]][[x]]
    )

    names(maps_10_km_s2) <- c("clay", "silt", "sand_f", "SOM", "CaCO3")

    out <- lapp(maps_10_km_s2, classify_soil_JB)
    return(out)
  }
) %>%
  rast()

levels(maps_10km_jb) <- rep(
  list(
    data.frame(
      id = 1:12,
      Class = paste0("JB", 1:12)
    )
  ),
  nlyr(maps_10km_jb)
)

names(maps_10km_jb) <- paste0("JB class, ", uppers, " - ", lowers, " cm")

myrgb <- col2rgb(mycolors)
tsp <- as.TSP(dist(t(myrgb)))
set.seed(1)
sol <- solve_TSP(tsp, control = list(repetitions = 1e3))
ordered_cols <- mycolors[sol]

classes_in_maps <- values(maps_10km_jb) %>%
  unlist() %>%
  matrix(ncol = 1) %>%
  unique() %>%
  sort()

cols_in_maps <- ordered_cols[classes_in_maps]

plot_jb <- autoplot(maps_10km_jb) +
  scale_fill_discrete(type = cols_in_maps)

tiff(
  paste0(dir_results, "/JB_test_", testn, ".tiff"),
  width = 16,
  height = 14,
  units = "cm",
  res = 300
)

print(plot_jb)

try(dev.off())

# END