# 07: Train texture models

# xgboost
# Use all observations, including NA (OK)

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
library(spatstat)  # weights

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

# To do:
# Accuracy by depth
# Maps for depths (ok)
# Effects of som removal
# Pdp with depth
# Profile examples
# Adaptive kernel for point densities


# Test 1 - 8: Cubist
# Test 8: New covariates (chelsa, river valley bottoms, hillyness)
# Test 9: xgboost
# Test 10: Predicted row and column (poor accuracy)
# Test 11: Fever data, more 
testn <- 12
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
  paste0(., "/cov_categories_20230501.csv") %>%
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

# bounds_lower <- c(0, 0, 0, 0, NA, NA)
bounds_lower <- rep(0, 6)
# bounds_upper <- c(100, 100, 100, 100, log(100), log(100))
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
# Use all observations?

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
  obs_top_v, "clay", breaks = 5, breakby = "cases", col = cividis(5),
  cex = 0.2
)

try(dev.off())
try(dev.off())

plot(
  obs_top_v, "clay", breaks = 5, breakby = "cases", col = cividis(5),
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
get_RMSEw <- function(d, w)
{
  sqe <- w*(d[, 1] - d[, 2])^2
  msqe <- sum(sqe)/sum(w)
  out <- sqrt(msqe)
  return(out)
}

# Weighted R^2
get_R2w <- function(d, w)
{
  require(boot)
  out <- boot::corr(d[, 1:2], w)^2
  return(out)
}

# Weighted summary function
WeightedSummary <- function (
    data,
    lev = NULL,
    model = NULL,
    ...
) {
  out <- numeric()
  # Weighted RMSE
  get_RMSEw <- function(d, w) {
    sqe <- w*(d[, 1] - d[, 2])^2
    msqe <- sum(sqe)/sum(w)
    out <- sqrt(msqe)
    return(out)
  }
  # Weighted R^2
  get_R2w <- function(d, w)
  {
    require(boot)
    out <- boot::corr(d[, 1:2], w)^2
    return(out)
  }
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c('RMSEw', 'R2w')
  return(out)
}

# Tuning grid

tgrid <- expand.grid(
  nrounds = 100,
  eta = seq(0.1, 1, 0.1),
  max_depth = 6,
  min_child_weight = 1,
  gamma = 0,
  colsample_bytree = 0.5,
  subsample = 0.3
)

max_depth_test <- seq(1, 20, 3)
min_child_weight_test <- c(1, 2, 4, 8, 16, 32)
gamma_test <- seq(0, 0.5, 0.1)

objectives <- c(rep("reg:squarederror", 4), rep("reg:tweedie", 2))

trees_per_round <- 10


# Small random sample for testing
# Remember to include full dataset in the final model
n <- 1000

use_all_points <- TRUE
# use_all_points <- FALSE

# 9: Train models

extra_tuning_xgb <- TRUE

models <- list()

for (i in 1:length(fractions))
{
  frac <- fractions[i]

  print(frac)

  cov_c_i <- cov_selected %>%
    c("upper", "lower") %>%
    paste0(collapse = " + ")
  if (i %in% 1:4) {
    cov_c_i <- cov_selected %>%
      c("upper", "lower") %>%
      c(., "SOM_removed") %>%
      paste0(collapse = " + ")
  } else {
    if (i == 5) {
      cov_c_i <- cov_selected %>%
        c("upper", "lower") %>%
        c(., "year") %>%
        paste0(collapse = " + ")
    }
  }

  formula_i <- paste0(frac, " ~ ", cov_c_i) %>%
    as.formula()

  # trdat <- obs_top %>%
  #   filter(is.finite(.data[[frac]]))
  trdat <- obs %>%
    filter(is.finite(.data[[frac]])) %>%
    filter(!is.na(UTMX) & !is.na(UTMX)) %>%
    filter(lower > 0, upper < 200)

  # Three folds (placeholder)
  trdat %<>% mutate(
    fold = ceiling(fold/3)
  )
  holdout_i <- trdat %>%
    filter(fold == 4)
  trdat %<>% filter(fold < 4)
  
  if (!use_all_points) {
    trdat %<>% sample_n(n)
  }

  # Calculate weights
  # Calculate densities for depth intervals
  # dens <- ppp(
  #   trdat$UTMX,
  #   trdat$UTMY,
  #   c(441000, 894000),
  #   c(6049000, 6403000)
  # ) %>%
  #   density(
  #     sigma = 250,
  #     at = 'points',
  #     leaveoneout = FALSE
  #   )
  # 
  # attributes(dens) <- NULL
  # 
  # min(dens)
  # max(dens)
  # 
  # trdat %<>%
  #   mutate(
  #     density = dens,
  #     w = min(dens) / dens
  #   )
  
  # Weighting by depth intervals
  
  w_interval <- 10
  w_increment <- 1
  w_startdepth <- 0
  w_maxdepth <- 200
  w_iterations <- round(w_maxdepth / w_increment, digits = 0)
  
  dens_mat <- matrix(numeric(), nrow = nrow(trdat), ncol = w_iterations)
  
  for(j in 1:w_iterations)
  {
    upper_j <- w_startdepth + w_increment*(j - 1)
    lower_j <- upper_j + w_interval
    
    trdat_ind <- trdat$lower > upper_j & trdat$upper < lower_j
    trdat_ind[is.na(trdat_ind)] <- FALSE
    
    trdat_j <- trdat[trdat_ind, ]
    
    #sigma equal to the radius of a cirkle with an equal area per sample
    sigma_j <- sqrt(42951/(nrow(trdat_j)*pi))*1000
    
    # separate densities for wetlands and uplands
    
    dens_j <- numeric(nrow(trdat_j))
    
    for(k in 0:1) {
      trdat_j_wl_ind <- trdat_j$wetlands_10m == k
      
      trdat_jk <- trdat_j[trdat_j_wl_ind, ]
      
      dens_jk <- ppp(
        trdat_jk$UTMX,
        trdat_jk$UTMY,
        c(441000, 894000),
        c(6049000, 6403000)
      ) %>%
        density(
          sigma = sigma_j,
          at = 'points',
          leaveoneout = FALSE
        )
      
      attributes(dens_jk) <- NULL
      
      dens_j[trdat_j_wl_ind] <- dens_jk
    }
    # dens_j <- ppp(
    #   trdat_j$UTMX,
    #   trdat_j$UTMY,
    #   c(441000, 894000),
    #   c(6049000, 6403000)
    # ) %>%
    #   density(
    #     sigma = sigma_j,
    #     at = 'points',
    #     leaveoneout = FALSE
    #   )
    # 
    # attributes(dens_j) <- NULL
    
    dens_mat[trdat_ind, j] <- dens_j
  }
  
  dens_depth <- apply(
    dens_mat,
    1,
    function(x) {
      out <- mean(x, na.rm = TRUE)
      return(out)
    }
  )
  
  w_depth <- min(dens_depth, na.rm = TRUE) / dens_depth
  
  w_depth[!is.finite(w_depth)] <- 1
  
  trdat$w <- w_depth
  
  # List of folds
  
  folds_i <- lapply(
    unique(trdat$fold),
    function(x) {
      out <- trdat %>%
        mutate(
          is_j = fold != x,
          rnum = row_number(),
          ind_j = is_j*rnum
        ) %>%
        filter(ind_j != 0) %>%
        dplyr::select(., ind_j) %>%
        unlist() %>%
        unname()
    }
  )

  showConnections()

  # cl <- makePSOCKcluster(10)
  # registerDoParallel(cl)
  
  # xgboost optimization
  # 1: Fit learning rate (eta) and nrounds
  print("Step 1")

  set.seed(1)

  models[[i]] <- caret::train(
    form = formula_i,
    data = trdat,
    method = "xgbTree",
    na.action = na.pass,
    tuneGrid = tgrid,
    trControl = trainControl(
      index = folds_i,
      savePredictions = "final",
      predictionBounds = c(bounds_lower[i], bounds_upper[i]),
      summaryFunction = WeightedSummary,
      allowParallel = FALSE
    ),
    metric = 'RMSEw',
    maximize = FALSE,
    weights = trdat$w,
    num_parallel_tree = trees_per_round,
    objective = objectives[i]
  )

  # registerDoSEQ()
  # rm(cl)
  
  if (extra_tuning_xgb) {
    # 2: Fit max_depth and min_child_weight
    print("Step 2")
    
    set.seed(1)
    
    model2 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = max_depth_test,  # NB
        min_child_weight = min_child_weight_test,  # NB
        gamma = models[[i]]$bestTune$gamma,        
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = WeightedSummary,
        allowParallel = FALSE
      ),
      metric = 'RMSEw',
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )
    
    # 3: Tune gamma
    print("Step 3")
    
    set.seed(1)
    
    model3 <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = model2$bestTune$max_depth,
        min_child_weight = model2$bestTune$min_child_weight,
        gamma = gamma_test,  # NB
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = WeightedSummary,
        allowParallel = FALSE
      ),
      metric = 'RMSEw',
      maximize = FALSE,
      weights = trdat$w,
      num_parallel_tree = trees_per_round,
      objective = objectives[i]
    )
    
    models[[i]] <- model3
  }
  print(models[[i]])
  
  saveRDS(
    models[[i]],
    paste0(dir_results, "/model_", frac, ".rds")
  )
}

# 4: Reduce learning rate, use cv to optimize nrounds, early stopping

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
  out <- x$results %>%
    filter(RMSEw == min(RMSEw))
  return(out)
}
) %>%
  bind_rows() %>%
  mutate(
    Fraction = fractions,
    .before = 1
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

getpred <- function(x2, i2) {
  df <- x2$pred %>%
    arrange(rowIndex) %>%
    distinct(rowIndex, .keep_all = TRUE) %>%
    dplyr::select(., c(pred, obs)) %>%
    mutate(
      pred = ifelse(pred < 0, 0, pred)
    )
  # if (i2 > 4) df %<>% exp
  df %<>% mutate(
    fraction = fractions[i2],
    upper = quantile(obs, 0.99)
    ) %>%
    filter(obs < upper) %>%
    filter(pred < upper) %>%
    filter(obs >= 0)
  return(df)
}

allpred <- foreach(i = 1:6, .combine=rbind) %do%
  getpred(models[[i]], i)

allpred$fraction %<>% factor(levels = fractions)

levels(allpred$fraction) <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
  )

tiff(
  paste0(dir_results, "/accuracy_test", testn, ".tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

allpred %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(alpha = .01, shape = 16) +
  facet_wrap(~ fraction, nrow = 2, scales = "free") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(col = "red") +
  geom_blank(aes(y = upper)) +
  geom_blank(aes(x = upper)) +
  geom_blank(aes(y = 0)) +
  geom_blank(aes(x = 0)) +
  xlab("Observation (%)") +
  ylab("Prediction (%)")

dev.off()

# Accuracy vs depth

models[[5]]$trainingData
m5_w <- models[[5]]$pred %>% arrange(rowIndex) %>% select(weights)

plot(m5_w$weights, -models[[5]]$trainingData$upper)

i <- 5


d_acc <- sapply(
  0:199,
  function(x) {
    mdata <- models[[i]]$pred %>%
      bind_cols(models[[i]]$trainingData)
    
    ddat <- mdata %>% filter(
      upper < x + 1 & lower > x
    )
    
    out <- get_R2w(
      select(ddat, pred, obs),
      ddat$weights
    )
    
    return(out)
  }
)

d_wsum <- sapply(
  0:199,
  function(x) {
    mdata <- models[[i]]$pred %>%
      bind_cols(models[[i]]$trainingData)
    
    ddat <- mdata %>% filter(
      upper < x + 1 & lower > x
    )
    
    out <- sum(ddat$weights)
    
    return(out)
  }
)

# Covariate importance

library(colorRamps)
library(rcartocolor) # for colorblind palette

mycolors <- carto_pal(12, "Safe") %>% sort()

library(TSP)

l <- list()

ntop <- 20

for(i in 1:length(models))
{
  l[[i]] <- varImp(models[[i]])$importance %>%
    as_tibble(rownames = "covariate") %>%
    drop_na %>%
    arrange(- Overall) %>%
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
    category = ifelse(
      category == "basic",
      scorpan,
      category
    )
  )

l %<>%
  left_join(l_cat)

l %<>%
  ungroup() %>%
  arrange(target, Overall) %>%
  mutate(order = row_number())

l$category %<>% as.factor()

# levels(l$category) <- c(
#   "Bare soil",
#   "Spatial position",
#   "Parent materials",
#   "Topography",
#   "S2 time series",
#   "Soil"
# )

catcolors <- l$category %>%
  levels() %>%
  length() %>%
  carto_pal(., "Safe")
names(catcolors) <- levels(l$category)
colScale <- scale_fill_manual(name = "category", values = catcolors)


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
    ~ target,
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

dev.off()

# 10: Make maps for the test area

dir_cov_10km <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

predfolder <- dir_results %>%
  paste0(., "/predictions_testarea/") %T>%
  dir.create()

source("f_predict_passna.R")

# Make the maps

uppers <- c(0, 25, 50, 100)
lowers <- c(25, 50, 100, 200)

map_spec <- expand_grid(
  fraction_i = 1:6,
  interval = 1:4
  )

showConnections()

numCores <- 20

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
  c("uppers",
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
        year = 2010,
        upper = uppers[map_spec$interval[x]],
        lower = lowers[map_spec$interval[x]]
      ),
      n_const = 4,
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

for(i in 1:length(fractions)) {
  maps_10_km[[i]] <- c(1:4) %>%
    paste0(
      predfolder, "/", fractions[i],
      "_depth", .,
      ".tif"
    ) %>% rast()
  names(maps_10_km[[i]]) <- paste0(
    fraction_names[i], " ", uppers, " - ", lowers, " cm"
    )
}

# SOC depth distribution is very obviously wrong. I will need to fix it.
# I should probably apply a wider spatial filter to reduce the influence of
# the sinks data points.
# Maybe I should use a model to estimate density using covariates, as there is
# a higher concentration in some areas.


# Looking at 10 km maps

library(viridisLite)
library(tidyterra)

try(dev.off())

autoplot(maps_10_km[[1]]) +
  scale_fill_gradientn(colours = viridis(100), na.value = NA)

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
}
)
dev.off()

# maps_10km_s2 <- c(maps_10km[[1]], maps_10km[[2]], maps_10km[[3]], exp(maps_10km[[5]])/0.568, exp(maps_10km[[6]]))

maps_10km_s2 <- c(maps_10km[[1]], maps_10km[[2]], maps_10km[[3]], maps_10km[[5]]/0.568, maps_10km[[6]])

maps_10km_jb <- lapp(maps_10km_s2, JB) %>% as.factor()

myrgb <- col2rgb(mycolors)
tsp <- as.TSP(dist(t(myrgb)))
set.seed(1)
sol <- solve_TSP(tsp, control = list(repetitions = 1e3))
ordered_cols <- mycolors[sol]

# ggplot2::qplot(x = 1:12, y = 1, fill = I(ordered_cols), geom = "col", width = 1) + ggplot2::theme_void()

tiff(
  paste0(dir_results, "/JB_test", testn, ".tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

plot(
  maps_10km_jb,
  col = ordered_cols[levels(maps_10km_jb)[[1]]$ID],
  main = "JB-nummer"
)

dev.off()

# To do:
# Further refine xgboost use
# Use all observations
# Implement forward feature selection?
# Include forest samples
# Include profile database
# Include depth

# END