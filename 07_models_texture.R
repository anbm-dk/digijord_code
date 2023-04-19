# 07: Train texture models

# Cubist
# Use all observations, including NA (OK)
# Fractions:
# - Clay
# - Silt
# - Fine sand
# - Coarse sand
# - SOC (log)
# - CaCO3 (log)
# Topsoil
# Start with a small map for the test area
# Expand from there

# 1: Start up

library(Cubist)
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

# Test 1 - 8: Cubist
# Test 9: xgboost
testn <- 9
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


# 4: Load covariate data

dir_cov <- dir_dat %>% paste0(., "/covariates")

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230323.csv") %>%
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
    paste0(., "/dsc_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )
  
  SEGES_extr <- dir_extr %>%
    paste0(., "/SEGES_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )
}

SINKS_extr <- dir_extr %>%
  paste0(., "/SINKS_extr.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )


# 6: Merge data and transform the target variables

obs <- list(dsc, SEGES, SINKS) %>%
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

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)

bounds_lower <- c(0, 0, 0, 0, NA, NA)
bounds_upper <- c(100, 100, 100, 100, log(100), log(100))


# 7: Make training data

folds <- bind_rows(
  dsc_folds,
  SEGES_folds,
  SINKS_folds
)

names(folds) <- "fold"

extr <- bind_rows(
  dsc_extr,
  SEGES_extr,
  SINKS_extr
)

obs <- cbind(obs, extr, folds)

obs_top <- obs %>%
  filter(
    upper == 0,
    is.finite(fold)
    )

obs_top_v <- obs_top %>% vect(geom = c("UTMX", "UTMY"))

library(viridisLite)

tiff(
  paste0(dir_results, "/obs_map_test", testn, ".tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

plot(obs_top_v, "clay", breaks = 5, breakby = "cases", col = cividis(5))

dev.off()

plot(obs_top_v, "clay", breaks = 5, breakby = "cases", col = cividis(5))

# 8: Set up models

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  select(name) %>%
  unlist() %>%
  unname()

# Tuning grid

# For cubist:
# tgrid <- data.frame(
#   committees = 20,
#   neighbors = 0
# )

# For xgboost

tgrid <- expand.grid(
  nrounds = seq(100, 200, 50),
  eta = seq(0.1, 0.3, 0.1),
  max_depth = 6,
  min_child_weight = 1,
  gamma = 0,
  colsample_bytree = 0.8,
  subsample = 0.8
)

# Template for custom eval
# evalerror <- function(preds, dtrain) {
#   labels <- getinfo(dtrain, "label")
#   err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
#   return(list(metric = "error", value = err))
# }

# Weighted RMSE
RMSEw <- function(d, w)
{
  sqe <- w*(d[, 1] - d[, 2])^2
  msqe <- sum(sqe)/sum(w)
  out <- sqrt(msqe)
  return(out)
}

# Weighted R^2
R2w <- function(d, w)
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
  RMSEw <- function(d, w) {
    sqe <- w*(d[, 1] - d[, 2])^2
    msqe <- sum(sqe)/sum(w)
    out <- sqrt(msqe)
    return(out)
  }
  out[1] <- RMSEw(data[, 1:2], data$weights)
  
  # Weighted R^2
  require(boot)
  out[2] <- boot::corr(
    data[, 1:2],
    data$weights
  )^2
  names(out) <- c('RMSEw', 'R2w')
  return(out)
}

# Small random sample for testing
# Remember to include full dataset in the final model
n <- 1000

use_all_points <- TRUE
# use_all_points <- FALSE

# NB: Weighted cubist

l <- getModelInfo("cubist")

cubist_weighted <- l$cubist

cubist_weighted$label <- "cubist_weighted"

cubist_weighted$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  if(!is.null(wts)) {
    out <- Cubist::cubist(x,
                          y,
                          committees = param$committees,
                          weights = wts,
                          ...)
  } else {
    out <- Cubist::cubist(x, y, committees =  param$committees,  ...)
  }
  if(last) out$tuneValue$neighbors <- param$neighbors
  out
}

cubist_weighted$predict <- function(modelFit, newdata, submodels = NULL) {
  out <- predict(modelFit,
                 as.data.frame(newdata),
                 neighbors = modelFit$tuneValue$neighbors)
  if(!is.null(submodels)) {
    tmp <- vector(mode = "list", length = nrow(submodels) + 1)
    tmp[[1]] <- out
    
    for(j in seq(along = submodels$neighbors))
      tmp[[j+1]] <- predict(modelFit,
                            as.data.frame(newdata),
                            neighbors = submodels$neighbors[j])
    
    out <- tmp
  }
  out
}

cubist_weighted$tags <- c(l$cubist$tags, "Accepts Case Weights")


# 9: Train models

extra_tuning_xgb <- TRUE

models <- list()
models2 <- list()
models3 <- list()

for (i in 1:length(fractions))
{
  frac <- fractions[i]

  print(frac)

  cov_c_i <- cov_selected %>% paste0(collapse = " + ")

  if (i %in% 1:4) {
    cov_c_i <- cov_selected %>%
      c(., "SOM_removed") %>%
      paste0(collapse = " + ")
  } else {
    if (frac == "logSOC") {
      cov_c_i <- cov_selected %>%
        c(., "year") %>%
        paste0(collapse = " + ")
    }
  }

  formula_i <- paste0(frac, " ~ ", cov_c_i) %>%
    as.formula()

  trdat <- obs_top %>%
    filter(is.finite(.data[[frac]]))

  if (!use_all_points) {
    trdat %<>% sample_n(n)
  }

  # Calculate weights
  dens <- ppp(
    trdat$UTMX,
    trdat$UTMY,
    c(441000, 894000),
    c(6049000, 6403000)
  ) %>%
    density(
      sigma = 250,
      at = 'points',
      leaveoneout = FALSE
    )
  
  attributes(dens) <- NULL
  
  min(dens)
  max(dens)

  trdat %<>%
    mutate(
      density = dens,
      w = min(dens) / dens
    )
  
  # Three folds (placeholder)
  trdat %<>% mutate(
    fold = ceiling(fold/3)
  )
  holdout_i <- trdat %>%
    filter(fold == 4)
  trdat %<>% filter(fold < 4)
  
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
        select(ind_j) %>%
        unlist() %>%
        unname()
    }
  )

  showConnections()

  # cl <- makePSOCKcluster(10)
  # registerDoParallel(cl)
  
  # xgboost optimization
  # 1: Fit learning rate (eta) and nrounds

  set.seed(1)

  models[[i]] <- caret::train(
    form = formula_i,
    data = trdat,
    method = "xgbTree",
    # method = cubist_weighted,
    # method = "cubist",
    na.action = na.pass,
    tuneGrid = tgrid,
    trControl = trainControl(
      index = folds_i,
      savePredictions = "final",
      predictionBounds = c(bounds_lower[i], bounds_upper[i]),
      summaryFunction = WeightedSummary
    ),
    metric = 'RMSEw',
    maximize = FALSE,
    weights = trdat$w,
    allowParallel = FALSE  # for xgboost
  )

  # registerDoSEQ()
  # rm(cl)
  
  saveRDS(
    models[[i]],
    paste0(dir_results, "/model_", frac, ".rds")
  )
  
  if (extra_tuning_xgb) {
    
    # 2: Fit max_depth and min_child_weight
    set.seed(1)
    
    models2[[i]] <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = seq(1, 11, 2),
        min_child_weight = c(1, 2, 4, 8, 16, 32),
        gamma = 0,
        colsample_bytree = 0.8,
        subsample = 0.8
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = WeightedSummary
      ),
      metric = 'RMSEw',
      maximize = FALSE,
      weights = trdat$w,
      allowParallel = FALSE  # for xgboost
    )
    
    # 3: Tune gamma, 0 to 0.4
    set.seed(1)
    
    models3[[i]] <- caret::train(
      form = formula_i,
      data = trdat,
      method = "xgbTree",
      na.action = na.pass,
      tuneGrid = expand.grid(
        nrounds = models[[i]]$bestTune$nrounds,
        eta = models[[i]]$bestTune$eta,
        max_depth = models2[[i]]$bestTune$max_depth,
        min_child_weight = models2[[i]]$bestTune$min_child_weight,
        gamma = seq(0, 0.4, 0.1),
        colsample_bytree = 0.8,
        subsample = 0.8
      ),
      trControl = trainControl(
        index = folds_i,
        savePredictions = "final",
        predictionBounds = c(bounds_lower[i], bounds_upper[i]),
        summaryFunction = WeightedSummary
      ),
      metric = 'RMSEw',
      maximize = FALSE,
      weights = trdat$w,
      allowParallel = FALSE  # for xgboost
    )
    
    saveRDS(
      models3[[i]],
      paste0(dir_results, "/model_opt_", frac, ".rds")
    )
    
    models3[[i]]
  }
}

models_loaded <- lapply(
  1:6,
  function(x) {
    out <- fractions[x] %>%
      paste0(dir_results, "/model_", ., ".rds") %>%
      readRDS()
    return(out)
  }
)

# models
models3
models <- models3

# 4: Reduce learning rate, use cv to optimize nrounds, early stopping

# models <- models_loaded

names(models) <- fractions

models %>%
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

# Inspect models

get_acc <- function(x2, i2) {
  df <- x2$pred %>%
    arrange(rowIndex) %>%
    distinct(rowIndex, .keep_all = TRUE) %>%
    select(c(pred, obs, weights))
  
  if (i2 > 4) df %<>% exp
  
  df %<>% bind_cols(x2$trainingData)
  
  r2_all <- df %$% R2w(cbind(pred, obs), weights)
  
  r2_bare <- df %>%
    filter(!is.na(s2_geomedian_b2)) %$%
    R2w(cbind(pred, obs), weights)
  
  r2_covered <- df %>%
    filter(is.na(s2_geomedian_b2)) %$%
    R2w(cbind(pred, obs), weights)
  
  rmse_all <- df %$% RMSEw(cbind(pred, obs), weights)
  
  rmse_bare <- df %>%
    filter(!is.na(s2_geomedian_b2)) %$%
    RMSEw(cbind(pred, obs), weights)
  
  rmse_covered <- df %>%
    filter(is.na(s2_geomedian_b2)) %$%
    RMSEw(cbind(pred, obs), weights)
  
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
    select(c(pred, obs))
  if (i2 > 4) df %<>% exp
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

outfolder <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

cov_10km <- outfolder %>%
  list.files(full.names = TRUE) %>%
  rast()

predfolder <- dir_dat %>%
  paste0(., "/testarea_10km/predictions_", testn, "/") %T>%
  dir.create()

source("f_predict_passna.R")

# Make the maps

maps_10km <- list()

showConnections()

for(i in 1:length(fractions))
{
  frac <- fractions[i]
  
  maps_10km[[i]] <- predict(
    cov_10km,
    models[[i]],
    fun = predict_passna,
    na.rm = FALSE,
    filename = paste0(predfolder, frac,  "_10km.tif"),
    overwrite = TRUE,
    const = data.frame(
      SOM_removed = TRUE,
      year = 2010
    )
  )
}

maps_10km <- predfolder %>%
  paste0(., fractions,  "_10km.tif") %>%
  rast()

names(maps_10km) <- fractions

# Looking at 10 km maps

library(viridisLite)

plot(maps_10km, col = cividis(100))

maps_10km_stack2 <- c(
  maps_10km[[1:4]],
  exp(maps_10km[[5]]),
  exp(maps_10km[[6]])
)

names(maps_10km_stack2) <- fraction_names

tiff(
  paste0(dir_results, "/maps_test", testn, ".tiff"),
  width = 24,
  height = 16,
  units = "cm",
  res = 300
)

plot(maps_10km_stack2, col = cividis(100))

dev.off()

JB <- function(clay, silt, sand_f, SOM, CaCO3)
{
  out <- rep(0, length(clay))
  out[CaCO3 > 10] <- 12
  out[out == 0 & SOM > 10] <- 11
  out[out == 0 & clay < 5 & silt < 20 & sand_f < 50] <- 1
  out[out == 0 & clay < 5 & silt < 20] <- 2
  out[out == 0 & clay < 10 & silt < 25 & sand_f < 40] <- 3
  out[out == 0 & clay < 10 & silt < 25]<-4
  out[out == 0 & clay < 15 & silt < 30 & sand_f < 40] <- 5
  out[out == 0 & clay < 15 & silt < 30] <- 6
  out[out == 0 & clay < 25 & silt < 35] <- 7
  out[out == 0 & clay < 45 & silt < 45] <- 8
  out[out == 0 & silt < 50] <- 9
  out[out == 0] <- 10
  return(out)
}

maps_10km_s2 <- c(maps_10km[[1]], maps_10km[[2]], maps_10km[[3]], exp(maps_10km[[5]])/0.568, exp(maps_10km[[6]]))

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
# Implement forward feature selection
# Include forest samples
# Include profile database
# Include depth

# END