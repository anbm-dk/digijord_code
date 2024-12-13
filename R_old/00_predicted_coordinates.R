# Predicted coordinates and residuals

library(terra)
library(caret)
library(xgboost)
library(dplyr)
library(magrittr)
library(foreach)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

mycrs <- "EPSG:25832"

tmpfolder <- paste0(dir_dat, "/Temp/")

terraOptions(tempdir = tmpfolder)

# Load covariates

dir_cov <- dir_dat %>% paste0(., "/covariates")
cov_files <- dir_cov %>% list.files(full.names = TRUE)
cov_names <- cov_files %>%
  basename() %>%
  tools::file_path_sans_ext()

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230323.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  select(name) %>%
  unlist() %>%
  unname()

source("f_predict_passna.R")

# Make layers with and col

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- rast(cov_files[dem_ind])

# row <- init(
#   dem,
#   fun = "row",
#   filename = paste0(tmpfolder, "row_raw.tif"),
#   names = "row",
#   overwrite = TRUE,
#   datatype = "INT2U"
# )
# col <- init(
#   dem,
#   fun = "col",
#   filename = paste0(tmpfolder, "col_raw.tif"),
#   names = "col",
#   overwrite = TRUE,
#   datatype = "INT2U",
#   gdal = c("TILED=YES", "BLOCKYSIZE=8000", "BLOCKXSIZE=16")
# )

row_raw <- rast(paste0(tmpfolder, "row_raw.tif"))
col_raw <- rast(paste0(tmpfolder, "col_raw.tif"))

# Deselect oblique geographic coordinates

cov_selected <- grep(
  paste0("ogc", collapse = "|"),
  cov_selected,
  invert = TRUE,
  value = TRUE
)

cov <- rast(cov_files)
cov <- subset(cov, cov_selected)

# Sample dem for row and col

row_col_dem <- c(row, col, dem)

set.seed(1234)

# pts <- spatSample(row_col_dem, 10^5, na.rm = TRUE, cells = TRUE, xy = TRUE)

# pts <- select(pts, !dhm2015_terraen_10m)
#
# extr <- extract(
#   cov,
#   pts$cell
# )
#
# pts_2 <- bind_cols(pts, extr)

dir_extr <- dir_dat %>%
  paste0(., "/extracts/")

# write.table(
#   pts_2,
#   file = paste0(dir_extr, "/random_100k.csv"),
#   sep = ";",
#   row.names = FALSE
# )

pts_2 <- read.table(
  paste0(dir_extr, "/random_100k.csv"),
  header = TRUE,
  sep = ";"
)

# Models for predicting row and col

dir_xy_models <- dir_dat %>%
  paste0(., "/models_xy/") %T>%
  dir.create()

# Small random sample for testing
# Remember to include full dataset in the final models
n <- 1000

use_all_points <- TRUE
# use_all_points <- FALSE

# Tuning grid
# For xgboost

tgrid <- expand.grid(
  nrounds = 20,
  eta = seq(0.1, 1, 0.1),
  max_depth = 6,
  min_child_weight = 1,
  gamma = 0,
  colsample_bytree = 0.5,
  subsample = 0.5
)

gamma_test <- seq(0, 0.5, 0.1)
max_depth_test <- seq(1, 20, 3)
min_child_weight_test <- c(1, 2, 4, 8, 16, 32)

trees_per_round <- 10

# For xgblinear

# 9: Train models

extra_tuning_xgb <- FALSE

models <- list()

targets <- c("row", "col")

for (i in 1:length(targets))
{
  targ <- targets[i]

  print(targ)

  cov_c_i <- cov_selected %>% paste0(collapse = " + ")

  formula_i <- paste0(targ, " ~ ", cov_c_i) %>%
    as.formula()

  trdat <- pts_2

  if (!use_all_points) {
    set.seed(11111)
    trdat %<>% sample_n(n)
  }

  set.seed(789)

  folds <- select(trdat, any_of(targ)) %>%
    unlist() %>%
    createFolds(., 3)

  showConnections()

  # xgboost optimization
  # 1: Fit learning rate (eta) and nrounds

  print("Step 1")
  set.seed(123)

  models[[i]] <- caret::train(
    form = formula_i,
    data = trdat,
    method = "xgbLinear",
    na.action = na.pass,
    # tuneGrid = tgrid,
    trControl = trainControl(
      index = folds,
      savePredictions = "final",
      allowParallel = FALSE
    )
    # ,
    # num_parallel_tree = trees_per_round
  )

  if (extra_tuning_xgb) {
    # 2: Fit max_depth and min_child_weight
    print("Step 2")
    set.seed(123)

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
        index = folds,
        savePredictions = "final"
      )
      # ,
      # num_parallel_tree = trees_per_round
    )

    # 3: Tune gamma
    print("Step 3")
    set.seed(123)

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
        gamma = gamma_test, # NB
        colsample_bytree = models[[i]]$bestTune$colsample_bytree,
        subsample = models[[i]]$bestTune$subsample
      ),
      trControl = trainControl(
        index = folds,
        savePredictions = "final"
      ),
      num_parallel_tree = trees_per_round
    )

    models[[i]] <- model3
  }
  print(models[[i]])

  saveRDS(
    models[[i]],
    paste0(dir_xy_models, "/model_", targ, ".rds")
  )
}

models_loaded <- lapply(
  1:2,
  function(x) {
    out <- targets[x] %>%
      paste0(dir_xy_models, "/model_", ., ".rds") %>%
      readRDS()
    return(out)
  }
)

# models <- models_loaded

names(models) <- targets

# Model summary

models_sum <- lapply(models, function(x) {
  out <- x$results %>%
    filter(RMSE == min(RMSE))
  return(out)
}) %>%
  bind_rows() %>%
  mutate(
    target = targets,
    .before = 1
  ) %T>%
  write.table(
    file = paste0(dir_xy_models, "/models_sum.csv"),
    sep = ";",
    row.names = FALSE
  )

models_sum

# Covariate importance

library(tidyr)
library(tibble)

imp_all <- models %>%
  seq_along() %>%
  lapply(
    function(x) {
      out <- varImp(models[[x]])$importance %>%
        as.data.frame() %>%
        rownames_to_column(var = "covariate") %>%
        mutate(target = names(models)[x])
      return(out)
    }
  ) %>%
  bind_rows() %>%
  pivot_wider(
    id_cols = covariate,
    names_from = target,
    values_from = Overall
  ) %>%
  rowwise() %>%
  mutate(mean_imp = mean(c_across(-covariate))) %>%
  arrange(-mean_imp) %T>%
  write.table(
    file = paste0(dir_xy_models, "/var_imp.csv"),
    sep = ";",
    row.names = FALSE
  )

imp_all

# Inspect models

# get_acc <- function(x2, i2) {
#   df <- x2$pred %>%
#     arrange(rowIndex) %>%
#     distinct(rowIndex, .keep_all = TRUE) %>%
#     select(c(pred, obs, weights)) %>%
#     mutate(
#       pred = ifelse(pred < 0, 0, pred)
#     )
#
#   # if (i2 > 4) df %<>% exp
#
#   df %<>% bind_cols(x2$trainingData)
#
#   r2_all <- df %$% get_R2w(cbind(pred, obs), weights)
#
#   r2_bare <- df %>%
#     filter(!is.na(s2_geomedian_b2)) %$%
#     get_R2w(cbind(pred, obs), weights)
#
#   r2_covered <- df %>%
#     filter(is.na(s2_geomedian_b2)) %$%
#     get_R2w(cbind(pred, obs), weights)
#
#   rmse_all <- df %$% get_RMSEw(cbind(pred, obs), weights)
#
#   rmse_bare <- df %>%
#     filter(!is.na(s2_geomedian_b2)) %$%
#     get_RMSEw(cbind(pred, obs), weights)
#
#   rmse_covered <- df %>%
#     filter(is.na(s2_geomedian_b2)) %$%
#     get_RMSEw(cbind(pred, obs), weights)
#
#   out <- data.frame(
#     r2_all,
#     r2_bare,
#     r2_covered,
#     rmse_all,
#     rmse_bare,
#     rmse_covered
#   )
#
#   return(out)
# }
#
# library(foreach)
#
# acc_all <- foreach(i = 1:2, .combine = rbind) %do%
#   get_acc(models[[i]], i)
#
# acc_all %<>% mutate(target = targets, .before = 1)
#
# write.table(
#   acc_all,
#   paste0(dir_xy_models, "/acc_xy.csv"),
#   sep = ";",
#   row.names = FALSE
# )

getpred <- function(x2, i2) {
  df <- x2$pred %>%
    arrange(rowIndex) %>%
    distinct(rowIndex, .keep_all = TRUE) %>%
    select(c(pred, obs)) %>%
    mutate(
      pred = ifelse(pred < 0, 0, pred)
    )
  df %<>% mutate(
    target = targets[i2],
    upper = quantile(obs, 0.99)
  ) %>%
    filter(obs < upper) %>%
    filter(pred < upper) %>%
    filter(obs >= 0)
  return(df)
}

library(foreach)

allpred <- foreach(i = 1:2, .combine = rbind) %do%
  getpred(models[[i]], i)

allpred$target %<>% factor(levels = targets)

tiff(
  paste0(dir_xy_models, "/accuracy_xy.tiff"),
  width = 15,
  height = 10,
  units = "cm",
  res = 300
)

allpred %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(alpha = .01, shape = 16) +
  facet_wrap(~target, nrow = 1, scales = "free") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(col = "red") +
  geom_blank(aes(y = upper)) +
  geom_blank(aes(x = upper)) +
  geom_blank(aes(y = 0)) +
  geom_blank(aes(x = 0)) +
  xlab("Observation (m)") +
  ylab("Prediction (m)")

dev.off()

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
    mutate(target = targets[i]) %>%
    mutate(rank = 1:ntop)
}

l %<>% bind_rows() %>%
  mutate(
    target = factor(
      target,
      levels = targets
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
  paste0(dir_xy_models, "/importance_xy.tiff"),
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
    ncol = 2,
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
  paste0(., "/testarea_10km/predictions_xy/") %T>%
  dir.create()

source("f_predict_passna.R")

# Make the maps

maps_10km <- list()

showConnections()

for (i in 1:length(targets))
{
  targ <- targets[i]

  maps_10km[[i]] <- predict(
    cov_10km,
    models[[i]],
    fun = predict_passna,
    na.rm = FALSE,
    filename = paste0(predfolder, targ, "_10km.tif"),
    overwrite = TRUE,
    n_digits = 0
  )
}

maps_10km <- predfolder %>%
  paste0(., targets, "_10km.tif") %>%
  rast()

names(maps_10km) <- targets

# Looking at 10 km maps

library(viridisLite)

plot(maps_10km, col = cividis(100), box = FALSE)

maps_10km_stack2 <- maps_10km

names(maps_10km_stack2) <- targets

tiff(
  paste0(dir_xy_models, "/maps_xy.tiff"),
  width = 24,
  height = 16,
  units = "cm",
  res = 300
)

plot(maps_10km_stack2, col = cividis(100))

dev.off()

# Make maps

library(foreach)
library(parallel)

numCores <- detectCores()
numCores

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

subdir_tiles <- dir_tiles %>%
  list.dirs() %>%
  .[-1]

dir_pred_all <- dir_xy_models %>%
  paste0(., "/predictions/") %T>%
  dir.create()

dir_pred_tiles <- dir_pred_all %>%
  paste0(., "/tiles/") %T>%
  dir.create()

n_digits <- 0

for (i in 1:length(targets)) {
  targ <- targets[i]

  dir_pred_tiles_targ <- dir_pred_tiles %>%
    paste0(., "/", targ, "/") %T>%
    dir.create()

  model_i <- models[[i]]

  showConnections()

  cl <- makeCluster(numCores)

  clusterEvalQ(
    cl,
    {
      library(terra)
      library(caret)
      library(xgboost)
      library(magrittr)
      library(dplyr)
    }
  )

  clusterExport(
    cl,
    c(
      "i",
      "model_i",
      "subdir_tiles",
      "dir_pred_tiles_targ",
      "targ",
      "cov_names",
      "cov_selected",
      "predict_passna",
      "dir_dat",
      "n_digits"
    )
  )

  parSapplyLB(
    cl,
    1:length(subdir_tiles),
    function(x) {
      tmpfolder <- paste0(dir_dat, "/Temp/")

      terraOptions(memfrac = 0.02, tempdir = tmpfolder)

      cov_x <- subdir_tiles[x] %>%
        list.files(full.names = TRUE) %>%
        rast()

      names(cov_x) <- cov_names

      cov_x2 <- subset(cov_x, cov_selected)

      tilename_x <- basename(subdir_tiles[x])

      outname_x <- dir_pred_tiles_targ %>%
        paste0(., "/", targ, "_", tilename_x, ".tif")

      outmap <- predict(
        cov_x2,
        model_i,
        fun = predict_passna,
        na.rm = FALSE,
        filename = outname_x,
        overwrite = TRUE,
        n_digits = 0
      )

      # if (i > 4) {
      #   outmap2 <- terra::exp(outmap)
      #   outmap <- outmap2
      # }

      # terra::math(
      #   outmap,
      #   "round",
      #   digits = n_digits,
      #   filename = outname_x,
      #   overwrite = TRUE
      # )
    }
  )

  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)

  outtiles_targ <- dir_pred_tiles_targ %>%
    list.files(full.names = TRUE) %>%
    sprc()

  if (i == 2) {
    merge(
      outtiles_targ,
      filename = paste0(dir_pred_all, targ, "_merged.tif"),
      overwrite = TRUE,
      datatype = "INT2U",
      gdal = c("TILED=YES", "BLOCKYSIZE=8000", "BLOCKXSIZE=16")
    )
  } else {
    merge(
      outtiles_targ,
      filename = paste0(dir_pred_all, targ, "_merged.tif"),
      overwrite = TRUE,
      datatype = "INT2U"
    )
  }
}

# Subtract predicted row and col from real coordinates for residuals

row_predicted <- rast(paste0(dir_pred_all, "/row_merged.tif"))
col_predicted <- rast(paste0(dir_pred_all, "/col_merged.tif"))

s_row <- c(row_raw, row_predicted)
s_col <- c(col_raw, col_predicted)

f <- function(i) {
  out <- i[1] - i[2]
  return(out)
}

# terraOptions(memfrac = 0.4, tempdir = tmpfolder)
#
# app(s_row, f,  filename = paste0(dir_pred_all, "/row_residual.tif"),
#     overwrite = TRUE)
# app(s_col, f,  filename = paste0(dir_pred_all, "/col_residual.tif"),
#     overwrite = TRUE)

tile_shapes <- dir_tiles %>%
  paste0(., "/tiles.shp") %>%
  vect()

max_char <- length(tile_shapes) %>%
  1:. %>%
  as.character() %>%
  nchar() %>%
  max()

library(stringr)

tile_numbers <- length(tile_shapes) %>%
  1:. %>%
  str_pad(
    .,
    max_char,
    pad = "0"
  )

dir_row_tiles <- dir_pred_tiles %>% paste0(., "/row_residual/") %T>% dir.create()
dir_col_tiles <- dir_pred_tiles %>% paste0(., "/col_residual/") %T>% dir.create()

dir_residual_tiles <- c(dir_row_tiles, dir_col_tiles)

# Before: Pass file names for raw and predicted row and col to clusters.

files_raw <- paste0(tmpfolder, targets, "_raw.tif")

dir_pred_tiles_targets <- dir_pred_tiles %>%
  paste0(., "/", targets, "/")

# In parallel:
# 1: crop raw col and row rasters
# 2: Calculate residuals of predictions, based on same tile

# After: Merge residual rasters

library(parallel)

numCores <- detectCores()
numCores

for (i in 1:length(targets)) {
  pred_files_i <- dir_pred_tiles_targets[[i]] %>%
    list.files(full.names = TRUE)

  showConnections()

  cl <- makeCluster(numCores)

  clusterEvalQ(
    cl,
    {
      library(raster)
      library(terra)
      library(magrittr)
      library(dplyr)
      library(tools)
    }
  )

  clusterExport(
    cl,
    c(
      "dir_tiles",
      "tile_numbers",
      "files_raw",
      "pred_files_i",
      "tmpfolder",
      "f",
      "i",
      "dir_residual_tiles",
      "targets"
    )
  )

  parSapplyLB(
    cl,
    1:length(tile_shapes),
    function(j) {
      terraOptions(memfrac = 0.02, tempdir = tmpfolder)

      pred_tile <- rast(pred_files_i[j])

      tile_shapes <- dir_tiles %>%
        base::paste0(., "/tiles.shp") %>%
        terra::vect()

      r_raw <- rast(files_raw[i])

      raw_tile <- crop(
        r_raw,
        tile_shapes[j]
      )

      s <- c(raw_tile, pred_tile)

      outname <- dir_residual_tiles[i] %>%
        paste0(., "/", targets[i], "_residual_tile_", tile_numbers[j], ".tif")

      app(
        s,
        f,
        filename = outname,
        overwrite = TRUE,
        wopt = list(datatype = "INT2S")
      )
    }
  )
  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)

  outtiles_targ_res <- dir_residual_tiles[i] %>%
    list.files(full.names = TRUE) %>%
    sprc()

  merge(
    outtiles_targ_res,
    filename = paste0(dir_pred_all, targets[i], "_residual.tif"),
    overwrite = TRUE,
    datatype = "INT2S"
  )
}

# Write rasters

rs_residuals <- paste0(dir_pred_all, targets, "_residual.tif") %>% rast()

outrasters <- c(row_predicted, col_predicted, rs_residuals)

outnames <- c("row_predicted", "col_predicted", "row_residual", "col_residual")

dtypes <- c("INT2U", "INT2U", "INT2S", "INT2S")

names(outrasters) <- outnames

for (i in 1:length(outnames)) {
  writeRaster(
    outrasters[[i]],
    filename = paste0(dir_cov, "/", outnames[i], ".tif"),
    datatype = dtypes[i],
    names = outnames[i],
    gdal = c("TILED=YES", "BLOCKYSIZE=992", "BLOCKXSIZE=992"),
    overwrite = TRUE
  )
}

# End
