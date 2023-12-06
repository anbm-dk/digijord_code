# 09: Models for soil drainage

# TO DO:

# First: Predict soil drainage classes:
# Import soil drainage classes (get correct coordinates first) [ok]
# Extract covariates [ok]
# Include texture predictions as covariates [ok]
# Rearrange tiles for texture predictions to fit covariate structure
# Make summary function with weighted MAE for accuracy [ok]
# Calculate weights [ok]
# Train xgboost regression model [ok, poor accuracy]
# Train xgboost classification model [ok]
# Analyse results
# Make map for test area
# Make national map

# Secondly: Predict artificially drained areas
# Import and merge data:
# - 7 km grid (original coordinates) (or new LR pts?)
# - data from orbicon
# - data from Landskontoret for Planteavl
# Aggregate relevant covariates at field scale (including new texture maps and
# soil drainage classes)
# Split into tiles
# Extract covariates
# Make summary function for weighted accuracy (or weighted AUC?)
# Calculate weights (sum of weights should be equal for drained and undrained
# points)
# Train xgboost classification model
# Analyse results
# Make map for test area
# Make national map

# 1: Start up

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)
library(sf)
library(exactextractr)
library(party)
library(rpart)
library(doParallel)
library(spatstat) # weights
library(RODBC)
library(stringr) 

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

fraction_names_underscore <- c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand", "SOC", "CaCO3"
)

# Results folder

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") %T>%
  dir.create()

# Part 1: Soil drainage classes

profiles_shp <- dir_dat %>%
  paste0(
    .,
    "/observations/profiles/Profiles_coordinates_new/",
    "Profiles_coordinates_new.shp"
  ) %>%
  vect()

# 1.1: Correct coordinates for profiles

profiles_db <- dir_dat %>%
  paste0(., "/observations/profiles/DDJD2023.accdb")

con3 <- odbcConnectAccess2007(profiles_db)

profiles_DC <- sqlFetch(con3, "PROFIL") %>%
  select(c(PROFILNR, DRAENKL)) %>%
  drop_na() %>%
  filter(DRAENKL %in% 1:5)

profiles_DC %<>% inner_join(
  x = values(profiles_shp),
  y = .,
  by = "PROFILNR"
) %>% arrange(
  PROFILNR
) %>% vect(
  geom = c("x", "y"),
  crs = mycrs,
  keepgeom = TRUE
)

library(viridisLite)

plot(profiles_DC, "DRAENKL", col = rev(cividis(5)))


# 1.2:  Load covariates

cov_dir <- dir_dat %>% paste0(., "/covariates")
cov_files <- cov_dir %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()
cov <- paste0(cov_dir, "/", cov_files) %>%
  rast()
names(cov) <- cov_names
crs(cov) <- mycrs

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  dplyr::select(., name) %>%
  unlist() %>%
  unname()

cov %<>% subset(cov_selected)

# 1.3 Load texture predictions 

dir_pred_all <- dir_results %>%
  paste0(., "/predictions/")

tex_pred <- dir_pred_all %>% list.files(
  pattern = ".tif",
  full.names = TRUE
) %>%
  grep(".ovr", ., names(cov), value = TRUE, invert = TRUE) %>%
  grep(".aux.xml", ., names(cov), value = TRUE, invert = TRUE) %>%
  grep(".vat.cpg", ., names(cov), value = TRUE, invert = TRUE) %>%
  grep(".vat.dbf", ., names(cov), value = TRUE, invert = TRUE) %>%
  grep(
    pattern = paste(fraction_names_underscore, collapse = "|"),
    .,
    value = TRUE
  ) %>%
  rast()

cov_DC <- c(cov, tex_pred)

# obs_DC <- terra::extract(
#   cov_DC,
#   profiles_DC,
#   bind = TRUE,
#   ID = FALSE
# )

dir_extr <- dir_dat %>%
  paste0(., "/extracts/")

# obs_DC %>%
#   values() %>%
#   write.table(
#     file = paste0(dir_extr, "obs_DC_extr.csv"),
#     row.names = FALSE,
#     col.names = TRUE,
#     sep = ";"
#   )

obs_DC <- read.table(
  paste0(dir_extr, "obs_DC_extr.csv"),
  header = TRUE,
  sep = ";"
)

cov_DC_names <- names(cov_DC) %>%
  match(., names(obs_DC)) %>%
  names(obs_DC)[.]

# 3: Extract folds and mask

dir_folds <- dir_dat %>%
  paste0(., "/folds/")

file_folds_10_100m <- paste0(dir_folds, "/folds_10_100m.tif")

folds_10_100m <- rast(file_folds_10_100m)

folds_DC <- terra::extract(
  x = folds_10_100m,
  y = profiles_DC,
  ID = FALSE,
) %>% unlist() %>% unname()

mask_LU <- paste0(dir_dat, "/layers/Mask_LU.tif") %>% rast()

mask_LU_DC <- terra::extract(
  mask_LU,
  profiles_DC,
  ID = FALSE
) %>% unlist() %>% unname()

obs_DC %<>%
  values() %>%
  mutate(
    fold = folds_DC,
    mask_LU = mask_LU_DC,
    UTMX = x,
    UTMY = y
  ) %>%
  filter(
    !is.na(mask_LU)
  )

# Calculate weights
source("f_get_dens.R")
mean_dens_DC <- nrow(obs_DC) / (43107 * 10^6)
sigma_DC <- sqrt(43107 / (nrow(obs_DC) * pi)) * 1000
dens_DC <- get_dens(obs_DC, sigma_DC)
w_DC <- mean_dens_DC / dens_DC
w_DC[w_DC > 1] <- 1

obs_DC$w <- w_DC

# Three folds (placeholder)
trdat_DC <- obs_DC %>% 
  mutate(
    fold = ceiling(fold / 3)
  ) %>%
  filter(fold < 4)

trdat_indices_DC <- which(obs_DC$PROFILNR %in% trdat_DC$PROFILNR)

holdout_DC <- obs_DC %>% filter(fold == 10)
holdout_indices_DC <- which(obs_DC$PROFILNR %in% holdout_DC$PROFILNR)

# List of folds

folds_DC <- lapply(
  unique(trdat_DC$fold),
  function(fold_i) {
    out <- trdat_DC %>%
      mutate(
        rnum = row_number(),
      ) %>%
      filter(!(fold %in% fold_i)) %>%
      dplyr::select(., rnum) %>%
      unlist() %>%
      unname()
  }
)

# Set up model
source("f_optimize_xgboost.R")
source("f_weighted_summaries.R")

bounds_DC <- list(
  eta = c(0.1, 1),
  min_child_weight_sqrt = c(1, sqrt(100)),
  gamma_sqrt = c(0, sqrt(100)),
  colsample_bytree = c(0.1, 1),
  subsample = c(0.1, 1),
  colsample_bynode = c(0.1, 1),
  ogcs_index = c(1L, 7L),
  total_imp = c(0.5, 1)
)

tgrid <- expand.grid(
  nrounds = 10,
  eta = 0.3,
  max_depth = 6,
  min_child_weight = 1,
  gamma = 0,
  colsample_bytree = 0.75,
  subsample = 0.75
)

# Drainage classes as factor

trdat_DC %<>% mutate(
  DCfac = factor(DRAENKL, labels = paste0("DC", 1:5))
)

foreach::registerDoSEQ()
showConnections()

model_DCfac <- optimize_xgboost(
  target = "DCfac",  # character vector (length 1), target variable.
  cov_names = cov_DC_names,  # Character vector, covariate names,
  data = trdat_DC, # data frame, input data
  bounds_bayes = bounds_DC, # named list with bounds for bayesian opt.
  bounds_pred = c(FALSE, FALSE), # numeric, length 2, bounds for predicted values
  cores = 19, # number cores for parallelization
  trgrid = tgrid, # data frame with tuning parameters to be tested in basic model
  folds = folds_DC, # list with indices, folds for cross validation
  sumfun = WeightedSummary_DCfac, # summary function for accuracy assessment
  metric = "OAw", # character, length 1, name of evaluation metric
  max_metric = TRUE, # logical, should the evaluation metric be maximized
  weights = trdat_DC$w, # numeric, weights for model training and evaluation
  trees_per_round = 10, # numeric, length 1, number of trees that xgboost should train in each round
  obj_xgb = "multi:softprob", # character, length 1, objective function for xgboost
  colsample_bynode_basic = 0.75, # numeric, colsample_bylevel for basic model
  cov_keep = NULL, # Character vector, covariates that should always be present
  final_round_mult = 10,  # Multiplier for the number of rounds in the final model
  seed = 321,  # Random seed for model training,
  classprob = TRUE
)

model_DCfac$model %>% varImp()

# Add more iterations for bayes optimization?

# Load covariates for test area

dir_cov_10km <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

predfolder <- dir_results %>%
  paste0(., "/predictions_testarea/")

breaks <- c(0, 30, 60, 100, 200)

uppers <- breaks %>% rev() %>% .[-1] %>% rev()
lowers <- breaks %>% .[-1]

max_char <- breaks %>%
  as.character() %>%
  nchar() %>%
  max()

breaks_chr <- breaks %>%
  str_pad(
    .,
    max_char,
    pad = "0"
  )

depth_int_chr <- paste0(breaks_chr[1:4], "_", breaks_chr[2:5], "_cm")

cov_10km <- dir_cov_10km %>%
  list.files(full.names = TRUE) %>%
  rast() %>%
  subset(cov_selected)

# Load texture maps for test area

maps_10_km <- list()

for (i in 1:length(uppers)) {
  maps_10_km[[i]] <- paste0(
      predfolder, "/SOM_remov/depth_", i, "/sum/") %>%
    list.files(full.names = TRUE) %>%
    rast()
  names(maps_10_km[[i]]) <- paste0(
    fraction_names_underscore, "_", breaks_chr[i], "_", breaks_chr[i + 1], "_cm"
  )
}

texture_10km <- rast(maps_10_km)

cov_DCmodel_names <- rownames(varImp(model_DCfac$model)$importance)

cov_10km_DC <- c(cov_10km, texture_10km)

cov_10km_DC <- cov_10km_DC %>%
  terra::subset(., cov_DCmodel_names)

# Make predictions for the test area

source("f_predict_passna.R")

pred_DC <- predict(
  cov_10km_DC,
  model_DCfac$model,
  fun = predict_passna,
  na.rm = FALSE,
  const = data.frame(
    dummy = NA
  ),
  n_const = 0,
  n_digits = 0,
  # ,
  # filename = outname,
  # overwrite = TRUE
)

pred_DC_prob <- predict(
  cov_10km_DC,
  model_DCfac$model,
  fun = predict_passna_prob,
  na.rm = FALSE,
  const = data.frame(
    dummy = NA
  ),
  n_const = 0,
  n_digits = 2
  # ,
  # filename = outname,
  # overwrite = TRUE
)

pred_DC_mean <- app(
  pred_DC_prob,
  function(x) {
    out <- stats::weighted.mean(1:5, x, na.rm = TRUE)
    return(out)
  }
)

pred_DC_median <- app(
  pred_DC_prob,
  function(x) {
    out <- weighted.median(1:5, x, na.rm = TRUE, type = 4)
    return(out)
  }
)

plot(pred_DC, col = rev(cividis(100)))
plot(pred_DC_mean, col = rev(cividis(100)))
plot(pred_DC_median, col = rev(cividis(100)))
plot(pred_DC_prob, col = viridis(100))

# Part 2: Artificially drained areas



# 2 Aggregate covariates to field level
# Load field data (first used repair geometry in ArcGIS)

DK_fields <- dir_dat %>%
  paste0(., "fields_2022/Marker_2022_slut.shp") %>%
  vect()

# Deselect covariates that should not be aggregated at field level:
# - field data (imk)
# - climatic data (chelsa)
# - land use (lu)
# - oblique geographic coordinates (ogc_pi)
# - wetlands (wetlands_10m)
# - landscape elements (landscape_)
# - georegions (georeg_)
# - river valley bottoms (rvb_bios)
# - central wetlands (cwl_10m_fuzzy)

cov_exclude <- c(
  grep('imk_', names(cov), value = TRUE),
  grep('chelsa_', names(cov), value = TRUE),
  grep('lu_', names(cov), value = TRUE),
  grep('ogc_pi', names(cov), value = TRUE),
  grep('wetlands_10m', names(cov), value = TRUE),
  grep('landscape_', names(cov), value = TRUE),
  grep('georeg_', names(cov), value = TRUE),
  grep('rvb_bios', names(cov), value = TRUE),
  grep('cwl_10m_fuzzy', names(cov), value = TRUE)
)

cov_for_aggr <- setdiff(names(cov), cov_exclude)

# Extract to field level (mean)
# Count NA per field
# Rasterize
# Split into tiles
# Extract for both drainage datasets


# END of new script

# Old script

# library(devtools)
# 
# devtools::install_github("topepo/C5.0")
# library(Cubist)
# library(C50)

# 2 Load observations

obs_drain <- dir_dat %>%
  paste0(., "/observations/drainage_data/drainpts_2022.shp") %>%
  vect()

plot(obs_drain, "DR_N")


# 3: Extract folds

dir_folds <- dir_dat %>%
  paste0(., "/folds/")

file_folds_10_100m <- paste0(dir_folds, "/folds_10_100m.tif")

folds_10_100m <- rast(file_folds_10_100m)

folds_drain <- terra::extract(
  x = folds_10_100m,
  y = obs_drain,
  ID = FALSE,
)

write.table(
  folds_drain,
  paste0(dir_folds, "/folds_drain.csv"),
  row.names = FALSE,
  sep = ";"
)

folds_drain <- dir_folds %>%
  paste0(., "folds_drain.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  )

names(folds_drain) <- "fold"

# 4: Load covariates

cov_dir <- dir_dat %>% paste0(., "/covariates")

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ";",
    header = TRUE,
    fileEncoding = "latin1"
  )

cov_use <- cov_cats %>%
  filter(anbm_use == 1) %>%
  select(name) %>%
  unlist()

cov_files <- cov_dir %>% list.files()

cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_names[!cov_names %in% cov_cats$name]

cov <- paste0(cov_dir, "/", cov_files) %>%
  rast()

names(cov) <- cov_names

dir_extr <- dir_dat %>%
  paste0(., "/extracts/")

# drain_extr <- terra::extract(
#   x = cov,
#   y = obs_drain,
#   ID = FALSE,
# )
#
# write.table(
#   drain_extr,
#   paste0(dir_extr, "drain_extr.csv"),
#   row.names = FALSE,
#   sep = ";"
# )
#
# buffer_drain <- terra::buffer(
#   obs_drain,
#   width = 40
# ) %>%
#   st_as_sf()
#
# buffer_drain_extr <- exact_extract(
#   x = cov,
#   y = buffer_drain,
#   fun = "mean",
#   progress = TRUE
# )
#
# names(buffer_drain_extr) <- names(cov)
#
# write.table(
#   buffer_drain_extr,
#   paste0(dir_extr, "/buffer_drain_extr.csv"),
#   row.names = FALSE,
#   sep = ";"
# )

usebuffer <- TRUE

if (usebuffer) {
  drain_extr <- dir_extr %>%
    paste0(., "/buffer_drain_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )
} else {
  drain_extr <- dir_extr %>%
    paste0(., "/drain_extr.csv") %>%
    read.table(
      header = TRUE,
      sep = ";",
    )
}

# drain_extr1 <- drain_extr

# 6: Merge data

coords_drain <- geom(obs_drain) %>%
  as.data.frame() %>%
  select(c(x, y))

names(coords_drain) <- c("UTMX", "UTMY")

trdat_drain <- values(obs_drain) %>%
  cbind(., coords_drain, folds_drain, drain_extr)

trdat_drain$drained <- as.factor(trdat_drain$DR_N)

levels(trdat_drain$drained) <- c("No", "Yes")

# 7: Train models

cov_c <- cov_use %>%
  paste0(collapse = " + ")

formula_drain <- paste0("drained ~ ", cov_c) %>%
  as.formula()

trdat_drain %<>% filter(is.finite(fold))

# Calculate weights

dens_drain <- ppp(
  trdat_drain$UTMX,
  trdat_drain$UTMY,
  c(441000, 894000),
  c(6049000, 6403000)
) %>%
  density(
    sigma = 1000,
    at = "points",
    leaveoneout = FALSE
  )

attributes(dens_drain) <- NULL

trdat_drain %<>%
  mutate(
    density = dens_drain,
    w = min(dens_drain) / dens_drain
  )

folds_drain_list <- lapply(
  1:10,
  function(x) {
    out <- trdat_drain %>%
      mutate(
        is_j = fold != x,
        rnum = row_number(),
        ind_j = is_j * rnum
      ) %>%
      filter(ind_j != 0) %>%
      select(ind_j) %>%
      unlist() %>%
      unname()
  }
)

tgrid <- expand.grid(
  list(
    model = c("rules", "tree"),
    winnow = FALSE,
    trials = c(1, seq(10, 100, 10))
  ),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

showConnections()

cl <- makePSOCKcluster(10)
registerDoParallel(cl)

set.seed(1)

model_drain <- caret::train(
  form = formula_drain,
  data = trdat_drain,
  method = "C5.0",
  na.action = na.pass,
  tuneGrid = tgrid,
  trControl = trainControl(
    index = folds_drain_list,
    savePredictions = "final"
  )
)

registerDoSEQ()
rm(cl)

model_drain

# Without buffers
# C5.0
# 745 samples
# 121 predictors
# Resampling: Bootstrapped (10 reps)
#   model  trials  Accuracy   Kappa
# rules  100     0.6987109  0.3931187

# With buffers
# C5.0
# 745 samples
# 121 predictors
# Resampling: Bootstrapped (10 reps)
# Summary of sample sizes: 667, 668, 672, 658, 680, 667, ...
# Resampling results across tuning parameters:
#   model  trials  Accuracy   Kappa
# tree    70     0.6968081  0.3854105

model_drain2 <- C5.0(
  form = formula_drain,
  data = trdat_drain,
  trials = 100,
  model = "tree"
)

varImp(model_drain)
C5imp(model_drain$finalModel, metric = "splits", pct = FALSE)

saveRDS(
  model_drain,
  paste0(dir_results, "/model_drain.rds")
)

# END
