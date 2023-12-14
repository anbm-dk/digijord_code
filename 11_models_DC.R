# 09: Models for soil drainage classes

# TO DO:
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

# 1: Observations
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
  # values() %>%
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

model_DCfac

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

pred_DC_sd <- app(
  pred_DC_prob,
  function(x) {
    out <- spatstat.geom::weighted.var(1:5, x, na.rm = TRUE)^0.5
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


library(tidyterra)

names(pred_DC_prob) <- c(
  "Very well-drained",
  "Well-drained",
  "Moderately well-drained",
  "Poorly drained",
  "Very poorly drained"
)

tiff(
  paste0(dir_results, "/pred_DC_prob_10km_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

autoplot(pred_DC_prob) +
  facet_wrap(~lyr, ncol = 3) + 
  scale_fill_gradientn(colours = cividis(100), na.value = NA) +
  ggtitle("Soil drainage class probability")

try(dev.off())

tiff(
  paste0(dir_results, "/pred_DC_10km_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

autoplot(pred_DC) +
  scale_fill_gradientn(colours = rev(cividis(5)), na.value = NA) +
  ggtitle("Soil drainage class")

try(dev.off())

tiff(
  paste0(dir_results, "/pred_DC_mean_10km_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

autoplot(pred_DC_mean) +
  scale_fill_gradientn(colours = rev(cividis(100)), na.value = NA) +
  ggtitle("Mean soil drainage class")

try(dev.off())

tiff(
  paste0(dir_results, "/pred_DC_sd_10km_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

autoplot(pred_DC_sd) +
  scale_fill_gradientn(colours = inferno(100), na.value = NA) +
  ggtitle("Soil drainage class uncertainty (SD)")

try(dev.off())

tiff(
  paste0(dir_results, "/pred_DC_median_10km_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

autoplot(pred_DC_median) +
  scale_fill_gradientn(colours = rev(cividis(100)), na.value = NA) +
  ggtitle("Median soil drainage class")

try(dev.off())

tiff(
  paste0(dir_results, "/varImp_modelDC_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)

plot(varImp(model_DCfac$model), top = 15)

try(dev.off())

# Plot for covariates

cov_plot_names <- c(
  "chelsa_bio12_1981_2010_10m",
  "cwl_10m_fuzzy",
  "dhm2015_terraen_10m",
  "filled_s1_baresoil_composite_vh_8_days",
  "filled_s2_geomedian_b12",
  "fuzzy_geology_7",
  "slope_height",
  "vdtochn"
)

cov_mapplot_10km <- cov_10km %>%
  subset(cov_plot_names)

tiff(
  paste0(dir_results, "/cov_10km_test_", testn, ".tiff"),
  width = 23,
  height = 10,
  units = "cm",
  res = 300
)

plot(
  cov_mapplot_10km,
  col = viridis(100),
  nr = 2,
  axes = FALSE,
  box = TRUE,
  mar = c(0.1, 0.1, 1, 4),
  cex.main = 0.8
  )

try(dev.off())

# END