# 09: Models for soil drainage

# TO DO:

# First: Predict soil drainage classes:
# Import soil drainage classes (get correct coordinates first)
# Extract coordinates
# Include texture predictions as covariates
# Rearrange tiles for texture predictions to fit covariate structure
# Make summary function with weighted MAE for accuracy
# Calculate weights
# Train xgboost regression model
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

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 13
mycrs <- "EPSG:25832"

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

# 1.4.2: Profiles: Texture

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

plot(profiles_DC, "DRAENKL")

# 2 Aggregate covariates to field level
# Load field data (first used repair geometry in ArcGIS)

DK_fields <- dir_dat %>%
  paste0(., "fields_2022/Marker_2022_slut.shp") %>%
  vect()
  
# Load covariates

cov_dir <- dir_dat %>% paste0(., "/covariates")
cov_files <- cov_dir %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()
cov <- paste0(cov_dir, "/", cov_files) %>%
  rast()
names(cov) <- cov_names
crs(cov) <- mycrs

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230712.csv") %>%
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
  paste0(., "/cov_categories_20230227.csv") %>%
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
