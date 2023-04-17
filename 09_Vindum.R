# 07: Vindum comparison

# 1: Start up

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 8
mycrs <- "EPSG:25832"

# Results folder

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

outfolder <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

predfolder <- dir_dat %>%
  paste0(., "/testarea_10km/predictions_", testn, "/")

# Load Vindum data

vindum_obs <- dir_dat %>%
  paste0(., "/observations/Vindum/Vindum_everything.csv") %>%
  read.table(
    header = TRUE,
    sep = ";"
  ) %>%
  filter(DEPTH == 25) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
    )

# Load predictions

predictions <- predfolder %>%
  list.files(full.names = TRUE) %>%
  rast()

names(predictions) <- predfolder %>%
  list.files() %>%
  file_path_sans_ext()

# Extract predictions

vindum_extr <- terra::extract(predictions, vindum_obs)

# Plot

plot(predictions[[1]], ext = ext(vindum_obs))
plot(vindum_obs, "LER", add = TRUE)

plot(vindum_extr$clay_10km, vindum_obs$LER)
abline(1,1)

cor(vindum_extr$clay_10km, vindum_obs$LER, use =  "pairwise.complete.obs")^2

plot(vindum_extr$logSOC_10km, log(vindum_obs$SOC))
abline(1,1)

plot(exp(vindum_extr$logSOC_10km), vindum_obs$SOC)
abline(1,1)

cor(exp(vindum_extr$logSOC_10km), vindum_obs$SOC, use =  "pairwise.complete.obs")^2

plot(exp(predictions[[5]]), ext = ext(vindum_obs))
plot(vindum_obs, "SOC", add = TRUE)

plot(vindum_extr$logSOC_10km, vindum_obs$SOC)
abline(1,1)

# END