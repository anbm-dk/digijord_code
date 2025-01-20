# Residuals fo SOC topsoil peat (>6%)

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)
library(xgboost)
library(stringr)

library(doParallel)
library(spatstat) # weights

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

source("f_predict_passna.R")

testn <- 14
mycrs <- "EPSG:25832"

# Results folder

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

obs_all <- readRDS(file = paste0(dir_results, "/observations_texture.rds"))

obs_top <- obs_all %>%
  filter(
    upper < 30,
    is.finite(SOC),
    db != "SINKS"
  ) %>%
  mutate(
    ispeat = SOC > 6
  ) %>%
  vect(geom=c("UTMX", "UTMY"), crs=mycrs, keepgeom=TRUE)

soc_map <- rast("C:/Users/au542768/digijord/Soil_maps_10m_new/Kulstof2022/SOC_000_030_cm/SOC_mean_000_030_cm_combined.tif")

landscapes <- vect("O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Landscapeelements/nyland.shp")
	
library(tidyterra)

obs_soc_extr <- terra::extract(x = soc_map, y = obs_top, ID = FALSE, bind = TRUE) %>%
  mutate(SOC_pred = Band_1) %>%
  filter(SOC_pred > 6)
  
plot(obs_soc_extr, "ispeat")

obs_soc_extr2 <- terra::extract(x = landscapes, y = obs_soc_extr)

obs_soc_extr3 <- bind_cols(values(obs_soc_extr), obs_soc_extr2)

values(obs_soc_extr) <- obs_soc_extr3

obs_soc_extr3 %>%
  group_by(M_FO_TEKST) %>%
  summarise(mean = mean(ispeat)) %>%
  ggplot(aes(y = mean, x = M_FO_TEKST)) +
  geom_col()

obs_soc_extr3 %>%
  group_by(db) %>%
  summarise(mean = mean(ispeat)) %>%
  ggplot(aes(y = mean, x = db)) +
  geom_col()

obs_soc_extr3 %>%
  mutate(centrallavbund = cwl_10m_fuzzy > 0.5) %>%
  group_by(centrallavbund, M_FO_TEKST) %>%
  summarise(mean = mean(ispeat)) %>%
  ggplot(aes(y = mean, x = M_FO_TEKST, fill = centrallavbund)) +
  geom_col(position = position_dodge())
  

# Log residual

obs_soc_extr %<>%
  mutate(logSOC_res = log(SOC) - log(SOC_pred)) %>%
  mutate(logSOC_res = case_when(
    is.finite(logSOC_res) ~ logSOC_res,
    .default = 0
  ))


obs_soc_extr %>%
  values() %>%
  group_by(M_FO_TEKST) %>%
  summarise(mean = mean(logSOC_res)) %>%
  ggplot(aes(y = mean, x = M_FO_TEKST)) +
  geom_col()

obs_soc_extr %>%
  values() %>%
  group_by(db) %>%
  summarise(mean = mean(logSOC_res)) %>%
  ggplot(aes(y = mean, x = db)) +
  geom_col()

obs_soc_extr %>%
  values() %>%
  mutate(centrallavbund = cwl_10m_fuzzy > 0.5) %>%
  group_by(centrallavbund, M_FO_TEKST) %>%
  summarise(mean = mean(logSOC_res)) %>%
  ggplot(aes(y = mean, x = M_FO_TEKST, fill = centrallavbund)) +
  geom_col(position = position_dodge())

writeVector(obs_soc_extr, paste0(dir_results, "/soc_residuals_peat.shp"))

obs_soc_extr %>%
  values() %>%
  mutate(centrallavbund = cwl_10m_fuzzy > 0.5) %>%
  group_by(M_FO_TEKST, centrallavbund) %>%
  summarise(n = n())

obs_soc_extr %>%
  values() %>%
  mutate(centrallavbund = cwl_10m_fuzzy > 0.5) %>%
  group_by(db) %>%
  summarise(n = n())

# END