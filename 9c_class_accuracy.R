# 9c: Analyze class accuracies

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)
library(xgboost)
library(stringr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")

dir_boot_models <- dir_boot %>%
  paste0(., "/models/")

# Load observations used in models

obs_texture <- paste0(dir_results, "/observations_texture.rds") %>%
  readRDS()

# Load mean bootstrap predictions

boot_predictions_mean <- dir_boot %>%
  paste0(., "/models_boot_predictions_mean.rds") %>%
  readRDS()

# Load model weights

fractions_alt  <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")
fractions      <- fractions_alt

models_weights <- list()
models_indices <- list()

for (i in 1:length(fractions)) {
  frac <- fractions[i]
  
  models_weights[[i]] <- dir_results %>%
    paste0(., "/models_weights_", frac, ".rds") %>%
    readRDS()
  
  models_indices[[i]] <- dir_results %>%
    paste0(., "/models_indices_", frac, ".csv") %>%
    readRDS()
}

# 1: Accuracy for SOC classes

# Load combination raster showing parts covered by Tørv2022 and my own map.
SOC_combination_map <- root %>%
  paste0(., "/Tekstur2024_kort/Kulstof_arealkombination.tif") %>%
  rast()

# Extract combination

SOC_combination <- obs_texture %>%
  vect(geom = c("UTMX", "UTMY"), crs = mycrs) %>%
  terra::extract(x = SOC_combination_map, y = .)

# Combine predictions, observations, weights, and combination

SOC_df <- obs_texture %>%
  select(c(ID_new, db, ID_old, upper, lower, SOC, imputed, fold)) %>%
  mutate(
    observed = SOC,
    predicted = boot_predictions_mean[, 5],
    w = models_weights[[5]],
    combination = SOC_combination$Band_1,
    mean_depth = (upper + lower) / 2,
    class_pred = cut(predicted, breaks = c(0, 6, 12, 100)),
    class_obs = cut(observed, breaks = c(0, 6, 12, 100)),
  ) %>%
  filter(
    is.finite(SOC),
    is.finite(predicted),
    is.finite(w),
    imputed == FALSE,
    upper < 31,
    combination == 2,
    fold == 10,
    is.finite(mean_depth),
    mean_depth <= 30
  )

plot(SOC_df$observed, SOC_df$predicted)


library(MachineShop)  # weighted confusion matrix

SOC_acc_w <- SOC_df %>%
  select(c(class_pred, class_obs, w)) %>%
  na.omit() %>%
  confusion(
    x = .[, 1],
    y = .[, 2],
    weights = .[, 3]
  )

# 2: Compare accuracy with 2014 maps

texture_validation_obs <- obs_texture %>%
  mutate(mean_depth = (upper + lower) / 2) %>%
  filter(
    db != "Danish Soil Classification",
    db != "SINKS",
    db != "Profile database",
    fold == 10,
    upper < 31,
    is.finite(mean_depth),
    mean_depth <= 30
  ) %>%
  select(
    c(ID_new, db, ID_old, UTMX, UTMY, upper, lower, clay, silt, fine_sand,
      coarse_sand, SOC, mask_LU
    )
  )
  
texture_2014_maps <- c(
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/aclaynor.tif",
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/asiltnor.tif",
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/afsandno.tif",
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/agsandno.tif",
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/akulstof.tif"
) %>%
  rast()

names(texture_2014_maps) <- fractions_alt[1:5]

extract_2014 <- texture_validation_obs %>%
  vect(geom = c("UTMX", "UTMY"), crs = mycrs) %>%
  terra::extract(x = texture_2014_maps, y = .)

predictions_2024_val <- boot_predictions_mean[obs_texture$ID_new %in% texture_validation_obs$ID_new, ]

weights_val_2024 <- models_weights %>%
  bind_cols()

colnames(weights_val_2024) <- fractions

weights_val_2024 <- weights_val_2024[obs_texture$ID_new %in% texture_validation_obs$ID_new, ]

df_tex_obs_val <- texture_validation_obs %>%
  select(c(ID_new, db, clay, silt, fine_sand, coarse_sand, SOC)) %>%
  mutate(type = "observed")
df_tex_val_w <- weights_val_2024[1:5] %>%
  mutate(
    type = "w",
    db = df_tex_obs_val$db,
    ID_new = df_tex_obs_val$ID_new
  )
df_tex_pred_2014 <- extract_2014 %>%
  select(c(clay, silt, fine_sand, coarse_sand, SOC)) %>%
  mutate(
    type = "predicted",
    db = df_tex_obs_val$db,
    ID_new = df_tex_obs_val$ID_new
    )
df_tex_pred_2024 <- predictions_2024_val[, 1:5] %>%
  as.data.frame() %>%
  mutate(
    type = "predicted",
    db = df_tex_obs_val$db,
    ID_new = df_tex_obs_val$ID_new
  )

df_val_2014 <- bind_rows(
    df_tex_obs_val,
    df_tex_val_w,
    df_tex_pred_2014
  ) %>%
  mutate(year = 2014)
df_val_2024 <- bind_rows(
    df_tex_obs_val,
    df_tex_val_w,
    df_tex_pred_2024
  ) %>%
  mutate(year = 2024)

df_val_all <- bind_rows(df_val_2014, df_val_2024) %>%
  pivot_longer(
    c(clay, silt, fine_sand, coarse_sand, SOC),
    values_to = "value",
    names_to = "fraction"
    ) %>%
  pivot_wider(
    names_from = "type",
    values_from = "value",
    id_cols = c(db, year, fraction, ID_new)
  ) %>%
  mutate(
    db = case_when(
      db != "SEGES" ~ "Forests",
      .default = db
      )
  )

source("f_weighted_summaries.R")

df_val_all %>%
  drop_na() %>%
  group_by(db, year, fraction) %>%
  summarise(
    n = n(),
    r2w = round(get_R2w(cbind(predicted, observed), w), digits = 3),
    rmsew = round(get_RMSEw(cbind(predicted, observed), w), digits = 3)
  ) %>%
  mutate(
    fraction = factor(
      fraction,
      levels = fractions[1:5]
    )
  ) %>%
  arrange(db, fraction)


# 3: Uncertainty vs accuracy for texture classes

source("f_classify_soil_JB.R")

JB_observed <- classify_soil_JB(
  clay = obs_texture$clay,
  silt = obs_texture$silt,
  sand_f = obs_texture$fine_sand,
  SOM = obs_texture$SOC / 0.6,
  CaCO3 = obs_texture$CaCO3
)

JB_predicted <- classify_soil_JB(
  clay = boot_predictions_mean[, 1],
  silt = boot_predictions_mean[, 2],
  sand_f = boot_predictions_mean[, 3],
  SOM = boot_predictions_mean[, 5] / 0.6,
  CaCO3 = boot_predictions_mean[, 6]
)

mean_weights <- models_weights %>% 
  bind_cols() %>%
  apply(., 1, function(x) { mean(x, na.rm = TRUE) })

# Load bootstrap predictions

models_boot_predictions <- list()

for (i in 1:length(fractions)) {
  frac <- fractions[i]
  models_boot_predictions[[i]] <- dir_boot %>%
    paste0(., "/models_boot_predictions_", frac, ".rds") %>%
    readRDS()
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

JB_boot_unc <- lapply(
  1:100,
  function(x) {
    JB <- classify_soil_JB(
      clay = models_boot_predictions[[1]][, x],
      silt = models_boot_predictions[[2]][, x],
      sand_f = models_boot_predictions[[3]][, x],
      SOM = models_boot_predictions[[5]][, x] / 0.6,
      CaCO3 = models_boot_predictions[[6]][, x]
    )
    return(JB)
  }
) %>%
  bind_cols() %>%
  as.matrix() %>%
  apply(
    .,
    1,
    function(x2) {
      out <- sum(x2 != Mode(x2))
      return(out)
    }
    )

# Load table with numeric texture class differences

Texture_class_num_diff <- dir_dat %>%
  paste0(., "/Texture_class_num_diff.csv") %>%
  read.table(
    header = TRUE,
    sep = ";"
  ) %>%
  as.matrix()

Texture_class_num_diff <- Texture_class_num_diff[, -1]

JB_correct <- JB_predicted == JB_observed

JB_dev <- sapply(
  1:length(JB_predicted),
  function(x) {
    out <- Texture_class_num_diff[JB_predicted[x], JB_observed[x]]
    return(out)
  }
  ) %>%
  unname()
  

JB_unc_df <- data.frame(
  predicted = JB_predicted,
  observed = JB_observed,
  w = mean_weights,
  correct = JB_correct,
  dev = JB_dev,
  unc = JB_boot_unc,
  dataset = obs_texture$fold == 10
) %>%
  mutate(
    unc_int = cut(unc, breaks = c(0:9)*10 - 0.5) %>% as.numeric()
  )

# Percent correct by uncertainty interval

JB_unc_df %>%
  group_by(dataset, unc_int) %>%
  summarise(
    acc_w = weighted.mean(correct, w, na.rm = TRUE),
    mean_dev_w = weighted.mean(dev, w, na.rm = TRUE)
  )

JB_unc_df %>%
  group_by(dataset, dev, unc_int) %>%
  summarise(
    sum_w = sum(w, na.rm = TRUE)
  ) %>%
  pivot_wider(
    id_cols = c(dataset, unc_int),
    names_from = dev,
    values_from = sum_w
  )

# END