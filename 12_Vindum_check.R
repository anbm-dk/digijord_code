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

testn <- 14
mycrs <- "EPSG:25832"

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

# Results folder

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

predfolder <- dir_results %>%
  paste0(., "/predictions_testarea/SOM_remov/depth_1/sum/")

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

names(predictions) <- fractions

# Extract predictions

vindum_extr <- terra::extract(predictions, vindum_obs, ID = FALSE)[, 1:5]

vindum_extr %<>% mutate(
  source = "prediction"
)

# Plot

vindum_fractions <- vindum_obs %>%
  values() %>%
  mutate(
    clay = LER,
    silt = SILT,
    fine_sand = C_SILT + C_FSAND + FSAND,
    coarse_sand = C_SAND
  ) %>%
  select(any_of(fractions)) %>%
  mutate(
    source = "observation"
  )

library(dplyr)

vindum_comparison <- bind_rows(
  vindum_fractions,
  vindum_extr
) %>%
  mutate(
    id = rep(1:nrow(vindum_fractions), times = 2)
  ) %>%
  pivot_longer(
    cols = -c(source, id),
    names_to = "fraction",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols = c(fraction, id),
    names_from = source
  ) %>%
  mutate(
    fraction = factor(fraction, levels = fractions[1:5])
    ) %>%
  as.data.frame()

tiff(
  paste0(dir_results, "/vindum_test_", testn, ".tiff"),
  width = 23,
  height = 10,
  units = "cm",
  res = 300
)

ggplot(vindum_comparison, aes(x = observation, y = prediction)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~ fraction) +
  coord_equal() +
  geom_abline()

try(dev.off())

vindum_comparison %>%
  drop_na() %>%
  group_by(fraction) %>%
  summarise(
    r2 = cor(observation, prediction, use = "pairwise.complete.obs")^2,
    rmse = RMSE(prediction, observation)
  ) %>%
  as.data.frame()

library(viridisLite)

tiff(
  paste0(dir_results, "/vindum_pred_test_", testn, ".tiff"),
  width = 23,
  height = 10,
  units = "cm",
  res = 300
)

plot(predictions, y = 1:5, ext = ext(vindum_obs), col = cividis(100))

try(dev.off())

plot(predictions[[1]], ext = ext(vindum_obs))
plot(vindum_obs, "LER", add = TRUE)

plot(vindum_extr$clay_10km, vindum_obs$LER)
abline(1, 1)

cor(vindum_extr$clay_10km, vindum_obs$LER, use = "pairwise.complete.obs")^2

plot(vindum_extr$logSOC_10km, log(vindum_obs$SOC))
abline(1, 1)

plot(exp(vindum_extr$logSOC_10km), vindum_obs$SOC)
abline(1, 1)

cor(exp(vindum_extr$logSOC_10km), vindum_obs$SOC, use = "pairwise.complete.obs")^2

plot(exp(predictions[[5]]), ext = ext(vindum_obs))
plot(vindum_obs, "SOC", add = TRUE)

plot(vindum_extr$logSOC_10km, log(vindum_obs$SOC))
abline(1, 1)

# END
