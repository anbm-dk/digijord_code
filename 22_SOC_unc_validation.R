# Validate SOC uncertainties


library(parallel)
library(caret)
library(terra)
library(magrittr)
library(dplyr)
library(xgboost)
library(foreach)
library(stringr)
library(tools) # file_path_sans_ext

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

frac_ind_mineral <- c(1:4)
frac_ind_predict <- c(1:length(fractions))[-3]  # Exclude fine sand

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)

fraction_names_underscore <- c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand", "SOC", "CaCO3"
)

dir_cov <- dir_dat %>% paste0(., "/covariates")
cov_files <- dir_cov %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

# For peat probabilities and confidence intervals

soc6pct_log <- log(6)

prob_q_out <- c(0.5, 2.5, 5.0, 16.0, 84.0, 95.0, 97.5, 99.5)/100
prob_q_out8 <- prob_q_out
prob_q_out_chr <- prob_q_out %>%
  multiply_by(1000) %>%
  formatC(width = 4, flag = "0") %>%
  paste0("p", .)

# Calculate mean logSOC prediction and MSE for logsoc

obs_texture <- paste0(dir_results, "/observations_texture.rds") %>%
  readRDS()

dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")

models_boot_predictions_soc <- dir_boot %>%
  paste0(., "/models_boot_predictions_SOC.rds") %>%
  readRDS()

models_weights_soc <- dir_results %>%
  paste0(., "/models_weights_SOC.rds") %>%
  readRDS()

logsoc_mean_prediction <- models_boot_predictions_soc %>%
  log() %>%
  apply(., 1, mean)

# Load combination raster showing parts covered by Tørv2022 and my own map.
# Combination == 1 shows peat2022 extent
SOC_combination_map <- root %>%
  paste0(
    ., "/Soil_maps_10m_new/Kulstof2022/SOC_000_030_cm/",
    "SOC_combination_000_030_cm.tif") %>%
  rast()

# Extract combination

SOC_combination <- obs_texture %>%
  vect(geom = c("UTMX", "UTMY"), crs = mycrs) %>%
  terra::extract(x = SOC_combination_map, y = .)

breaks <- c(0, 30, 60, 100, 200)

logSOC_df <- obs_texture %>%
  select(c(ID_new, db, ID_old, upper, lower, SOC, imputed, fold)) %>%
  mutate(
    SOC = case_when(
      SOC == 0 ~ 0.01831564,
      .default = SOC
    ),
    observed = log(SOC),
    predicted = logsoc_mean_prediction,
    w = models_weights_soc,
    combination = SOC_combination$Band_1,
    indices = factor((!fold == 10) + 1, labels = c("Holdout", "CV")),
    mean_d = (upper + lower)/2,
    depth = cut(mean_d, breaks, include.lowest = TRUE)
  ) %>%
  filter(
    is.finite(SOC),
    is.finite(predicted),
    is.finite(w),
    imputed == FALSE,
    is.finite(mean_d),
    is.finite(depth)
  )

# MSE by depth - all samples

library(MetricsWeighted)

logsoc_mse_all <- logSOC_df %>%
  group_by(
    indices, depth
  ) %>%
  summarise(
    # r2w = round(get_R2w(cbind(predicted, observed), w), digits = 3),
    # rmsew = round(get_RMSEw(cbind(predicted, observed), w), digits = 3),
    msew = MetricsWeighted::mse(observed, predicted, w = w)
  ) %>%
  arrange(desc(indices))
logsoc_mse_all
#   indices depth      msew
#   <fct>   <fct>     <dbl>
#   1 CV    [0,30]    0.247
# 2 CV      (30,60]   0.834
# 3 CV      (60,100]  1.81 
# 4 CV      (100,200] 1.73 
# 5 Holdout [0,30]    0.283
# 6 Holdout (30,60]   0.898
# 7 Holdout (60,100]  1.57 
# 8 Holdout (100,200] 2.10

logsoc_mse_splitpeat2022 <- logSOC_df %>%
  group_by(
    indices, depth, combination
  ) %>%
  summarise(
    # r2w = round(get_R2w(cbind(predicted, observed), w), digits = 3),
    # rmsew = round(get_RMSEw(cbind(predicted, observed), w), digits = 3),
    msew = MetricsWeighted::mse(observed, predicted, w = w)
  ) %>%
  arrange(combination, desc(indices))
logsoc_mse_splitpeat2022
#    indices depth     combination  msew
#    <fct>   <fct>           <int> <dbl>
# 1  CV      [0,30]              1 0.656
# 2  CV      (30,60]             1 1.92 
# 3  CV      (60,100]            1 3.01 
# 4  CV      (100,200]           1 2.59 
# 5  Holdout [0,30]              1 0.779
# 6  Holdout (30,60]             1 2.44 
# 7  Holdout (60,100]            1 2.85 
# 8  Holdout (100,200]           1 4.35 
# 9  CV      [0,30]              2 0.150
# 10 CV      (30,60]             2 0.618
# 11 CV      (60,100]            2 1.57 
# 12 CV      (100,200]           2 1.61 
# 13 Holdout [0,30]              2 0.164
# 14 Holdout (30,60]             2 0.617
# 15 Holdout (60,100]            2 1.26 
# 16 Holdout (100,200]           2 1.69 
# 17 CV      [0,30]             NA 0.304
# 18 CV      (60,100]           NA 1.10 

logsoc_mse_depths <- logsoc_mse_all %>%
  ungroup() %>%
  filter(indices == "CV") %>%
  select(msew) %>%
  unlist() %>%
  unname()
logsoc_mse_depths

# Accidentally used these values, when creating the prob and q maps for the
# topsoil (discovered 2024-11-08):
# [1] 0.2830142 0.8975956 1.5693857 2.1013693

# Should have used these ones:
# [1] 0.2472204 0.8343481 1.8066069 1.7327846

# Follow-up (2024-11-11): I re-ran the analyses, and apparently this difference
# is so small that the maps do not change (at least within the specified 
# number of decimal places)

# Function for log mean and variance

calc_log_mean_pse <- function(x, mse_log) {
  x[x == 0] <- 0.01831564
  x %<>% log()
  out <- c(
    mean(x, na.rm = TRUE),
    sqrt(var(x, na.rm = TRUE) + mse_log)
  )
  return(out)
}

calc_prob_q <- function(x, q) {
  out <- pnorm(
    q = q,
    mean = x[1],
    sd = x[2],
    lower.tail = FALSE
  ) %>%
    round(., digits = 3) %>%
    multiply_by(100)
  return(out)
}

calc_q_log <- function (x) {
  out <- qnorm(
    p = prob_q_out,
    mean = x[1],
    sd = x[2]
  ) %>%
    exp() %>%
    round(., digits = 1) %>%
    set_names(prob_q_out_chr)
  out[out < 0] <- 0
  out[out > 60] <- 60
  return(out)
}

# Prediction variance for logsoc

logsoc_pvar <- models_boot_predictions_soc %>%
  log() %>%
  apply(., 1, var)

logSOC_df %<>% mutate(
  mse = logsoc_mse_depths[as.numeric(depth)],
  pvar = logsoc_pvar[obs_texture$ID_new %in% logSOC_df$ID_new],
  pse = sqrt(mse + pvar)
)

# Peat probability for validation sample

peat_probability_all <- logSOC_df %>%
  select(predicted, pse) %>%
  apply(., 1, function(x) {
    out <- calc_prob_q(x, q = log(6))
    return(out)
  }
  )

prob_breaks <- seq(0, 100, 10)

logSOC_df %<>% mutate(
  peat_probability = peat_probability_all,
  prob_interval = cut(
    peat_probability,
    prob_breaks, 
    include.lowest = TRUE, 
    right = FALSE
  ),
  ispeat = observed > log(6)
)

logSOC_df %>%
  group_by(depth, indices, prob_interval) %>%
  summarise(
    f_peat = weighted.mean(ispeat, w = w, na.rm = TRUE)*100,
    n = n()
    ) %>%
  arrange(indices) %>%
  print(n = 100)

#QQ for 100 quantiles

qs100 <- seq(0.005, 0.995, 0.01)

prob_q_out <- qs100
prob_q_out_chr <- prob_q_out %>%
  multiply_by(1000) %>%
  formatC(width = 4, flag = "0") %>%
  paste0("p", .)

soc_qs100 <- logSOC_df %>%
  select(predicted, pse) %>%
  apply(., 1, function(x) {
    out <- calc_q_log(x)
    return(out)
  }
  ) %>%
  t()

library(tidyr)

soc_p_under100 <- logSOC_df %>%
  bind_cols(soc_qs100) %>%
  pivot_longer(
    cols = all_of(prob_q_out_chr),
    names_to = "quantile"
    ) %>%
  group_by(indices, depth, quantile) %>%
  summarise(
    p_under = weighted.mean(exp(observed) < value, w = w)*100,
    n = n()
    ) %>%
  arrange(indices, depth)
soc_p_under100


tiff(
  paste0(dir_results, "/SOC_quantiles_test", testn, ".tiff"),
  width = 10,
  height = 10,
  units = "cm",
  res = 300
)

ggplot(
  soc_p_under100,
  aes(x = rep(qs100*100, times = 8), y = p_under, col = indices)
) + 
  geom_line() +
  facet_wrap(~depth) +
  coord_equal() +
  geom_abline(slope = 1, intercept = 1) +
  ylab("Percent of data") +
  xlab("Prediction quantile (%)")

try(dev.off())
try(dev.off())

# Only eight specific quantiles

prob_q_out <- prob_q_out8
prob_q_out_chr <- prob_q_out %>%
  multiply_by(1000) %>%
  formatC(width = 4, flag = "0") %>%
  paste0("p", .)

soc_qs8 <- logSOC_df %>%
  select(predicted, pse) %>%
  apply(., 1, function(x) {
    out <- calc_q_log(x)
    return(out)
  }
  ) %>%
  t()

soc_p_under8 <- logSOC_df %>%
  bind_cols(soc_qs8) %>%
  pivot_longer(
    cols = all_of(prob_q_out_chr),
    names_to = "quantile"
  ) %>%
  group_by(indices, depth, quantile) %>%
  summarise(
    p_under = weighted.mean(exp(observed) < value, w = w)*100,
    n = n()
  ) %>%
  pivot_wider(names_from = quantile, values_from = p_under) %>%
  arrange(indices, depth)
soc_p_under8




# END