# Soil type majority percentage

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(ggridges)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") 

# Load table with zonal statistics

zonal_df <- root %>%
  paste0("/JB_majority.txt") %>%
  read.table(header = TRUE, sep = "\t")

# Recalculate zonal majority percentage

breaks_ha <- seq(0, 30, 10) %>% c(., Inf)

zonal_df %<>%
  mutate(
    majority_perc = MAJORITY_COUNT*100 / COUNT,
    ha = AREA / 10000,
    ha_class = cut(ha, breaks_ha, right = FALSE)
  )

zonal_summary <- zonal_df %>%
  group_by(ha_class) %>%
  summarise(ha_sum = sum(ha))

zonal_df <- zonal_df %>%
  left_join(zonal_summary, by = "ha_class") %>%
  mutate(w = ha/ha_sum)

mean(zonal_df$majority_perc)
# [1] 89.54849

weighted.mean(zonal_df$majority_perc, zonal_df$COUNT)
# [1] 88.22608

# Size interval plot

zonal_df %>%
  ggplot(
    aes(
      y = ha_class,
      x = majority_perc,
      fill = ha_class
    )
  ) +
  xlim(c(NA, 100)) +
  geom_density_ridges(
    aes(height = after_stat(density),
        weight = w),
    stat = "density"
  )

tiff(
  paste0(dir_results, "/texture_class_freq_size_", testn, ".tiff"),
  width = 16,
  height = 13,
  units = "cm",
  res = 600
)

zonal_df %>%
  ggplot(
    aes(
      x = majority_perc,
      # y = after_stat(density),
      weight = w
    )
  ) +
  xlim(c(NA, 100)) +
  geom_histogram(breaks = seq(0, 100, 10), color = "black", fill = "#1B9E77") +
  facet_wrap(~ ha_class) +
  labs(title = "Marker inddelt efter størrelse (ha)",
       y = "Hyppighed",
       x = "Dominerende jordbundstype (Procent af marken)"
       )

try(dev.off())

tiff(
  paste0(dir_results, "/texture_class_freq_", testn, ".tiff"),
  width = 16,
  height = 13,
  units = "cm",
  res = 600
)

zonal_df %>%
  mutate(w2 = ha/sum(ha)) %>%
  ggplot(
    aes(
      x = majority_perc,
      # y = after_stat(density),
      weight = w2
    )
  ) +
  xlim(c(NA, 100)) +
  geom_histogram(breaks = seq(0, 100, 10), color = "black", fill = "#1B9E77") +
  labs(title = "Alle marker",
       y = "Del af samlet areal",
       x = "Dominerende jordbundstype (Procent af marken)"
  )

try(dev.off())

# END