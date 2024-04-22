# 16: Compare maps of the effect of SOM removal

library(terra)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

# Load SOC layer

SOC_predicted <- paste0(
  root,
  "/Texture_maps_10m/depth_000_030_cm/Kulstof2022_000_030_cm/",
  "Kulstof_kombineret.tif"
) %>%
  rast()

names(SOC_predicted) <- "SOC"

# Load clay with SOM removed

clay_SOMremoved <- paste0(
  root,
  "/Texture_maps_10m/depth_000_030_cm/Tekstur2024_000_030_cm/",
  "Clay_mean.tif"
  ) %>%
  rast()

names(clay_SOMremoved) <- "clay_SOMremoved"

# Load clay without SOM removal

clay_noSOMremoval <- paste0(
  dir_dat,
  "/results_test_14/bootstrap/predictions_difference/final/depth_000_030_cm/",
  "Clay_norem_mean_000_030_cm.tif"
) %>%
  rast()

names(clay_noSOMremoval) <- "clay_noremoval"

# Extract points

pts <- c(
  clay_SOMremoved,
  clay_noSOMremoval,
  SOC_predicted
) %>%
  spatSample(
    size = 15000,
    na.rm = TRUE
  )

pts_10k <- pts %>%
  filter(
    SOC <= 6
  ) %>%
  mutate(
    clay_dif = clay_SOMremoved - clay_noremoval
  ) %>%
  sample_n(10000)

plot(pts_10k)

my_breaks_dens <- c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1)
mycolors2 <- rgb(
  seq(0, 0, length.out = 20),
  seq(1, 0, length.out = 20),
  seq(1, 1, length.out = 20)
)
mycolors3 <- rgb(
  seq(0, 0, length.out = 21),
  seq(0, 0, length.out = 21),
  seq(1, 0, length.out = 21)
)

mycolorgradient <- c(
  # mycolors1,
  mycolors2[-1],
  mycolors3[-1])

# Clay vs SOC

pts_10k %>%
  pivot_longer(
    cols = c(clay_SOMremoved, clay_noremoval),
    names_to = "method",
    values_to = "clay"
  ) %>%
  ggplot(aes(x = SOC, y = clay)) +
  facet_wrap(~ method) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    ),
    bins = 30
  ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits = c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
  ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# SOM removal effect vs SOC

df_line <- data.frame(
  x = c(1, 2.27, (max(pts_10k$clay_dif) + 12.6294) / 5.88),
  y = c(-0.12, 0.7182, max(pts_10k$clay_dif))
)

# 12.6294
# (max(pts_10k$clay_dif) + 12.6294) / 5.88
# 2.641054
# 2.641054*0.66 - 0.78 + 5.22*(2.641054 - 2.27)

pts_10k %>%
  ggplot(aes(x = SOC, y = clay_dif)) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    ),
    bins = 30
  ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits = c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
  ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(
      min(pts_10k$clay_dif),
      max(pts_10k$clay_dif)
    )
  ) +
  geom_smooth() +
  geom_line(
    aes(x = x, y = y),
    data = df_line,
    inherit.aes = FALSE
  )

# Clay difference vs clay


pts_10k %>%
  ggplot(aes(x = clay_noremoval, y = clay_dif)) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    ),
    bins = 30
  ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits = c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
  ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# Clay difference vs dexter

pts_10k %>%
  ggplot(
    aes(
      x = SOC / clay_SOMremoved,
      y = clay_dif/ clay_SOMremoved
    )
  ) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    ),
    bins = 30
  ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits = c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
  ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_smooth()


pts_10k %>%
  ggplot(
    aes(
      x = SOC / (SOC + clay_SOMremoved),
      y = clay_dif / clay_SOMremoved
    )
  ) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    ),
    bins = 30
  ) +
  scale_fill_gradientn(
    colours = mycolorgradient,
    aesthetics = "fill",
    breaks = my_breaks_dens,
    limits = c(0, 1),
    na.value = 0,
    trans = "sqrt"
  ) +
  guides(
    fill = "legend"
    , alpha = "legend"
  ) +
  scale_alpha_continuous(
    limits= c(0, 1),
    range = c(0.05, 1),
    breaks = my_breaks_dens,
    na.value = 0,
    trans = "sqrt"
  ) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_smooth(method = "lm")


# END