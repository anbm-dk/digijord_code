# 16: Compare maps of the effect of SOM removal

library(terra)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

dir_fig <- dir_dat %>% paste0(., "/figures/") %T>% dir.create()

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

set.seed(1)

pts <- c(
  clay_SOMremoved,
  clay_noSOMremoval,
  SOC_predicted
) %>%
  spatSample(
    size = 15000,
    na.rm = TRUE,
    xy = TRUE
  )

set.seed(1)

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
      # x = log(SOC / clay_SOMremoved),
      # x = SOC / clay_SOMremoved,
      # x = log(SOC) / clay_SOMremoved,
      # y = clay_dif / clay_SOMremoved
      # y = clay_noremoval / clay_SOMremoved
      y = SOC / clay_SOMremoved - SOC / clay_noremoval
      # y = log(clay_noremoval) - log(clay_SOMremoved)
    )
  ) +
  geom_hex(
    aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    )
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

# Load raster with 2014 clay contents

texture_2014_clay <- paste0(
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/",
  "geotiffs/aclaynor.tif"
) %>%
  rast()

pts_15k_v <- vect(
  pts,
  geom = c("x", "y"),
  crs = mycrs
)

pts_15k_2014extr <- terra::extract(
  texture_2014_clay,
  pts_15k_v
)

# Plot distributions

library(ggridges)

clay_breaks2 <- c(5, 10, 15, 25, 45)

set.seed(1)

df_3maps <- pts_15k_v %>%
  values() %>%
  bind_cols(pts_15k_2014extr) %>%
  mutate(
    clay_2014 = aclaynor
  ) %>%
  select(
    clay_SOMremoved, clay_noremoval, clay_2014
  ) %>%
  drop_na() %>%
  sample_n(10000) %>%
  pivot_longer(
    cols = c(clay_SOMremoved, clay_noremoval, clay_2014)
  ) 

tiff(
  paste0(dir_fig, "/Clay_2014_2024_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
df_3maps %>%
  ggplot(aes(x = value,
             y = name,
             fill = stat(x))) +
  geom_density_ridges_gradient(
    rel_min_height = 0.001,
    scale = 0.95,
    panel_scaling = FALSE,
    bandwidth = 1
  ) +
  scale_fill_viridis_b(breaks = clay_breaks2) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1)), limits = rev)

try(dev.off())

library(modeest)

df_3maps %>%
  group_by(name) %>% summarise(
    mode = mlv(value, method = "meanshift")
  )

# Modes
# 1 clay_2014        3.60
# 3 clay_noremoval   4.56
# 2 clay_SOMremoved  6.41


# END