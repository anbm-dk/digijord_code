# A small map of Denmark

library(terra)
library(tidyterra)
library(ggplot2)
library(tidyverse)
library(magrittr)


dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") 

dem1km <- dir_dat %>%
  paste0(., "/layers/dhm2015_terraen_1km.tif") %>%
  rast()

dem100m <- dir_dat %>%
  paste0(., "/layers/DHM2015_terraen_100m_corrected.tif") %>%
  rast()

europe <- vect("O:/Tech_AGRO/Jord/anbm/GIS/europe_shp/europe_1km.shp")

tiff(
  paste0(dir_results, "/Denmark_map_test_", testn, ".tiff"),
  width = 12,
  height = 9,
  units = "cm",
  res = 600
)

autoplot(dem100m) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    legend.position = c(0.85, 0.725),
    legend.background = element_rect(
      fill="white",
      size = 0.25,
      linetype = "solid", 
      colour ="black"
    )
  ) +
  scale_fill_viridis_c(na.value = NA, name = "Elevation (m)") +
  geom_spatvector(data = europe, fill = NA) +
  coord_sf(xlim = c(441000, 894000), ylim = c(6049000, 6403000))


try(dev.off())


# END