# 15: Figure to illustrate texture maps

library(terra)
library(tidyterra)
library(magrittr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(viridis)

mycoords_center <- c(562710, 6317250)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

my_ext_figure <- mycoords_center %>%
  t() %>% 
  as.data.frame() %>% 
  as_spatvector(
    geom = c("V1", "V2"), 
    crs = mycrs
    ) %>%
  ext() %>%
  extend(500)

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)

fraction_names_underscore <- c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand", "SOC", "CaCO3"
)

fraction_names_DK <- c(
  "Ler", "Silt", "Finsand", "Grovsand", "Organisk kulstof", "Kalk"
)

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")


dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")


# Directory for saving bootstrap predictions

dir_pred_boot <- dir_boot %>%
  paste0(., "/predictions/")

# Set up loop for predicting each soil depth

breaks <- c(0, 30, 60, 100, 200)

max_char <- breaks %>%
  as.character() %>%
  nchar() %>%
  max()

breaks_chr <- breaks %>%
  str_pad(
    .,
    max_char,
    pad = "0"
  )


# Extract values

map_df_list <- list()
depth_names <- character()

for(i in 1:length(fractions)) {
  maps_i <- list()
  
  for (j in 1:4) {
    breaks_j <- breaks[j:(j + 1)]
    breaks_j_chr <- breaks_chr[j:(j + 1)]
    
    depth_name_j <- paste0(breaks_j[1], " - ", breaks_j[2], " cm")
    depth_names[j] <- depth_name_j
    
    dir_merged_depth <- dir_pred_boot %>%
      paste0(
        ., "/final_maps/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
      )
    
    rasters_j_big <- dir_merged_depth %>%
      paste0(., fraction_names_underscore[i], "_mean.tif") %>%
      rast()
    
    rasters_j_small <- crop(rasters_j_big, my_ext_figure)
    
    maps_i[[j]] <- rasters_j_small
  }
  
  maps_i <- rast(maps_i)
  
  names(maps_i) <- depth_names
  
  # nx <- c(min(minmax(maps_i)), max(minmax(maps_i)))
  # nx[1] <- 0
  
  
  
  # maps_i <- (maps_i - nx[1]) / (nx[2] - nx[1])
  
  map_df_list[[i]] <- maps_i %>%
    as.data.frame(xy = TRUE) %>%
    mutate(
      Fraktion = fraction_names_DK[i],
    ) %>%
    pivot_longer(
      -c(Fraktion, x, y),
      names_to = "Dybde"
    ) %>%
    mutate(
      value = scale(value)
    )
}

map_df_list %<>%
  bind_rows() %>%
  mutate(
    Dybde = factor(Dybde, depth_names, depth_names),
    Fraktion = factor(Fraktion, fraction_names_DK, fraction_names_DK)
  )

# Plotting

tiff(
  paste0(dir_results, "/example_maps_1km.tiff"),
  width = 18,
  height = 12,
  units = "cm",
  res = 300
)

ggplot(
  map_df_list,
  aes(
    x = x,
    y = y,
    fill = value
    )
  ) + 
  geom_raster(show.legend = FALSE) +
  facet_grid(Dybde ~ Fraktion, switch = "y") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),
        panel.spacing.y = unit(0.1, "lines")
  ) +
  coord_equal() +
  scale_fill_viridis(option = "E") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

try(dev.off())

# Tørv2022

kulstof_combo <- paste0(
  root,
  "/Texture_maps_10m/depth_000_030_cm/",
  "Kulstof2022_000_030_cm/Kulstof_arealkombination.tif"
) %>%
  rast() %>%
  crop(my_ext_figure)

peat2022_full <- paste0(
  "O:/Tech_AGRO/Jord/ambe/sinks/feb_2024_full_maps/",
  "cubist_pred_soc_all_final_feb2024.tif"
) %>%
  rast()

peat2022_small <- crop(
  peat2022_full,
  my_ext_figure
  ) %>%
  mask(
    kulstof_combo,
    maskvalue = 1,
    inverse = TRUE
    )

maxvalueforplot <- peat2022_small %>% as.data.frame() %>% max()

tiff(
  paste0(dir_results, "/example_maps_1km_peat2022.tiff"),
  width = 4,
  height = 4,
  units = "cm",
  res = 300
)

plot(
  peat2022_small, 
  range = c(0, maxvalueforplot), 
  col = cividis(255),
  axes = FALSE, 
  box = TRUE, 
  legend = FALSE,
  mar = 0.1
)

try(dev.off())

# Kulstof2022

kulstof2022_full <- paste0(
  root,
  "/Texture_maps_10m/depth_000_030_cm/",
  "Kulstof2022_000_030_cm/Kulstof_kombineret.tif"
) %>%
  rast()

kulstof2022_small <- crop(kulstof2022_full, my_ext_figure)

maxvalueforplot <- kulstof2022_small %>% as.data.frame() %>% max()

tiff(
  paste0(dir_results, "/example_maps_1km_kulstof2022.tiff"),
  width = 4,
  height = 4,
  units = "cm",
  res = 300
)

plot(
  kulstof2022_small, 
  range = c(0, maxvalueforplot), 
  col = cividis(255),
  axes = FALSE, 
  box = TRUE, 
  legend = FALSE,
  mar = 0.1
)

try(dev.off())

# JB2024

clt <- c(
  1,255,240,166,
  2,255,206,181,
  3,252,187,96,
  4,250,168,15,
  5,201,156,87,
  6,161,106,3,
  7,130,85,17,
  8,91,150,54,
  9,76,99,0,
  10,74,74,6,
  11,200,230,80,
  12,130,130,130
) %>%
  matrix(nrow = 4) %>%
  t() %>%
  .[, -1] %>%
  rgb(maxColorValue = 255) %>%
  data.frame(value = 1:12, color = .)

jb2024_full <- paste0(
  root,
  "/Texture_maps_10m/depth_000_030_cm/",
  "/JB2024_000_030_cm/JB_klasse.tif"
) %>%
  rast()

jb2024_small <- crop(jb2024_full, my_ext_figure)

tiff(
  paste0(dir_results, "/example_maps_1km_jb2024.tiff"),
  width = 4,
  height = 4,
  units = "cm",
  res = 300
)

coltab(jb2024_small) <- clt

plot(
  jb2024_small, 
  axes = FALSE, 
  box = TRUE, 
  legend = FALSE,
  mar = 0.1
)

try(dev.off())


# END