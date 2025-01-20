# Map for the test area
# Make maps of the test area for each fraction
# Clip the final maps to the area
# Make maps for mean and SD

library(terra)
library(tidyr)
library(tidyterra)
library(viridisLite)
library(stringr)
library(magrittr)
library(ggplot2)
library(extrafont)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")
dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")
dir_pred_boot <- dir_boot %>%
  paste0(., "/predictions/")

tmpfolder <- paste0(dir_dat, "/Temp/")
terraOptions(tempdir = tmpfolder)

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

fraction_names_underscore <- c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand", "SOC", "CaCO3"
)

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

breaks_titles <- paste0(breaks[1:4], " - ", breaks[2:5], " cm")

squareshape <- dir_dat %>%
  paste0(., "/testarea_10km/square10km.shp") %>%
  vect()

square_ext <- squareshape %>%
  ext() %>%
  round(-1)

dir_test_area_figs <- paste0(dir_results, "/figures_testarea_2025/") %T>%
  dir.create()

layer_spec <- c("mean", "sd")

labels_title_mean <- c(
  expression("Clay"~"(%)"),
  expression("Silt"~"(%)"),
  expression("Fine"~"sand"~"(%)"),
  expression("Coarse"~"sand"~"(%)"),
  expression("SOC"~"(%)"),
  expression("CaCO"[3]~"(%)")
)

labels_title_SD <- c(
  expression("Clay"~"SD"~"(%)"),
  expression("Silt"~"SD"~"(%)"),
  expression("Fine"~"sand"~"SD"~"(%)"),
  expression("Coarse"~"sand"~"SD"~"(%)"),
  expression("SOC"~"SD"~"(%)"),
  expression("CaCO"[3]~"SD"~"(%)")
)

for(i in 1:length(fractions)) {
  for(k in 1:2) {  # layer
    
    layers_cropped <- list()
    
    for (j in 1:4) {   # depth
      breaks_j_chr <- breaks_chr[j:(j + 1)]
      dir_merged_depth <- dir_pred_boot %>%
        paste0(
          ., "/final_maps/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
        )
      
      layer_ij <- dir_merged_depth %>%
        paste0(., fraction_names_underscore[i], "_", layer_spec[k], ".tif") %>%
        rast()
      
      layers_cropped[[j]] <- layer_ij %>%
        terra::crop(x = ., y = squareshape)
    }
    
    layers_cropped %<>%
      rast()
    
    names(layers_cropped) <- breaks_titles
    
    if (k == 1) {
      colorscale_k <- cividis(500)
      labels_title <- labels_title_mean
    } else {
      colorscale_k <- plasma(500)
      labels_title <- labels_title_SD
    }
    
    myplot <- autoplot(layers_cropped) +
      scale_fill_gradientn(
        colors = colorscale_k,
        na.value = NA,
        name = labels_title[i]) +
      theme(text = element_text(family = "serif"))
    
    fname <- paste0(dir_test_area_figs, "/", fraction_names_underscore[i], "_",
                    layer_spec[k], "_10km_test", testn, ".tiff")
    
    tiff(
      fname,
      width = 18,
      height = 14,
      units = "cm",
      res = 600
    )
    
    print(myplot)
    
    try(dev.off())
    try(dev.off())
    
  }
}

# END