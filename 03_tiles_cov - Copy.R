# 09: Split covariates using tiles

# 1: Start up

library(terra)
library(magrittr)
library(dplyr)
library(stringr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")
dir_cov <- dir_dat %>% paste0(., "/covariates")

testn <- 6
mycrs <- "EPSG:25832"

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")


# Load tile shape polygons

tile_shapes <- dir_tiles %>%
  paste0(., "/tiles.shp") %>%
  vect()

# Find covariates

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

# Make names for tiles

max_char <- length(tile_shapes) %>%
  1:. %>%
  as.character() %>%
  nchar() %>%
  max()

tile_numbers <- length(tile_shapes) %>%
  1:. %>%
  str_pad(
    .,
    max_char,
    pad = "0"
  )

# Cropping function

cropstack <- function(
    x, # list of files
    y, # extent
    folder # target folder
    ) {
  for (i in 1:length(x)) {
    r <- x[i] %>% rast()
    dtype <- r %>%
      sources() %>%
      raster::raster(.) %>%
      raster::dataType(.)
    outname <- r %>%
      sources() %>%
      basename() %>%
      tools::file_path_sans_ext(.) %>%
      make.names() %>%
      paste0(folder, ., ".tif")
    names(r) <- r %>%
      sources() %>%
      basename() %>%
      tools::file_path_sans_ext(.) %>%
      make.names()
    r %>%
      crop(y = y) %>%
      writeRaster(
        datatype = dtype,
        filename = outname,
        overwrite = TRUE
      )
  }
}

# Loop for tile creation

for (i in 1:length(tile_shapes)) {
  print(i)

  dir_tile_i <- dir_tiles %>%
    paste0(., "/tile_", tile_numbers[i], "/") %T>%
    dir.create()

  cov_files %>%
    cropstack(
      y = tile_shapes[i],
      folder = dir_tile_i
    )
}

# END
