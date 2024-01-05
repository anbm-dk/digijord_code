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

mycrs <- "EPSG:25832"

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

dir_mask_tiles <- dir_dat %>%
  paste0(., "/layers/Mask_LU_tiles/")

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

# Select relevant covariates

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

cov_names_all <- cov_files %>%
  basename() %>%
  tools::file_path_sans_ext(.)

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  select(name) %>%
  unlist() %>%
  unname()

cov_files_selected <- cov_files[cov_names_all %in% cov_selected]

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

source("f_cropstack.R")

# Process for tile creation
# Mask all covariates

cov_files_notogc <- cov_files_selected %>%
  grep('ogc_pi', ., value = TRUE, invert = TRUE)

library(parallel)

numCores <- detectCores()
numCores

showConnections()

cl <- makeCluster(numCores)

clusterEvalQ(
  cl,
  {
    library(terra)
    library(magrittr)
    library(dplyr)
    library(tools)
  }
)

clusterExport(
  cl,
  c(
    "dir_dat",
    "dir_tiles",
    "dir_code",
    "tile_numbers",
    "cov_files_selected",
    "dir_mask_tiles"
  )
)

parSapplyLB(
  cl,
  1:length(tile_shapes),
  function(j) {
    tmpfolder <- paste0(dir_dat, "/Temp/")

    terraOptions(memfrac = 0.02, tempdir = tmpfolder)

    dir_tile_j <- dir_tiles %>%
      paste0(., "/tile_", tile_numbers[j], "/") %T>%
      dir.create()
    
    my_ext <- paste0(
      dir_mask_tiles, "/Mask_LU_tile_", tile_numbers[j], ".tif"
    ) %>%
      rast()

    source(paste0(dir_code, "/f_cropstack.R"))

    cropstack(
      x = cov_files_selected,
      y = my_ext,
      folder = dir_tile_j,
      mask = TRUE
    )
  }
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)


# END
