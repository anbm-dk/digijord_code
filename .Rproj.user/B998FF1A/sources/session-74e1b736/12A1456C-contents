# Mask for predictions based on LU basemap

library(terra)
library(magrittr)
library(dplyr)
library(stringr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")
dir_cov <- dir_dat %>% paste0(., "/covariates")

# Load raster data

mycrs <- "EPSG:25832"

dem <- dir_cov %>%
  paste0("/dhm2015_terraen_10m.tif") %>%
  rast()

basemap <- dir_dat %>%
  paste0(., "/layers/Basemap/Geotiff/Basemap04_2021/lu_agg_2021.tif") %>%
  rast()

crs(dem) <- mycrs
crs(basemap) <- mycrs

# Table for reclassification

rcl_table_full <- dir_code %>%
  paste0(., "/LU_basemap_mask.csv") %>%
  read.table(
    ., header = TRUE, sep = ";"
  ) %>% mutate(
    Mask = case_when(
      Mask == 1 ~ 1,
      .default = NA
    )
  )

rcl_mat <- rcl_table_full %>%
  select(c(C_01, Mask)) %>%
  as.matrix()

# Reclassfiffy basemap to produce mask

tmpfolder <- paste0(dir_dat, "/Temp/")

terraOptions(tempdir = tmpfolder)

mask_full <- classify(
  basemap,
  rcl = rcl_mat,
  others = NA
)

names(mask_full) <- "Mask_LU"

# Crop mask and write to file

mask_crop <- crop(
  mask_full,
  dem,
  extend = TRUE,
  filename = paste0(dir_dat, "/layers/Mask_LU.tif"),
  datatype = "INT1U",
  overwrite = TRUE
)

# Split mask into tiles

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

tile_shapes <- dir_tiles %>%
  paste0(., "/tiles.shp") %>%
  vect()

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

dir_mask_tiles <- dir_dat %>%
  paste0(., "/layers/Mask_LU_tiles/") %T>%
  dir.create()

mask_file <- sources(mask_crop)

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
    "tile_numbers",
    "dir_mask_tiles",
    "mask_file"
  )
)

parSapplyLB(
  cl,
  1:length(tile_shapes),
  function(j) {
    tmpfolder <- paste0(dir_dat, "/Temp/")
    
    terraOptions(memfrac = 0.02, tempdir = tmpfolder)
    
    tile_shapes <- dir_tiles %>%
      base::paste0(., "/tiles.shp") %>%
      terra::vect()
    
    mask_r <- rast(mask_file)

    crop(
      mask_r,
      y = tile_shapes[j],
      filename = paste0(
        dir_mask_tiles, "/Mask_LU_tile_", tile_numbers[j], ".tif"
        ),
      datatype = "INT1U",
      overwrite = TRUE
    )
  }
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)


# End