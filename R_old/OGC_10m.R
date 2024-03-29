# Oblique geographic coordinates at 10 m resoltion

library(OGC)
library(raster)
library(tidyr)
library(magrittr)
library(terra)

if (!exists("started")) {
  wd <- getwd()
  setwd("..")
  root <- getwd()
  setwd(wd)
  rm(wd)
  started <- TRUE
}

dem_100_withsea_file <- root %>%
  paste0(., "/dem_100m_withsea.tif")

dem_100 <- dem_100_withsea_file %>% raster()

# Make coordinates at 100 m resolution

ogc_100m <- makeOGC(dem_100, 8)

ogc_100_file <- root %>%
  paste0(., "/ogc_100m.tif")

raster::writeRaster(
  ogc_100m,
  filename = ogc_100_file,
  datatype = "INT4S",
  overwrite = TRUE
)

tempfile()

# Load 10 m dem

dem <- "D:/anbm/digijord/cov_basic/DHM2015_terraen_10m.tif" %>%
  rast()

ogc_100m_t <- ogc_100_file %>% rast()

crs(ogc_100m_t) <- crs(dem)

ogc_10_raw_file <- root %>%
  paste0(., "/ogc_10_raw.tif")

ogc_10m_raw <- disagg(ogc_100m_t, 10,
  method = "bilinear",
  filename = ogc_10_raw_file,
  datatype = "INT4S",
  overwrite = TRUE
)

# Crop ogc layers to dem extent

ogc_cropped_file <- root %>%
  paste0(., "/ogc_10m_cropped.tif")

ogc_cropped <- ogc_10m_raw %>%
  crop(.,
    dem,
    filename = ogc_cropped_file,
    overwrite = TRUE
  )

# Mask the 10m ogc values

ogc_masked_file <- root %>%
  paste0(., "/ogc_10m_masked.tif")
#
# ogc_masked <- ogc_cropped %>%
#   mask(.
#        , dem
#        , filename = ogc_masked_file
#        , overwrite = TRUE
#   )

ogc_masked <- ogc_masked_file %>% rast()

ogc_r <- ogc_masked_file %>% raster()

# Write ogc to separate files

ogc_outdir <- root %>%
  paste0(., "/ogc_10m/") %T>%
  dir.create()

ogc_outfiles <- seq(0, 1, length.out = 9) %>%
  head(-1) %>%
  format(digits = 3) %>%
  substr(
    nchar(.) - 2,
    nchar(.)
  ) %>%
  paste0(
    ogc_outdir,
    "ogc_pi0",
    .,
    ".tif"
  )

# ogc_masked %>% writeRaster(
#   filename = ogc_outfiles
#   , overwrite = TRUE
#   , progress = TRUE
#   , datatype = 'INT4S'
#   , filetype = 'GTiff'
# )

ogc_final <- ogc_outfiles %>% rast()

# Create a light version



ogc_light_file <- root %>%
  paste0(., "/ogc_10m_lite.tif")

mins_ogc <- minmax(ogc_masked)[1, ]

lighten <- function(i, mins) {
  out <- base::trunc((i - mins) / 10)
  return(out)
}

# terra::app(
#   ogc_masked
#   , lighten
#   , mins = mins_ogc
#   , cores = 11
#   , filename = ogc_light_file
#   , wopt = list(overwrite = TRUE
#                 , progress = TRUE
#                 , datatype = 'INT2U'
#                 , filetype = 'GTiff')
# )

# Write light version to rasters

ogc_light_outdir <- root %>%
  paste0(., "/ogc_10m_light/") %T>%
  dir.create()

ogc_light_outfiles <- seq(0, 1, length.out = 9) %>%
  head(-1) %>%
  format(digits = 3) %>%
  substr(
    nchar(.) - 2,
    nchar(.)
  ) %>%
  paste0(
    ogc_light_outdir,
    "ogc_lite_pi0",
    .,
    ".tif"
  )

ogc_light <- ogc_light_file %>% rast()

ogc_light %>% writeRaster(
  filename = ogc_light_outfiles,
  overwrite = TRUE,
  progress = TRUE,
  datatype = "INT2U",
  filetype = "GTiff"
)

# END
