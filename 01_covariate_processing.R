# Crop and rename covariates

library(terra)
library(magrittr)
library(raster)
library(tools)

root <- getwd()


# 1: Mask to coastline (done)

mycrs <- "EPSG:25832"

dir_cov <- root %>%
  paste0(., "/covariates/")

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

# dem_ind <- grepl(
#   "DHM",
#   cov_files
# )
# 
# dem <- cov_files[dem_ind] %>% rast
# 
# crs(dem) <- mycrs
# 
# not_dem <- cov_files[!dem_ind]
# 
# dir_masked <- root %>%
#   paste0(., "/covariates_masked/") %T>%
#   dir.create()
# 
# not_dem_basenames <- basename(not_dem)
# 
# for (i in 1:length(not_dem)) {
#   dtyp <- not_dem[i] %>%
#     raster() %>%
#     dataType()
#   
#   cov_i <- not_dem[i] %>%
#     terra::rast()
#   
#   crs(cov_i) <- mycrs
#   
#   terra::mask(
#     cov_i,
#     mask = dem,
#     filename = paste0(
#       dir_masked,
#       not_dem_basenames[i]
#     ),
#     overwrite = TRUE
#   )
# }


# 2: Change names (done)
# only underscores
# all lowercase

basenames <- cov_files %>%
  basename() %>%
  file_path_sans_ext() %>%
  gsub('\\.', '_', .) %>%
  tolower()

newnames <- basenames %>%
  paste0(dir_cov, . , ".tif")

# file.rename(
#   cov_files,
#   newnames
# )

# 3: Crop to test area

squareshape <- root %>%
  paste0(., '/testarea_10km/square10km.shp') %>%
  vect

square_ext <- squareshape %>%
  ext %>%
  round(-1)

outfolder <- root %>%
  paste0(., '/testarea_10km/covariates/')

outfolder %>% dir.create

cropstack <- function(
    x,  # list of files
    y,  # extent
    folder # target folder
) {
  for(i in 1:length(x)) {
    r <- x[i] %>% rast
    dtype <- r %>%
      sources %>%
      raster::raster(.) %>%
      raster::dataType(.)
    outname <- r %>%
      sources %>%
      basename %>%
      tools::file_path_sans_ext(.) %>%
      make.names %>%
      paste0(folder, ., '.tif')
    r %>%
      crop(y = y) %>%
      writeRaster(
        datatype = dtype,
        filename = outname,
        overwrite = TRUE
      )
  }
}

cov_files %>%
  cropstack(
    y = square_ext,
    folder = outfolder
  )

# END