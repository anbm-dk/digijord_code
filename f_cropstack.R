# Function for cropping covariates

cropstack <- function(
    x, # list of files
    y, # extent
    folder # target folder
    ) {
  for (i in 1:length(x)) {
    require(raster)
    require(terra)
    require(magrittr)
    require(dplyr)
    require(tools)
    r <- x[i] %>%
      terra::rast(.)
    dtype <- x[i] %>%
      raster::raster(.) %>%
      raster::dataType(.)
    outname_base <- x[i] %>%
      base::basename(.) %>%
      tools::file_path_sans_ext(.) %>%
      base::make.names(.)
    outfile <- outname_base %>%
      base::paste0(folder, "/", ., ".tif")
    crop_r <- terra::crop(x = r, y = y)
    names(crop_r) <- outname_base
    crop_r %>%
      terra::writeRaster(
        .,
        datatype = dtype,
        filename = outfile,
        overwrite = TRUE
      )
  }
}

# END
