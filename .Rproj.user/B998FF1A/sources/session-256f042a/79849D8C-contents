# 04: Create folds

library(terra)
library(magrittr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

mycrs <- "EPSG:25832"

# 1 Load observations

dir_obs_proc <- dir_dat %>%
  paste0(., "/observations/processed/")

dsc <- dir_obs_proc %>%
  paste0(., "dsc.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

SEGES <- dir_obs_proc %>%
  paste0(., "SEGES.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

SINKS <- dir_obs_proc %>%
  paste0(., "SINKS.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

profiles_shp <- dir_dat %>%
  paste0(
    .,
    "/observations/profiles/Profiles_coordinates_new/Profiles_coordinates_new.shp"
  ) %>%
  vect()

forest_samples <- dir_obs_proc %>%
  paste0(., "forest_samples.csv") %>%
  read.table(
    header = TRUE,
    sep = ";",
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

# 2 Load and aggregate DEM

cov_dir <- dir_dat %>% paste0(., "/covariates/")

DEM_10m <- cov_dir %>%
  paste0(., "/dhm2015_terraen_10m.tif") %>%
  rast()

fun_agg <- function(x) {
  out <- sum(!is.na(x)) == 0
  return(out)
}

dir_folds <- dir_dat %>%
  paste0(., "/folds/") %T>%
  dir.create()

# file_dem_mask_100m <- paste0(dir_folds, "/dem_mask_100m.tif")
#
# terra::aggregate(
#   DEM_10m,
#   10,
#   fun = fun_agg,
#   cores = 19,
#   filename = file_dem_mask_100m
# )
#
# dem_mask_100m <- rast(file_dem_mask_100m)
#
# dem_mask_100m2 <- ifel(dem_mask_100m == 1, NA, 1)

file_dem_mask_100m_2 <- paste0(dir_folds, "/dem_mask_100m_2.tif")

# writeRaster(
#   dem_mask_100m2,
#   filename = file_dem_mask_100m_2,
#   overwrite = TRUE,
#   datatype = "INT1U"
#   )

dem_mask_100m2 <- rast(file_dem_mask_100m_2)

plot(dem_mask_100m2)

# 3: Calculate 10 folds

set.seed(1)

# folds_10_100m <- app(
#   dem_mask_100m2,
#   function(x) {
#     out <- x*0 + sample.int(10, 1)
#     return(out)
#   },
#   cores = 2
# )

file_folds_10_100m <- paste0(dir_folds, "/folds_10_100m.tif")

# writeRaster(
#   folds_10_100m,
#   filename = file_folds_10_100m,
#   overwrite = TRUE,
#   datatype = "INT1U"
#   )

folds_10_100m <- rast(file_folds_10_100m)

# 4 Extract

dsc_folds <- terra::extract(
  x = folds_10_100m,
  y = dsc,
  ID = FALSE,
)

SEGES_folds <- terra::extract(
  x = folds_10_100m,
  y = SEGES,
  ID = FALSE,
)

SINKS_folds <- terra::extract(
  x = folds_10_100m,
  y = SINKS,
  ID = FALSE,
)

profiles_folds <- terra::extract(
  x = folds_10_100m,
  y = profiles_shp,
  ID = FALSE,
)

profiles_folds$PROFILNR <- profiles_shp$PROFILNR

forests_folds <- terra::extract(
  x = folds_10_100m,
  y = forest_samples,
  ID = FALSE,
)

# 5 Write to file

write.table(
  dsc_folds,
  paste0(dir_folds, "/dsc_folds.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  SEGES_folds,
  paste0(dir_folds, "/SEGES_folds.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  SINKS_folds,
  paste0(dir_folds, "/SINKS_folds.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  profiles_folds,
  paste0(dir_folds, "/profiles_folds.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  forests_folds,
  paste0(dir_folds, "/forest_folds.csv"),
  row.names = FALSE,
  sep = ";"
)


# Raster layers for poisson bootstrap

fun1 <- function(x, n_layers = 1) {
  # if (is.matrix(x) & ncol(x) > 1) {
  #   x <- x[, 1]
  # }
  out <- matrix(numeric(), nrow = length(x), ncol = n_layers)
  sumNA <- sum(is.na(x))
  out[!is.na(x), ] <- 0 + rpois(n_layers*(length(x) - sumNA), 1)
  # out <- x * 0 + rpois(n_layers, 1)
  # out <- as.matrix(out, nrow = length(x))
  return(out)
}

nlyr_out <- 100

# set.seed(8863)
# 
# r_poisson <- app(
#   dem_mask_100m2,
#   fun = function(i, ff, outlayers) ff(i, n_layers = outlayers),
#   cores = 19,
#   ff = fun1,
#   outlayers = nlyr_out,
#   filename  = paste0(dir_folds, "/poisson_100m.tif"),
#   overwrite = TRUE,
#   wopt = list(
#     datatype = "INT1U"
#   )
# )

poisson_r <- paste0(dir_folds, "/poisson_100m.tif") %>% rast()

poisson_r

plot(poisson_r[[1]])

# END
