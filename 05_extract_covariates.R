# 03: Extract covariates

library(terra)
library(magrittr)
library(exactextractr)
library(sf)

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

# 2 Load covariates

cov_dir <- dir_dat %>% paste0(., "/covariates")

cov_files <- cov_dir %>% list.files()

cov_names <- cov_files %>% tools::file_path_sans_ext()

cov <- paste0(cov_dir, "/", cov_files) %>%
  rast()

names(cov) <- cov_names

crs(cov) <- mycrs

# 3 Create buffers (40 m = ~ 0.5 ha)

buffer_dsc <- terra::buffer(
  dsc,
  width = 40
) %>%
  st_as_sf()

buffer_SEGES <- terra::buffer(
  SEGES,
  width = 40
) %>%
  st_as_sf()

# 4 Extract

dsc_extr <- terra::extract(
  x = cov,
  y = dsc,
  ID = FALSE,
)

# buffer_dsc_extr <- exact_extract(
#     x = cov,
#     y = buffer_dsc,
#     fun = "mean",
#     progress = TRUE
#   )

# names(buffer_dsc_extr) <- names(cov)

SEGES_extr <- terra::extract(
    x = cov,
    y = SEGES,
    ID = FALSE,
  )

# buffer_SEGES_extr <- exact_extract(
#   x = cov,
#   y = buffer_SEGES,
#   fun = "mean",
#   progress = TRUE
# )

# names(buffer_SEGES_extr) <- names(cov)

SINKS_extr <- terra::extract(
    x = cov,
    y = SINKS,
    ID = FALSE,
  )

profiles_extr <- terra::extract(
  x = cov,
  y = profiles_shp,
  ID = FALSE,
)

profiles_extr$PROFILNR <- profiles_shp$PROFILNR

forests_extr <- terra::extract(
  x = cov,
  y = forest_samples,
  ID = FALSE,
)

# 5 Write to csv

dir_extr <- dir_dat %>%
  paste0(., "/extracts/")

write.table(
  dsc_extr,
  paste0(dir_extr, "dsc_extr.csv"),
  row.names = FALSE,
  sep = ";"
)

# write.table(
#   buffer_dsc_extr,
#   paste0(dir_extr, "buffer_dsc_extr.csv"),
#   row.names = FALSE,
#   sep = ";"
# )

write.table(
  SEGES_extr,
  paste0(dir_extr, "SEGES_extr.csv"),
  row.names = FALSE,
  sep = ";"
)

# write.table(
#   buffer_SEGES_extr,
#   paste0(dir_extr, "buffer_SEGES_extr.csv"),
#   row.names = FALSE,
#   sep = ";"
# )

write.table(
  SINKS_extr,
  paste0(dir_extr, "SINKS_extr.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  profiles_extr,
  paste0(dir_extr, "profiles_extr.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  forests_extr,
  paste0(dir_extr, "forests_extr.csv"),
  row.names = FALSE,
  sep = ";"
)

# Save as RDS

saveRDS(
  dsc_extr,
  paste0(dir_extr, "dsc_extr.rds")
)

saveRDS(
  SEGES_extr,
  paste0(dir_extr, "SEGES_extr.rds")
)

saveRDS(
  SINKS_extr,
  paste0(dir_extr, "SINKS_extr.rds")
)

saveRDS(
  profiles_extr,
  paste0(dir_extr, "profiles_extr.rds")
)

saveRDS(
  forests_extr,
  paste0(dir_extr, "forests_extr.rds")
)

# END