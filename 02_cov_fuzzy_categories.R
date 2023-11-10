# 02: Fuzzy boundaries for categorical covariates

# Process these layers:
# wetlands_10m.tif - 1:20,000   - smallest units about 15 m across [ok]
# geology_         - 1:25,000   - smallest units about 25 m across [ok]
# landscape_       - 1:100,000  - smallest units about 100 m across [ok]
# georeg_          - 1:100,000  - the uncertainty seems to reach 500 m in some cases [ok]
# lu_              - 10 m - (Corine LU has a scale of 1:00,000, but the basemap has 10 m resolution)
# Use a sigma of less than 1 for lu.  [ok]
# imk             - 10 m, use half sigma for fuzzification [ok]
# cwl_10m_  # Already processed in ArcGIS (original resolution 20 m) [ok]

# 1: Start up

library(terra)
library(magrittr)
library(tools)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates/")

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

tmpfolder <- paste0(dir_dat, "/Temp/")

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

# Process wetlands layer [ok]

wetland_ind <- grepl(
  "wetlands_10m.tif",
  cov_files
)

wetlands_crisp <- cov_files[wetland_ind] %>% rast()
# drylands_crisp <- 1 - wetlands_crisp
# 
# wl_twolayers_crisp <- c(wetlands_crisp, drylands_crisp)

my_focal_weights <- focalMat(
  wetlands_crisp,
  c(10, 20),
  type = c('Gauss')
)
# 
# wl_twolayers_fuzzy <- focal(
#   wl_twolayers_crisp,
#   w = my_focal_weights,
#   na.policy = "omit",
#   na.rm = TRUE
# )
# 
# wl_twolayers_fuzzy_sum <- sum(wl_twolayers_fuzzy)
# wl_fuzzy_norm <- wl_twolayers_fuzzy[[1]] / wl_twolayers_fuzzy_sum
# wl_fuzzy_norm_round <- round(wl_fuzzy_norm, digits = 2)
# names(wl_fuzzy_norm_round) <- "wetlands_10m_fuzzy"
# 
# writeRaster(
#   wl_fuzzy_norm_round,
#   filename = paste0(tmpfolder, "/wetlands_10m_fuzzy.tif"),
#   datatype = "FLT4S",
#   overwrite = TRUE,
#   gdal = "TILED=YES"
#   )

# Process geological map [ok]

# geology_ind <- grepl(
#   "geology",
#   cov_files
# )
# 
# geology_crisp <- cov_files[geology_ind] %>% rast()
# 
# geology_sum <- sum(geology_crisp)
# geology_res <- 1 - geology_sum
# names(geology_res) <- "geology_res"
# geology_crisp_full <- c(geology_crisp, geology_res)
# 
# geology_fuzzy <- focal(
#   geology_crisp_full,
#   w = my_focal_weights,
#   na.policy = "omit",
#   na.rm = TRUE
# )
# 
# geology_fuzzy_sum <- sum(geology_fuzzy)
# geology_fuzzy_norm <- geology_fuzzy / geology_fuzzy_sum
# geology_fuzzy_norm_round <- round(geology_fuzzy_norm, digits = 2)
# geology_names <- names(geology_fuzzy_norm_round)
# geology_names_fuzzy <- paste0("fuzzy_", geology_names)
# geology_files_fuzzy <- paste0(tmpfolder, geology_names_fuzzy, ".tif")
# names(geology_fuzzy_norm_round) <- geology_names_fuzzy
# 
# for (i in 1:nlyr(geology_fuzzy_norm_round)) {
#   writeRaster(
#     geology_fuzzy_norm_round[[i]],
#     filename = geology_files_fuzzy[[i]],
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# Function for remaining layers

fuzzify_indicators <- function(
  x = NULL,
  aggregation_factor = 1,
  residual_layer = TRUE,
  local_filter = NULL,
  n_digits = 3,
  outfolder = NULL
) {
  x_sum <- sum(x)
  
  if (residual_layer) {
    x_res <- 1 - x_sum
    names(x_res) <- "x_res"
    x_crisp_full <- c(x, x_res)
  } else {
    x_crisp_full <- x
  }
  
  if (aggregation_factor > 1) {
    x_crisp_full <- terra::aggregate(
      x_crisp_full,
      fact = aggregation_factor,
      fun = "sum"
    )
    
    x_fuzzy <- focal(
      x_crisp_full,
      w = local_filter,
      na.policy = "all",
      na.rm = TRUE
    )
    
    x_fuzzy <- terra::resample(
      x = x_fuzzy,
      y = x,
      method = "cubicspline"
    )
    
    x_fuzzy <- mask(
      x_fuzzy,
      mask = x[[1]]
    )
  } else {
    x_fuzzy <- focal(
      x_crisp_full,
      w = local_filter,
      na.policy = "omit",
      na.rm = TRUE
    )
  }

  x_fuzzy_sum <- sum(x_fuzzy)
  x_fuzzy_norm <- x_fuzzy / x_fuzzy_sum
  x_fuzzy_norm_round <- signif(x_fuzzy_norm, digits = n_digits)
  x_names <- names(x_fuzzy_norm_round)
  x_names_fuzzy <- paste0("fuzzy_", x_names)
  x_files_fuzzy <- paste0(outfolder, x_names_fuzzy, ".tif")
  names(x_fuzzy_norm_round) <- x_names_fuzzy

  for (i in 1:nlyr(x_fuzzy_norm_round)) {
    writeRaster(
      x_fuzzy_norm_round[[i]],
      filename = x_files_fuzzy[[i]],
      datatype = "FLT4S",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
  }
  
  invisible(NULL)
}

# Process landscape elements [ok]

# landscape_ind <- grepl(
#   "landscape",
#   cov_files
# )
# 
# landscape_crisp <- cov_files[landscape_ind] %>% rast()
# 
# fuzzify_indicators(
#   landscape_crisp,
#   aggregation_factor = 5,
#   local_filter = my_focal_weights,
#   n_digits = 2,
#   outfolder = tmpfolder
# )

# Process georegions

# georeg_ind <- grepl(
#   "georeg_",
#   cov_files
# )
# 
# georeg_crisp <- cov_files[georeg_ind] %>% rast()
# 
# fuzzify_indicators(
#   georeg_crisp,
#   aggregation_factor = 5,
#   local_filter = my_focal_weights,
#   n_digits = 2,
#   outfolder = tmpfolder
# )

# Sigma experiment

r1 <- matrix(
  sample(c(0,1), 400, replace = TRUE), nrow = 20
) %>%
  rast()

plot(r1)

r2 <- focal(r1, w = my_focal_weights, na.rm = TRUE)

plot(r2)

# halfsigma <- focalMat(r1, d = c(0.56, 2), type = 'Gauss')

halfsigma <- focalMat(r1, d = c(0.5, 1), type = 'Gauss')

halfsigma

r3 <- focal(r1, w = halfsigma, na.rm = TRUE)

library(viridis)

plot(r1, col = cividis(100))
plot(r3, col = cividis(100))

hist(r3, breaks = 100)

plot(r3 - r1)

r3 - r1

# Process LU [ok]

# lu_ind <- grepl(
#   "lu_",
#   cov_files
# )
# 
# lu_crisp <- cov_files[lu_ind] %>% rast()
# 
# fuzzify_indicators(
#   lu_crisp,
#   local_filter = halfsigma,
#   n_digits = 2,
#   outfolder = tmpfolder
# )

# Process imk [ok]

# imk_files <- grep(
#   "imk_",
#   cov_files,
#   value = TRUE
# )
# 
# for (i in 1:length(imk_files)) {
#   r <- rast(imk_files[i])
#   
#   newname <- names(r) %>% paste0("fuzzy_", .)
#   
#   r_fuzzy <- x_fuzzy <- focal(
#     r,
#     w = halfsigma,
#     na.policy = "omit",
#     na.rm = TRUE
#   )
#   
#   r_fuzzy_round <- round(
#     r_fuzzy,
#     digits = 2
#   )
#   
#   names(r_fuzzy_round) <- newname
#   
#   writeRaster(
#     r_fuzzy_round,
#     filename = paste0(tmpfolder, "/", newname, ".tif"),
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# Make bare soil count layer with a maximum value of 10
# Fuzzify it with a half sigma

s2_count <- dir_dat %>%
  paste0(
    ., "layers/s2_geomedian_count.tif"
  ) %>% rast()

names(s2_count) <- "s2_geomedian_count"

s2_count_masked <- mask(
  s2_count,
  dem
)

s2_count_max10 <- ifel(
  s2_count_masked > 10,
  10,
  s2_count_masked
)

names(s2_count_max10) <- "s2_count_max10"

writeRaster(
  s2_count_max10,
  filename = paste0(tmpfolder, "/s2_count_max10.tif"),
  datatype = "INT1U",
  overwrite = TRUE,
  gdal = "TILED=YES"
)
  
halfsigma_round <- round(halfsigma, digits = 2)

s2_count_max10_fuzzy <- focal(
  s2_count_max10,
  w = halfsigma_round,
  na.policy = "omit"
)

names(s2_count_max10_fuzzy) <- "s2_count_max10_fuzzy"

writeRaster(
  s2_count_max10_fuzzy,
  filename = paste0(tmpfolder, "/s2_count_max10_fuzzy.tif"),
  datatype = "FLT4S",
  overwrite = TRUE,
  gdal = "TILED=YES"
)

# END