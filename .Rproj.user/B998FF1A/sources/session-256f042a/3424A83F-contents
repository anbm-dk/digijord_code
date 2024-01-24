# Fill gaps in the bare soil layers

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(viridis)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

tmpfolder <- paste0(dir_dat, "/Temp/")
terraOptions(
  tempdir = tmpfolder,
  memfrac = 0.3
  )

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates/")

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

bare_count <- paste0(dir_cov, "s2_count_max10.tif") %>% rast()

# Use only cells with at least 10 bare soil observations
reclasser <- c(10, 1) %>% matrix(nrow = 1)

bare_min10 <- terra::classify(
  x = bare_count,
  rcl = reclasser,
  others = NA
)

# Mask by land use and remove field edges

mask_LU <- paste0(dir_dat, "/layers/Mask_LU.tif") %>% rast()
field_edges <- paste0(dir_dat, "/layers/field_edges.tif") %>% rast()

mask_LU_field_edge <- terra::mask(
  x = mask_LU,
  mask = field_edges,
  inverse = TRUE
)

bare_min10_masked <- terra::mask(
  x = bare_min10,
  mask = mask_LU_field_edge
)

# Removing edge cells
bare_mask <- focal(
  bare_min10_masked,
  w = 5,
  fun = "max",
  na.policy = "omit",
  filename = paste0(tmpfolder, "/bare_mask.tif"),
  datatype = "INT2U"
)

# Load covariate info

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

names_in <- cov_cats %>%
  filter(
    category == "bare_soil",
    scorpan == "S"
    ) %>%
  select(name) %>%
  unlist() %>%
  unname()

# Function to fill gaps

fill_gaps_gauss <- function(
    inrast,
    nsteps,
    include_list = FALSE
) {
  r1 <- rast(ncols = 180, nrows = 180, xmin = 0)
  myfilter1 <- round(
    focalMat(r1, c(1, 2), "Gauss"),
    3
  )
  myfilter2 <- myfilter1
  
  smooth_up_list <- list()
  aggregated_list <- list()
  aggregated_list[[1]] <- c(
    inrast*0 + 1,
    inrast
  )
  names(aggregated_list[[1]]) <- c("count", "mean")
  # Stepwise smoothing and aggregation
  for (i in 2:nsteps) {
    smoothed_down <- terra::focal(
      aggregated_list[[i - 1]],
      w = myfilter1,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    aggregated_list[[i]] <- terra::aggregate(
      smoothed_down,  
      fun = "mean",
      na.rm = TRUE
    )
  }
  # Stepwise disaggregation, merging and smoothing
  smooth_up_list[[nsteps]] <- aggregated_list[[nsteps]]
  for (i in (nsteps - 1):1) {
    # Disaggregate by 2
    splitted <- terra::project(
      x = smooth_up_list[[i + 1]],
      y = aggregated_list[[i]],
      method = "near"
    )
    # Merge with aggregated layers
    merged <- terra::merge(
      x = aggregated_list[[i]],
      y = splitted
    )
    # Smoothing
    smooth_up_list[[i]] <- terra::focal(
      merged,
      w = myfilter2,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
  }
  # Divide mean values by the number of cells, to get a weighted mean
  final_lyr <- smooth_up_list[[1]][[2]] / smooth_up_list[[1]][[1]]
  out <- list()
  # Merge with input layer
  out$final <- terra::merge(
    inrast,
    final_lyr,
    wopt = list(datatype = datatype(inrast))
  )
  out$aggregated_list <- aggregated_list
  out$smooth_up_list <- smooth_up_list
  return(out)
}

f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
plot(r)

filled <- fill_gaps_gauss(r, nsteps = 7, include_list = TRUE)

plot(filled$final)

# Fill gaps

for (j in 1:length(names_in)) {
  # Mask all layers (especially s1), to the same extent.
  # Mainly to reduce the effect from edge cells.
  r <- paste0(dir_cov, names_in[[j]], ".tif") %>% rast()
  
  r_masked <- mask(r, mask = bare_mask, datatype = datatype(r))
  
  r2 <- fill_gaps_gauss(
    r_masked,
    11
    )
  
  r3 <- mask(
    r2$final,
    dem,
    filename = paste0(tmpfolder, "filled_", names(r), ".tif"),
    names = paste0("filled_", names(r)),
    datatype = datatype(r),
    gdal = "TILED=YES"
  )
  
  tmpFiles(remove = TRUE)
}

# END