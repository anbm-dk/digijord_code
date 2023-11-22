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
bare_mask <- "s2_geomedian_b2" %>%
  paste0(dir_cov, ., ".tif") %>%
  rast()

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
    nsteps
) {
  r1 <- rast(ncols = 180, nrows = 180, xmin = 0)
  myfilter1 <- focalMat(r1, c(1, 2), "Gauss")
  myfilter2 <- focalMat(r1, c(1, 2), "Gauss")
  
  smooth_up_list <- list()
  aggregated_list <- list()
  aggregated_list[[1]] <- list()
  aggregated_list[[1]][[1]] <- inrast
  aggregated_list[[1]][[2]] <- !is.na(inrast)
  
  # Function for weighted aggregation
  agg_weight <- function(x, x_count) {
    prod_w <- x*x_count
    prod_agg <- terra::aggregate(
      prod_w,
      fun = "sum",
      na.rm = TRUE
      )
    count_agg <- terra::aggregate(
      x_count,
      fun = "sum",
      na.rm = TRUE
      )
    out2 <- list()
    out2$mean <- prod_agg / count_agg
    out2$count <- count_agg
    return(out2)
  }
  # Function for weighted smoothing
  smooth_weight <- function(x, x_count, filt) {
    prod_w <- x*x_count
    smooth_prod <- focal(
      prod_w,
      w = filt,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    smooth_count <- focal(
      x_count,
      w = filt,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    out2 <- list()
    out2$mean <- smooth_prod / smooth_count
    out2$count <- smooth_count
    return(out2)
  }
  # Function for projecting maps and weights
  project_weight <- function(x, x_count, targ, dtyp) {
    out2 <- list()
    out2$x <- terra::project(
      x = x,
      y = targ,
      method = "near",
      datatype = dtyp
    )
    out2$x_count <- terra::project(
      x = x_count / 4,  # NB
      y = targ,
      method = "near"
    )
    return(out2)
  }
  # Function for merging maps and weights
  merge_weight <- function(x, y, dtyp) {
    out2 <- list()
    out2$x <- terra::merge(
      x = x[[1]],
      y = y[[1]],
      wopt = list(datatype = dtyp)
    )
    out2$x_count <- terra::ifel(
      is.na(x[[1]]),
      yes = y[[2]],
      no = x[[2]]
    )
    return(out2)
  }
  # Stepwise aggregation
  for (i in 2:nsteps) {
    smoothed_down <- smooth_weight(
      x = aggregated_list[[i - 1]][[1]],
      x_count = aggregated_list[[i - 1]][[2]],
      myfilter1
    )
    aggregated_list[[i]] <- agg_weight(
      x = smoothed_down[[1]],
      x_count = smoothed_down[[2]]
    )
  }
  smooth_up_list[[nsteps]] <- aggregated_list[[nsteps]]
  # Stepwise disaggregation
  for (i in (nsteps - 1):1) {
    # Disaggregate by 2
    splitted <- project_weight(
      x = smooth_up_list[[i + 1]][[1]],
      x_count = smooth_up_list[[i + 1]][[2]],
      targ = aggregated_list[[i]][[1]],
      dtyp = datatype(inrast)
    )
    # Merge means and counts
    merged <- merge_weight(
      x = aggregated_list[[i]],
      y = splitted,
      dtyp = datatype(inrast)
    )
    # Weighted smoothing
    smooth_up_list[[i]] <- smooth_weight(
      x = merged[[1]],
      x_count = merged[[2]],
      filt = myfilter2
    )
  }
  final_lyr <- smooth_up_list[[1]][[1]]
  out <- list()
  out$final <- terra::merge(
    inrast,
    final_lyr,
    wopt = list(datatype = datatype(inrast))
  )
  out$aggregated_list <- aggregated_list
  out$smooth_up_list <- smooth_up_list
  out$splitted <- splitted
  out$merged <- merged
  return(out)
}


f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
plot(r)

filled <- fill_gaps_gauss(r, 4)

plot(filled$final)

# Fill gaps

for (j in 1:length(names_in)) {
  # Mask all layers (especially s1), to the same extent.
  # Mainly to reduce the effect from edge cells.
  r <- paste0(dir_cov, names_in[[j]], ".tif") %>% rast()
  
  r_masked <- mask(r, mask = bare_mask, datatype = datatype(r))
  
  r2 <- fill_gaps_gauss(
    r_masked,
    10,
    smooth_down = FALSE
    )
  
  r3 <- mask(
    r2,
    dem,
    filename = paste0(tmpfolder, "filled_", names(r), ".tif"),
    names = paste0("filled_", names(r)),
    datatype = datatype(r),
    gdal = "TILED=YES"
  )
  
  tmpFiles(remove = TRUE)
}



# END