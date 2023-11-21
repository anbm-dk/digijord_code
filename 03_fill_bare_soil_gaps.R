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
terraOptions(tempdir = tmpfolder)

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
  paste0(dir_cov, ., ".tif")

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
  r1 <- rast(ncols = 180, nrows = 180, xmin=0)
  myfilter1 <- focalMat(r1, c(0.5, 1), "Gauss")
  myfilter2 <- focalMat(r1, c(0.5, 1), "Gauss")
  # myfilter2 <- focalMat(r1, c(1, 2), "Gauss")
  
  smoothed_down_list <- list()
  aggregated_list <- list()
  
  split_list <- list()
  merged_list <- list()
  smooth_up_list <- list()
  
  aggregated_list[[1]] <- inrast
  
  smoothed <- focal(
    aggregated_list[[1]],
    w = myfilter1,
    fun = "sum",
    na.policy = "all",
    na.rm = TRUE,
    wopt = list(datatype = datatype(inrast))
  )
  
  summed <- focal(
    !is.na(aggregated_list[[1]]),
    w = myfilter1,
    fun = "sum",
    na.policy = "all",
    na.rm = TRUE
  )
  
  smoothed_down_list[[1]] <- smoothed / summed
  
  for (i in 2:nsteps) {
    aggregated_list[[i]] <- aggregate(
      smoothed_down_list[[i - 1]],
      wopt = list(datatype = datatype(inrast))
      ) 
    
    smoothed <- focal(
      aggregated_list[[i]],
      w = myfilter1,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE,
      wopt = list(datatype = datatype(inrast))
    )
    
    summed <- focal(
      !is.na(aggregated_list[[i]]),
      w = myfilter1,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    
    smoothed_down_list[[i]] <- smoothed / summed
  }
  
  split_list[[nsteps - 1]] <- project(
    smoothed_down_list[[nsteps]],
    aggregated_list[[nsteps - 1]],
    method = "near",
    datatype = datatype(inrast)
  )
  
  merged_list[[nsteps - 1]] <- terra::merge(
    x = smoothed_down_list[[nsteps - 1]],
    y = split_list[[nsteps - 1]],
    wopt = list(datatype = datatype(inrast))
  )
  
  smoothed <- focal(
    merged_list[[nsteps - 1]],
    w = myfilter2,
    fun = "sum",
    na.policy = "all",
    na.rm = TRUE,
    wopt = list(datatype = datatype(inrast))
  )
  
  summed <- focal(
    !is.na(merged_list[[nsteps - 1]]),
    w = myfilter2,
    fun = "sum",
    na.policy = "all",
    na.rm = TRUE
  )
  
  smooth_up_list[[nsteps - 1]] <- smoothed / summed
  
  for (i in (nsteps - 2):1) {
    split_list[[i]] <- project(
      smooth_up_list[[i + 1]],
      aggregated_list[[i]],
      method = "near",
      datatype = datatype(inrast)
    )
    
    merged_list[[i]] <- terra::merge(
      x = smoothed_down_list[[i]],
      y = split_list[[i]],
      wopt = list(datatype = datatype(inrast))
    )
    
    smoothed <- focal(
      merged_list[[i]],
      w = myfilter2,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE,
      wopt = list(datatype = datatype(inrast))
    )
    
    summed <- focal(
      !is.na(merged_list[[i]]),
      w = myfilter2,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    
    smooth_up_list[[i]] <- smoothed / summed
  }
  
  out <- merge(
    inrast,
    smooth_up_list[[1]],
    wopt = list(datatype = datatype(inrast))
  )
  
  return(out)
}

# Fill gaps

for (j in 1:length(names_in)) {
  # Mask all layers (especially s1), to the same extent.
  # Mainly to reduce the effect from edge cells.
  r <- paste0(dir_cov, names_in[[j]], ".tif") %>%
    rast() %>%
    mask(., mask = bare_mask)
  
  r2 <- fill_gaps_gauss(r, 10)
  
  r3 <- mask(
    r2,
    dem,
    filename = paste0(tmpfolder, "filled_", names(r), ".tif"),
    names = paste0("filled_", names(r)),
    datatype = datatype(r)
  )
  
  tmpFiles(remove = TRUE)
}



# END