# Sampling test

library(terra)
library(tidyverse)
library(ClusterR) # Required in function, not to be confused with the clusterR function from raster
library(Rcpp) # Required in function
library(future.apply) # Required in function
library(magrittr) # Required in function
library(fields) # Required in function
library(tools) # Required in function
library(parallel)

# Empty memory

rm(list = ls())

# Start up

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

mycrs <- "EPSG:25832"

try(rm(sample_kmeans))
source("f_sample_kmeans.R")

library(obliquer)

myvars <- unwrap(Vindum_covariates)

crs(myvars) <- mycrs

safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)

# Try it out for Vindum

set.seed(1)

letstry2 <- sample_kmeans(
  input = myvars,
  use_xy = TRUE,
  # pca = TRUE,
  # n_pcs = 5,
  clusters = 10,
  scale = TRUE,
  filename_cl = paste0(dir_dat, "vindum_clusters.tif"),
  args_cl = list(
    overwrite = TRUE,
    datatype = "INT2U"
  ),
  filename_d = paste0(dir_dat, "vindum_distances.tif"),
  args_d = list(overwrite = TRUE),
  filename_pts = paste0(dir_dat, "vindum_points.shp"),
  args_pts = list(overwrite = TRUE),
  cores = 2,
  verbose = TRUE,
  sp_pts = FALSE
)

letstry2

letstry2$clusters %>% as.factor() %>% plot(col = safe_colorblind_palette)
points(letstry2$points, pch = 21, bg = "white")
plot(letstry2$distances, col = cividis(100))
points(letstry2$points)



# END
