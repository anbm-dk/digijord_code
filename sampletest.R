# Sampling test

library(RSAGA)
library(raster)
library(tidyverse)
library(ClusterR) # Required in function, not to be confused with the clusterR function from raster
library(Rcpp) # Required in function
library(future.apply) # Required in function
library(magrittr) # Required in function
library(fields) # Required in function
library(snow) # Required in function
library(tools) # Required in function


# Empty memory

rm(list = ls())


# Start up

dirname <- dirname(rstudioapi::getSourceEditorContext()$path)
wd <- getwd()
setwd(dirname)
setwd("..")
root <- paste0(getwd())
setwd(wd)
rm(wd, dirname)

source(paste0(root, "/R/sample_kmeans_v2.R"))

paste0(root, "/sampletest/") %>% dir.create()

work_env <- rsaga.env()

dir <- paste0(root, "/sampletest/")

setwd(dir)

rsaga.geoprocessor("grid_calculus",
  module = 6,
  env = work_env,
  param = list(
    TARGET_OUT_GRID = "DEM_random.sgrd",
    ITERATIONS = 1000
  )
)

rsaga.geoprocessor("io_gdal",
  module = 2,
  env = work_env,
  param = list(
    GRIDS = list("DEM_random.sgrd"),
    FILE = "DEM_random.tif"
  )
)

dem <- paste0(root, "/sampletest/DEM_random.tif") %>% raster()

plot(dem, col = cividis(100))

pt <- SpatialPoints(cbind(50, 50))

range01 <- function(x) {
  1 - (x - min(x)) / (max(x) - min(x))
}

dist <- distanceFromPoints(dem, pt) %>% calc(fun = range01)

plot(dist, col = cividis(100))

n_clusters <- 100

rm(sample_kmeans)
source(paste0(root, "/R/sample_kmeans_v2.R"))

set.seed(1)

letstry <- sample_kmeans(
  input = dem,
  use_xy = TRUE,
  clusters = n_clusters,
  scale = TRUE,
  weights = dist
)

cl_colors <- expand.grid(
  r = seq(0, 1, 0.05),
  g = seq(0, 1, 0.05),
  b = seq(0, 1, 0.05)
) %>%
  KMeans_rcpp(clusters = n_clusters) %>%
  .$centroids %>%
  rgb()

plot(letstry$clusters, col = cl_colors)
points(letstry$points)
plot(letstry$distances, col = cividis(100))
points(letstry$points)
plot(dem, col = cividis(100))
points(letstry$points)


# Try it out for Vindum

install.packages("devtools")

library(devtools)

install_bitbucket("abmoeller/ogc/rPackage/OGC")

library(OGC)

data("Vindum_covariates")

# Names sorted as by importance for SOM

sortednames <- c(
  "bluespot", "MRVBF", "SAGAWI", "DEM", "valleydepth", "ECa",
  "slope_gradient", "midslope", "NDVI", "curvature_plan", "SL",
  "SAVI", "aspect_cos", "DVI", "TWI", "RVI", "flow_accu",
  "aspect_sin", "curvature_prof"
)

# Five most important covariates

rm(sample_kmeans)
source(paste0(root, "/R/sample_kmeans_v2.R"))

vin_sel <- subset(Vindum_covariates, sortednames[1:5])

set.seed(1)

letstry2 <- sample_kmeans(
  input = vin_sel,
  use_xy = TRUE,
  pca = TRUE,
  n_pcs = 5,
  clusters = 100,
  scale = TRUE,
  filename_cl = paste0(dir, "vindum_clusters.tif"),
  args_cl = list(
    overwrite = TRUE,
    datatype = "INT2U"
  ),
  filename_d = paste0(dir, "vindum_distances.tif"),
  args_d = list(overwrite = TRUE),
  filename_pts = paste0(dir, "vindum_points.shp"),
  args_pts = list(overwrite = TRUE),
  cores = 2,
  verbose = TRUE,
  sp_pts = FALSE
)

letstry2

plot(letstry2$clusters, col = cl_colors)
points(letstry2$points)
plot(letstry2$distances, col = cividis(100))
points(letstry2$points)
plot(Vindum_covariates$DEM, col = cividis(100))
points(letstry2$points)


# Try for wetland areas

# Source code

source(paste0(root, "/R/loadandstack.R"))


# Load covariates

covs_wl <- paste0(root, "/existing_data/30m_covariates/final_selection/") %>%
  loadandstack()

# Select five important covariates

top5 <- c("IMK_drain_yes_clip", "demdetrend", "valldepth", "elevation", "sagawi")

covs_top5 <- subset(covs_wl, top5)


# Load risk of misclassification

risk <- paste0(root, "/predictions_class_RF_v2/risk.tif") %>% raster()


# Now run the model

rm(sample_kmeans)
source(paste0(root, "/R/sample_kmeans.R"))

dir.create(paste0(root, "/tmp/"))
rasterOptions(tmpdir = paste0(root, "/tmp/"))


# 10,000 points, no weights

time1 <- Sys.time()
print(time1)

set.seed(1)

letstry_wl <- sample_kmeans(
  input = covs_top5,
  ncells = 10^6,
  use_xy = TRUE,
  clusters = 10^4,
  scale = TRUE,
  filename_cl = paste0(dir, "wl_10000_w_clusters.tif"),
  args_cl = list(
    overwrite = TRUE,
    datatype = "INT2U"
  ),
  filename_d = paste0(dir, "wl_10000_w_dist.tif"),
  args_d = list(overwrite = TRUE),
  filename_pts = paste0(dir, "wl_10000_w_points.shp"),
  args_pts = list(overwrite = TRUE),
  cores = 11,
  verbose = TRUE,
  m = 4
  # , weights = risk
)

time2 <- Sys.time()
print(time2)
print(time2 - time1)

endCluster()

# 10,000 points, with weights

time1 <- Sys.time()
print(time1)

set.seed(1)

letstry_wl <- sample_kmeans(
  input = covs_top5,
  ncells = 10^6,
  use_xy = TRUE,
  clusters = 10^4,
  scale = TRUE,
  filename_cl = paste0(dir, "wl_10000_clusters.tif"),
  args_cl = list(
    overwrite = TRUE,
    datatype = "INT2U"
  ),
  filename_d = paste0(dir, "wl_10000_dist.tif"),
  args_d = list(overwrite = TRUE),
  filename_pts = paste0(dir, "wl_10000_points.shp"),
  args_pts = list(overwrite = TRUE),
  cores = 11,
  verbose = TRUE,
  m = 4,
  weights = risk
)

time2 <- Sys.time()
print(time2)
print(time2 - time1)

endCluster()


# Timing:

# 10 points:
# Time difference of 31.67064 mins
#
# letstry_wl <- sample_kmeans(input = covs_top5
#                             , ncells = 10^3
#                             , use_xy = TRUE
#                             , clusters = 10
#                             , scale = TRUE
#                             , filename_cl = paste0(dir, 'wl_10_clusters.tif')
#                             , args_cl = list(overwrite = TRUE
#                                              , datatype = 'INT2U')
#                             , filename_d = paste0(dir, 'wl_10_dist.tif')
#                             , args_d = list(overwrite = TRUE)
#                             , filename_pts = paste0(dir, 'wl_points.shp')
#                             , args_pts = list(overwrite = TRUE)
#                             , cores = 11
#                             , verbose = TRUE
#                             , m = 4
# )


# 100 points:
# Time difference of 34.89091 mins
#
# letstry_wl <- sample_kmeans(input = covs_top5
#                             , ncells = 10^4
#                             , use_xy = TRUE
#                             , clusters = 100
#                             , scale = TRUE
#                             , filename_cl = paste0(dir, 'wl_100_clusters.tif')
#                             , args_cl = list(overwrite = TRUE
#                                              , datatype = 'INT2U')
#                             , filename_d = paste0(dir, 'wl_100_dist.tif')
#                             , args_d = list(overwrite = TRUE)
#                             , filename_pts = paste0(dir, 'wl_100_points.shp')
#                             , args_pts = list(overwrite = TRUE)
#                             , cores = 11
#                             , verbose = TRUE
#                             , m = 4
# )


# 100 points with weights:
# Time difference of 37.67463 mins
#
# letstry_wl <- sample_kmeans(input = covs_top5
#                             , ncells = 10^4
#                             , use_xy = TRUE
#                             , clusters = 100
#                             , scale = TRUE
#                             , filename_cl = paste0(dir, 'wl_100_w_clusters.tif')
#                             , args_cl = list(overwrite = TRUE
#                                              , datatype = 'INT2U')
#                             , filename_d = paste0(dir, 'wl_100_w_dist.tif')
#                             , args_d = list(overwrite = TRUE)
#                             , filename_pts = paste0(dir, 'wl_100_w_points.shp')
#                             , args_pts = list(overwrite = TRUE)
#                             , cores = 11
#                             , verbose = TRUE
#                             , m = 4
#                             , weights = risk
# )


# END
