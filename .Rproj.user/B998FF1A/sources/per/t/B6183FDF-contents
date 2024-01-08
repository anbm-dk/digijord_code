# 5b: Principal components analysis for covariates

library(terra)
library(caret)
library(magrittr)
library(dplyr)
library(tibble)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")
dir_cov <- dir_dat %>% paste0(., "/covariates")

mycrs <- "EPSG:25832"

# Load covariates

dir_cov <- dir_dat %>% paste0(., "/covariates")

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

cov_files <- dir_cov %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_names %>%
  write.table(
    paste0("cov_names_", Sys.Date(), ".csv")
  )

cov_names[!cov_names %in% cov_cats$name]

cov <- paste0(dir_cov, "/", cov_files) %>%
  rast()

names(cov) <- cov_names

crs(cov) <- mycrs

# Select relevant covariates

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  dplyr::select(., name) %>%
  unlist() %>%
  unname()

cov_use <- terra::subset(cov, cov_selected)

# 2023-12-19: Find covariates with missing values for some islands

# missing_islands <- dir_dat %>%
#   paste0(., "/layers/missing_islands_20231219.shp") %>%
#   vect()
# 
# terra::extract(cov_use, missing_islands) %>%
#   apply(., 2, function(x) sum(is.na(x))) %>%
#   .[. != 0]

# Drop cost_dist
# Fill holes in terodep
 
# cost_dist terodep10m 
# 110         54

# Extract points from tiles in parallel

cells_per_pt <- 2000

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

subdir_tiles <- dir_tiles %>%
  list.dirs() %>%
  .[-1]

library(parallel)

numCores <- detectCores()
numCores

showConnections()

cl <- makeCluster(numCores)

clusterEvalQ(
  cl,
  {
    library(terra)
    library(magrittr)
    library(dplyr)
    library(tools)
  }
)

clusterExport(
  cl,
  c(
    "subdir_tiles",
    "dir_dat",
    "cells_per_pt"
  )
)

pts_tiles <- parSapplyLB(
  cl,
  1:length(subdir_tiles),
  function(j) {
    set.seed(1)
    
    cov_x_files <- subdir_tiles[j] %>%
      list.files(full.names = TRUE)
    
    cov_x_names <- cov_x_files %>%
      basename() %>%
      file_path_sans_ext()
    
    cov_x <- cov_x_files %>% rast()
    
    names(cov_x) <- cov_x_names
    
    cov_pts_x <- terra::spatSample(
      x = cov_x,
      size = ncell(cov_x) / cells_per_pt,
      na.rm = FALSE,
      # as.points = TRUE,
      xy = TRUE,
      exp = 1
    )
    
    return(cov_pts_x)
  },
  simplify = FALSE,
  USE.NAMES = FALSE
)

stopCluster(cl)
foreach::registerDoSEQ()
rm(cl)

cov_pts <- bind_rows(pts_tiles) %>%
  na.omit() %>%
  sample_n(
    min(
      10^5,
      nrow(.)
    )
  )

# Extract random points (sequential)

# set.seed(1)

# cov_pts <- terra::spatSample(
#   x = cov_use,
#   size = 10^5,
#   na.rm = TRUE,
#   as.df = FALSE
# )

# saveRDS(
#   cov_pts,
#   paste0(dir_dat, "cov_pts_pca.rds")
# )

cov_pts <- readRDS(paste0(dir_dat, "cov_pts_pca.rds"))

cov_pts %<>%
  select(-c(x, y)) %>%
  as.matrix()

# Drop OGCs (and existing principal components)

covnames_dropogc <- cov_pts %>%
  colnames() %>%
  grep('ogc_pi', ., value = TRUE, invert = TRUE) %>%
  grep('PC', ., value = TRUE, invert = TRUE)

cov_pts_dropogc <- cov_pts[, colnames(cov_pts) %in% covnames_dropogc]

# Calculate raw principal components

pcs_raw <- prcomp(
  x = cov_pts_dropogc,
  scale. = TRUE
)

# How many PCs does it take to cover 99% of the variation?

summary(pcs_raw)$importance

num_pcs <- sum((cumsum(pcs_raw$sdev) / sum(pcs_raw$sdev)) < 0.99) + 1

# Final principal components

pcs <- prcomp(
  x = cov_pts_dropogc,
  scale. = TRUE,
  rank. = num_pcs
)

pcs %>%
  saveRDS(file = paste0(dir_dat, "pcs_cov.rds"))

# Write rotation to file

pcs_rotation <- pcs$rotation %>%
  as.data.frame() %>%
  rownames_to_column()

pcs_rotation %>%
  saveRDS(file = paste0(dir_dat, "pcs_rotation.rds"))

pcs_rotation %>%
  write.table(
    file = paste0(dir_dat, "pcs_rotation.csv"),
    sep = ";",
    row.names = FALSE
  )

test_pca_10km <- TRUE
# test_pca_10km <- FALSE

if (test_pca_10km) {
  # Load covariates for the test area
  
  dir_cov_10km <- dir_dat %>%
    paste0(., "/testarea_10km/covariates/")
  
  cov_10km <- dir_cov_10km %>%
    list.files(full.names = TRUE) %>%
    rast() %>%
    subset(covnames_dropogc)
  
  # spatSample(cov_10km, 100000) %>%
  #   apply(., 2, function(x) sum(is.na(x))) %>%
  #   .[. != 0]
  
  # Set NAs to zero for terodep10m and the sine and cosine of the aspect
  
  cov_10km$terodep10m %<>% terra::subst(., NA, 0)
  cov_10km$cos_aspect_radians %<>% terra::subst(., NA, 0)
  cov_10km$sin_aspect_radians  %<>% terra::subst(., NA, 0)
  
  # Predict PCs for the test area
  
  pcs_10km  <- terra::predict(cov_10km, pcs, na.rm = TRUE)
  
  library(tidyterra)
  
  tiff(
    paste0(dir_dat, "/cov_pca_10km.tiff"),
    width = 23,
    height = 10,
    units = "cm",
    res = 300
  )
  
  pcs_10km %>%
    terra::subset(1:10) %>%
    autoplot() +
    facet_wrap(~ lyr, ncol = 5)
  
  try(dev.off())
}


# # Predict pcs for the entire country
# 
# library(parallel)
# library(dplyr)
# library(foreach)
# library(stringr)
# 
# dir_cov <- dir_dat %>% paste0(., "/covariates")
# cov_files <- dir_cov %>% list.files()
# cov_names <- cov_files %>% tools::file_path_sans_ext()
# n_digits <- 3
# 
# # Function for predicting pcs
# 
# predict_pcs <- function(mod, dat, n_digits = NULL, ...) {
#   rfun2 <- function(mod2, dat2, n_digits2, ...) {
#     out2 <- predict(
#       object = mod2,
#       newdata = dat2,
#       ...
#     )
#     if (!is.null(n_digits2)) {
#       out2 <- signif(out2, digits = n_digits2)
#     }
#     return(out2)
#   }
#   out <- rfun2(mod, dat, n_digits, ...)
#   return(out)
# }
# 
# # Tiles for model prediction
# 
# numCores <- detectCores()
# numCores
# 
# dir_tiles <- dir_dat %>%
#   paste0(., "/tiles_591/")
# 
# subdir_tiles <- dir_tiles %>%
#   list.dirs() %>%
#   .[-1]
# 
# dir_pcs <- dir_dat %>%
#   paste0(., "/pcs/") %T>%
#   dir.create()
# 
# showConnections()
# 
# cl <- makeCluster(numCores)
# 
# clusterEvalQ(
#   cl,
#   {
#     library(terra)
#     library(caret)
#     library(magrittr)
#     library(dplyr)
#     library(tools)
#     library(stats)
#   }
# )
# 
# clusterExport(
#   cl,
#   c("pcs",
#     "subdir_tiles",
#     "cov_selected",
#     "dir_dat",
#     "n_digits",
#     "dir_pcs",
#     "predict_pcs"
#   )
# )
# 
# parSapplyLB(
#   cl,
#   1:length(subdir_tiles),
#   function(x) {
#     tmpfolder <- paste0(dir_dat, "/Temp/")
#     
#     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
#     
#     cov_x_files <- subdir_tiles[x] %>%
#       list.files(full.names = TRUE)
#     
#     cov_x_names <- cov_x_files %>%
#       basename() %>%
#       file_path_sans_ext()
#     
#     cov_x <- cov_x_files %>% rast()
#     
#     names(cov_x) <- cov_x_names
#     
#     cov_x2 <- subset(cov_x, rownames(pcs$rotation))
#     
#     tilename_x <- basename(subdir_tiles[x])
#     
#     pcs_tilex <- predict(
#       cov_x2,
#       pcs,
#       fun = predict_pcs,
#       na.rm = TRUE,
#       n_digits = n_digits,
#       overwrite = TRUE
#     )
#     
#     for (k in 1:nlyr(pcs_tilex)) {
#       writeRaster(
#         pcs_tilex[[k]],
#         filename = paste0(subdir_tiles[x], "/PC", k, ".tif"),
#         overwrite = TRUE
#       )
#     }
#     
#     write.table(1, file = paste0(dir_pcs,"/", tilename_x, "_done.csv"))
#     
#     return(NULL)
#   }
# )
# 
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)
# 
# for(k in 1:num_pcs) {
#   outtiles_pc_k <- subdir_tiles %>%
#     paste0(., "/PC", k, ".tif") %>%
#     sprc()
#   
#   merge(
#     outtiles_pc_k,
#     filename = paste0(
#       dir_pcs, "/PC", k, ".tif"),
#     overwrite = TRUE,
#     gdal = "TILED=YES",
#     names = paste0(
#       "PC", k
#     )
#   )
# }

# TO do:
# Write a function to convert covariate importance from PCs to original covs
# In this, take into account different sums for the columns.

# END