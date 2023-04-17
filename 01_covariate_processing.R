# 01: Crop and rename covariates

# 1: Start up

library(terra)
library(magrittr)
library(raster)
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

dem <- cov_files[dem_ind] %>% rast

crs(dem) <- mycrs

# 1: Mask to coastline (done)

# not_dem <- cov_files[!dem_ind]
# 
# dir_masked <- dir_dat %>%
#   paste0(., "/covariates_masked/") %T>%
#   dir.create()
# 
# not_dem_basenames <- basename(not_dem)
# 
# for (i in 1:length(not_dem)) {
#   dtyp <- not_dem[i] %>%
#     raster() %>%
#     dataType()
#   
#   cov_i <- not_dem[i] %>%
#     terra::rast()
#   
#   crs(cov_i) <- mycrs
#   
#   terra::mask(
#     cov_i,
#     mask = dem,
#     filename = paste0(
#       dir_masked,
#       not_dem_basenames[i]
#     ),
#     overwrite = TRUE
#   )
# }


# 2: Change names (done)
# only underscores
# all lowercase

# basenames <- cov_files %>%
#   basename() %>%
#   file_path_sans_ext() %>%
#   gsub("\\.", "_", .) %>%
#   gsub("-", "_", .) %>%
#   tolower()
# 
# newnames <- basenames %>%
#   paste0(dir_cov, . , ".tif")
# 
# file.rename(
#   cov_files,
#   newnames
# )


# 3: Processing new images

# 3.1: S1  (2023-02-22)

# Load files

# dir_s1 <- dir_dat %>%
#   paste0(., "/new_s1/")
# 
# files_s1 <- dir_s1 %>%
#   list.files(full.names = TRUE) 
# 
# names_s1 <- files_s1 %>%
#   basename()
# 
# for (i in 1:length(files_s1)) {
#   dtyp <- files_s1[i] %>%
#     raster() %>%
#     dataType()
# 
#   cov_i <- files_s1[i] %>%
#     terra::rast()
# 
#   crs(cov_i) <- mycrs
#   
#   # Crop to coast
#   r1 <- terra::mask(
#     cov_i,
#     mask = dem
#   )
#   
#   # Replace 0 and above
#   ifel(
#     r1 < 0,
#     r1,
#     NA,
#     filename = paste0(
#       dir_cov,
#       names_s1[i]
#     ),
#     overwrite = TRUE,
#     datatype = dtyp
#   )
# }

# 3.2: S2 (2023-02-24)

# Load files

# dir_s2 <- dir_dat %>%
#   paste0(., "/new_s2/")
# 
# files_s2 <- dir_s2 %>%
#   list.files(full.names = TRUE)
# 
# names_s2 <- files_s2 %>%
#   basename()
# 
# tmpfolder <- paste0(dir_dat, "/Temp/")
# 
# terraOptions(tempdir = tmpfolder)
# 
# for (i in 1:length(files_s2)) {
#   dtyp <- files_s2[i] %>%
#     raster() %>%
#     dataType()
# 
#   cov_i <- files_s2[i] %>%
#     terra::rast()
# 
#   crs(cov_i) <- mycrs
# 
#   # Crop to coast
#   r1 <- terra::mask(
#     cov_i,
#     mask = dem,
#     filename = paste0(
#       tmpfolder,
#       "temp1.tif"
#     ),
#     overwrite = TRUE
#   )
# 
#   # Keep only 0 - 10000
#   ifel(
#     r1 > 0 & r1 <= 10000,
#     r1,
#     NA,
#     filename = paste0(
#       dir_cov,
#       names_s2[i]
#     ),
#     overwrite = TRUE,
#     datatype = dtyp
#   )
# }

# 3.3: Cost distance layer (2023-03-23)

# costdist <- dir_dat %>%
#   paste0(., "/Cost_Distance/Cost_Dist_DHM2015_terraen_10m.tif") %>%
#   rast()
# 
# costdist_newname <- "cost_dist"
# 
# names(costdist) <- costdist_newname
# 
# outname <- dir_cov %>%
#   paste0(., "/", costdist_newname, ".tif")
# 
# r1 <- terra::math(
#   costdist,
#   "round",
#   digits = 2,
#   filename = paste0(
#     tmpfolder,
#     "temp1.tif"
#   ),
#   overwrite = TRUE
# )
#   
# # Crop to extent
# terra::crop(
#   r1,
#   y = dem,
#   filename = outname,
#   overwrite = TRUE,
#   datatype = "FLT4S"
# )

# 3.4: Detrended DEM (2023-03-23)

# detrended <- dir_dat %>%
#   paste0(., "/Detrended_DEM/") %>%
#   list.files(full.names = TRUE) %>%
#   rast()
# 
# detrended_newnames <- sources(detrended) %>%
#     basename() %>%
#     file_path_sans_ext() %>%
#     gsub("\\.", "_", .) %>%
#     gsub("-", "_", .) %>%
#     tolower()
# 
# names(detrended) <- detrended_newnames
# 
# r1 <- terra::math(
#   detrended,
#   "round",
#   digits = 2,
#   filename = paste0(
#     tmpfolder,
#     "temp1.tif"
#   ),
#   overwrite = TRUE
# )
# 
# # Crop to extent
# r2 <- terra::crop(
#   r1,
#   y = dem,
#   filename = paste0(
#     tmpfolder,
#     "temp2.tif"
#   ),
#   overwrite = TRUE
# )
# 
# for (i in 1:nlyr(r2)) {
#   outname <- dir_cov %>%
#     paste0(., "/", detrended_newnames[i], ".tif")
#   
#   writeRaster(
#     r2[[i]],
#     datatype = "FLT4S",
#     filename = outname,
#     overwrite = TRUE
#   )
# }

# 3.5: Hillyness (2023-03-23)

# hillyness <- dir_dat %>%
#   paste0(., "/hillyness/hillyness.tif") %>%
#   rast()
# 
# names(hillyness) <- "hillyness"
# NAflag(hillyness) <- 128
# 
# outname_hillyness <- dir_cov %>%
#   paste0(., "/hillyness.tif")
# 
# terra::crop(
#   hillyness,
#   dem,
#   filename = outname_hillyness,
#   datatype = "INT1U",
#   overwrite = TRUE
# )

# Edit: Mask hillyness with DEM

# r <- dir_dat %>%
#   paste0(., "/hillyness.tif") %>%
#   rast()
# 
# terra::mask(
#   r,
#   mask = dem,
#   filename = outname_hillyness,
#   datatype = "INT1U",
#   overwrite = TRUE
# )

# 4: Rename the layers for all covariates (2023-03-23)

# cov_files <- dir_cov %>%
#   list.files(
#     pattern = ".tif",
#     full.names = TRUE
#   )
# 
# dir_cov_renamed <- dir_dat %>%
#   paste0(., "/covariates_renamed/") %T>%
#   dir.create()
# 
# library(parallel)
# 
# numCores <- detectCores()
# numCores
# 
# showConnections()
# 
# cl <- makeCluster(numCores)
# 
# clusterEvalQ(
#   cl,
#   {
#     library(terra)
#     library(raster)
#     library(magrittr)
#     library(dplyr)
#     library(tools)
#   }
# )
# 
# clusterExport(
#   cl,
#   c("mycrs",
#     "dir_cov_renamed",
#     "cov_files",
#     "tmpfolder"
#   )
# )
# 
# parSapplyLB(
#   cl,
#   cov_files,
#   function(x) {
#     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
#     
#     r <- x %>% rast()
#     
#     dtyp <- x %>%
#       raster() %>%
#       dataType()
#     
#     newname_x <- sources(r) %>%
#       basename() %>%
#       file_path_sans_ext() %>%
#       gsub("\\.", "_", .) %>%
#       gsub("-", "_", .) %>%
#       tolower()
#     
#     crs(r) <- mycrs
#     names(r) <- newname_x
#     
#     outname_x <- dir_cov_renamed %>%
#       paste0(., "/", newname_x, ".tif")
#     
#     writeRaster(
#       r,
#       datatype = dtyp,
#       filename = outname_x,
#       overwrite = TRUE
#     )
#     
#     out <- rast(outname_x)
#     
#     return(out)
#   }
# )
# 
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)
# 
# dir_cov_renamed %>%
#   list.files(full.names = TRUE) %>%
#   rast() %>% names()


# 5: Update names in covariate table (2023-02-27)

# cov_cats <- dir_code %>%
#   paste0(., "/cov_categories_20230202.csv") %>%
#   read.table(
#     sep = ";",
#     header = TRUE,
#     encoding = "latin1"
#   )
# 
# newnames <- cov_cats$name %>%
#   gsub("\\.", "_", .) %>%
#   gsub("-", "_", .) %>%
#   tolower()
# 
# cov_cats$name <- newnames
# 
# write.table(
#   cov_cats,
#   file = paste0(dir_code, "/cov_categories_20230227.csv"),
#   row.names = FALSE,
#   sep = ";",
#   fileEncoding = "latin1"
# )

# 6: Identify covariates missing in the overview table

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230323.csv") %>%
  read.table(
    sep = ";",
    header = TRUE,
    encoding = "latin1"
  )

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

cov_names <- cov_files %>%
  basename() %>%
  file_path_sans_ext()

setdiff(cov_names, cov_cats$name)

# 7: Crop all covariates for the test area

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

squareshape <- dir_dat %>%
  paste0(., "/testarea_10km/square10km.shp") %>%
  vect

square_ext <- squareshape %>%
  ext %>%
  round(-1)

outfolder <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

outfolder %>% dir.create()

source("f_cropstack.R")

cov_files %>%
  cropstack(
    y = square_ext,
    folder = outfolder
  )

# END