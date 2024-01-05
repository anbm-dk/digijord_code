# 01: Crop and rename covariates

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
#     rast() %>%
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
#     rast() %>%
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
#     rast() %>%
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

# 3.6: Fill holes in terodep (2023-12-19)
# 
# terodep_files <- cov_files %>%
#   grep('terodep', ., value = TRUE)
# 
# dir_cov_renamed <- dir_dat %>%
#   paste0(., "/covariates_renamed/") %T>%
#   dir.create()
# 
# for (i in 1:length(terodep_files)) {
#   r <- terodep_files[i] %>% rast()
#   dtyp <- datatype(r)
#   newname_x <- sources(r) %>%
#     basename() %>%
#     file_path_sans_ext() %>%
#     gsub("\\.", "_", .) %>%
#     gsub("-", "_", .) %>%
#     tolower()
#   crs(r) <- mycrs
#   names(r) <- newname_x
#   outname_x <- dir_cov_renamed %>%
#     paste0(., "/", newname_x, ".tif")
#   ifel(
#     test = is.na(r),
#     yes = dem*0,
#     no = r,
#     datatype = dtyp,
#     filename = outname_x,
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# 3.6: Fill holes in cos_aspect_radians and sin_aspect_radians (2023-12-19)

# aspect_files <- cov_files %>%
#   grep('aspect_radians', ., value = TRUE)
# 
# dir_cov_renamed <- dir_dat %>%
#   paste0(., "/covariates_renamed/") %T>%
#   dir.create()
# 
# for (i in 1:length(aspect_files)) {
#   r <- aspect_files[i] %>% rast()
#   dtyp <- datatype(r)
#   newname_x <- sources(r) %>%
#     basename() %>%
#     file_path_sans_ext() %>%
#     gsub("\\.", "_", .) %>%
#     gsub("-", "_", .) %>%
#     tolower()
#   crs(r) <- mycrs
#   names(r) <- newname_x
#   outname_x <- dir_cov_renamed %>%
#     paste0(., "/", newname_x, ".tif")
#   ifel(
#     test = is.na(r),
#     yes = dem*0,
#     no = r,
#     datatype = dtyp,
#     filename = outname_x,
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# 4: Update names in covariate table (2023-02-27)

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

# 5: Make new sets of oblique geographic coordinates

# install.packages('devtools')

# library(devtools)

# install_github("anbm-dk/obliquer")

# library(obliquer)
#
# dir_tiles <- dir_dat %>%
#   paste0(., "/tiles_591/")
#
# tile_shapes <- dir_tiles %>%
#   paste0(., "/tiles.shp") %>%
#   vect()
#
# # split dem into tiles
#
# tmp_dem_tiles <- paste0(tmpfolder, "/dem/") %T>% dir.create()
# tmp_ogc_tiles <- paste0(tmpfolder, "/ogc/") %T>% dir.create()
#
# for (i in 1:length(tile_shapes)) {
#   terra::crop(
#     dem,
#     tile_shapes[i],
#     filename = paste0(tmp_dem_tiles, "/dem_tile_", i, ".tif")
#   )
# }
#
# dem_files <- tmp_dem_tiles %>% list.files(full.names = TRUE)
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
#     library(magrittr)
#     library(obliquer)
#   }
# )
#
# clusterExport(
#   cl,
#   c(
#     "dem_files",
#     "tmp_ogc_tiles",
#     "tmpfolder"
#   )
# )
#
# parSapplyLB(
#   cl,
#   1:length(dem_files),
#   function(j) {
#     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
#
#     dem_j <- dem_files[j] %>% rast()
#
#     obliquify(
#       dem_j,
#       n_angles = 64,
#       n_digits = 0,
#       digits_names = 3,
#       filename = paste0(tmp_ogc_tiles, "/ogcs_tile_", j, ".tif"),
#       datatype = "INT4S"
#       )
#
#     return(NULL)
#   }
# )
#
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)
#
# ogc_files <- tmp_ogc_tiles %>% list.files(full.names = TRUE)
#
# ogc_names <- ogc_files[1] %>% rast() %>% names()
#
# for (i in 1:length(ogc_names)) {
#   ogcs_i <- ogc_files %>% lapply(
#     function(x) {
#       out <- x %>% rast() %>% subset(i)
#       return(out)
#     }
#   )
#
#   ogcs <- sprc(ogcs_i)
#
#   ogcs_merged <- merge(
#     ogcs,
#     filename = paste0(dir_cov, "/ogc_", ogc_names[i], ".tif"),
#     datatype = "INT4S",
#     gdal = "TILED=YES"
#   )
# }

# 6: Rename the layers for all covariates (2023-03-23)

# First check in the names match

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

cov_origin <- cov_files %>% rast()
cov_origin_lyrnames <- names(cov_origin)

cov_origin_filenames <- cov_files %>% basename() %>% tools::file_path_sans_ext()

mismatches <- cbind(
  cov_origin_filenames,
  cov_origin_lyrnames
  )[cov_origin_filenames != cov_origin_lyrnames]

mismatches

changethese <- cov_files[cov_origin_filenames != cov_origin_lyrnames]

# Correct it if they do not match

if (length(changethese) > 0) {
  
  dir_cov_renamed <- dir_dat %>%
    paste0(., "/covariates_renamed/") %T>%
    dir.create()
  
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
    c("mycrs",
      "dir_cov_renamed",
      "changethese",
      "tmpfolder"
    )
  )
  
  parSapplyLB(
    cl,
    changethese,
    function(x) {
      terraOptions(memfrac = 0.02, tempdir = tmpfolder)
      
      r <- x %>% rast()
      
      dtyp <- datatype(r)
      
      newname_x <- sources(r) %>%
        basename() %>%
        file_path_sans_ext() %>%
        gsub("\\.", "_", .) %>%
        gsub("-", "_", .) %>%
        tolower()
      
      crs(r) <- mycrs
      names(r) <- newname_x
      
      outname_x <- dir_cov_renamed %>%
        paste0(., "/", newname_x, ".tif")
      
      writeRaster(
        r,
        datatype = dtyp,
        filename = outname_x,
        overwrite = TRUE,
        gdal = "TILED=YES"
      )
      
      return(NA)
    }
  )
  
  stopCluster(cl)
  foreach::registerDoSEQ()
  rm(cl)
  
  dir_cov_renamed %>%
    list.files(full.names = TRUE) %>%
    rast() %>%
    names()
}

# 7: Identify covariates missing in the overview table

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

cov_names <- cov_files %>%
  basename() %>%
  file_path_sans_ext()

data.frame(new = setdiff(cov_names, cov_cats$name))

# 8: Crop all covariates for the test area

cov_files <- dir_cov %>%
  list.files(
    pattern = ".tif",
    full.names = TRUE
  )

squareshape <- dir_dat %>%
  paste0(., "/testarea_10km/square10km.shp") %>%
  vect()

square_ext <- squareshape %>%
  ext() %>%
  round(-1)

mask_LU <- paste0(dir_dat, "/layers/Mask_LU.tif") %>% rast()

mask_LU_10km <- crop(
  mask_LU,
  square_ext
)

outfolder <- dir_dat %>%
  paste0(., "/testarea_10km/covariates/")

outfolder %>% dir.create()

source("f_cropstack.R")

cov_files %>%
  cropstack(
    y = mask_LU_10km,
    folder = outfolder,
    mask = TRUE
  )

# END