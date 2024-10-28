# Calculating uncertainties for peat2022
# Using cube root transformed data to calculate distribution curves
# Back-transforming results to linear scale
# It seems the final map is a combination from two sets of models.

library(terra)
library(caret)
library(Cubist)
library(magrittr)
library(ModelMetrics)
library(dplyr)

# Getting started

soc6pct_cbrt <- 6^(1/3)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

dir_peat <- paste0(root, "/peat_stuff_2022/")
dir_peat_incl_S2 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/",
  "cubist/incl_S2/incl_nirs/bst/feb_2024_additional_areas/"
  )
dir_peat_excl_S2 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/",
  "cubist/excl_S2/soc/incl_nirs/bst/feb_2024_additional_areas/"
)
dir_new_maps <- dir_peat %>% paste0(., "/new_maps/") %T>% dir.create()

testn <- 14
mycrs <- "EPSG:25832"

# Test if the map is a combination from two sets of models.
# "incl_S2" and "excl_S2"

combined_soc_map <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/cubist/",
  "/combinations/feb_2024/",
  "cubist_pred_soc_all_final_feb2024.tif") %>%
  rast()

# soc_map_incl_S2 <- dir_peat_incl_S2 %>%
#   paste0(., "/maps/res_maps/meanPred_back_tr.tif") %>%
#   rast()
# soc_map_excl_S2 <- dir_peat_excl_S2 %>%
#   paste0(., "/maps/res_maps/meanPred_back_tr.tif") %>%
#   rast()
# 
# stack1 <- c(combined_soc_map, soc_map_incl_S2, soc_map_excl_S2)
# names(stack1) <- c("combined", "incl_s2", "excl_s2")
# 
# sample1 <- spatSample(stack1, 10000)
# 
# plot(sample1)

# It seems like the final maps aren't based on the bootstrap repetitions?
# Try these maps instead:

dir_peat_incl_S2_overall2 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/cubist/",
  "/incl_S2/incl_nirs/feb_2024_additional_areas/"
)
dir_peat_excl_S2_overall2 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/cubist/",
  "/excl_S2/soc/incl_nirs/feb_2024_additional_areas/"
)

# soc_map_incl_S2_2 <- dir_peat_incl_S2_overall2 %>%
#   paste0(., "/cubist_prediction_cbrtsoc_nirs.tif") %>%
#   rast()
# soc_map_excl_S2_2 <- dir_peat_excl_S2_overall2 %>%
#   paste0(., "/cubist_prediction_cbrtsoc_nirs.tif") %>%
#   rast()
# 
# stack2 <- c(combined_soc_map, soc_map_incl_S2_2, soc_map_excl_S2_2)
# names(stack2) <- c("combined", "incl_s2", "excl_s2")
# 
# sample2 <- spatSample(stack2, 10000)
# 
# plot(sample2)
# 
# sample2 %>%
#   filter(is.na(incl_s2)) %>%
#   select(-incl_s2) %>%
#   plot()
# 
# sample2 %>%
#   filter(!is.na(incl_s2)) %>%
#   plot()

# The results above are better, but the alignment is still not quite perfect.


# Try the maps delivered on October 2023.
# (for some reason this folder is empty: [...]
# peat_stuff_2022\peat2022_raw\data_from_2009_2010_2020_2021\predictions\cubist\
# incl_S2\incl_nirs\oct_2023)
# Not known why the folder is empty.
# Apparently, this folder was last edited on 2023-11-02.

dir_peat_incl_S2_overall3 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/cubist/",
  "/incl_S2/incl_nirs/all_combined/"
)
dir_peat_excl_S2_overall3 <- paste0(
  dir_peat, "/peat2022_raw/data_from_2009_2010_2020_2021/predictions/cubist/",
  "/excl_S2/soc/incl_nirs/oct_2023/"
)

# soc_map_incl_S2_3 <- dir_peat_incl_S2_overall3 %>%
#   paste0(., "/cubist_prediction_cbrtsoc_nirs.tif") %>%
#   rast() %>%
#   extend(y = ext(combined_soc_map)) %>%
#   crop(y = ext(combined_soc_map))
# soc_map_excl_S2_3 <- dir_peat_excl_S2_overall3 %>%
#   paste0(., "/cubist_prediction_cbrtsoc_nirs.tif") %>%
#   rast() %>%
#   extend(y = ext(combined_soc_map)) %>%
#   crop(y = ext(combined_soc_map))
# 
# stack3 <- c(combined_soc_map, soc_map_incl_S2_3, soc_map_excl_S2_3)
# names(stack3) <- c("combined", "incl_s2", "excl_s2")
# 
# sample3 <- spatSample(stack3, 10000)
# 
# plot(sample3)
# 
# sample3 %>%
#   filter(is.na(incl_s2)) %>%
#   select(-incl_s2) %>%
#   plot()  # Perfect correlation here
# 
# sample3 %>%
#   filter(!is.na(incl_s2)) %>%
#   plot()
# 
# sample3 %>%
#   filter(!is.na(incl_s2)) %>%
#   mutate(
#     combined = combined^(1/3)
#   ) %>%
#   plot() # OK correlation here, probably not quite the the right data set.

# It seems like it is correct to assume that the maps are a combination of two
# sets of maps:
# - excl_S2: Covering the entire wetland areas.
# - incl_S2: Covering only areas with bare soil
# The final map uses "incl_S2" over "excl_S2".

# Two do:
# Get mean cube root SOC:
# - I will use the mean SOC values from the final map also used in Kulstof2022.
# -- This map will have to be converted to cube root scale.
# Get cube root MSE for "incl_S2" and "excl_S2".
# - The mean values were predicted without bootstrapping, so I will use the MSE
#   from the cross validation of these two overall models.
# - The cross validation results from October 2023 are missing for "incl_S2".
#   I will therefore calculate the MSE based on the results from February 2024
#   in both cases for sake of consistency.
# - The models were trained on cube root transformed data, so the MSE values
#   will be usable without need of transformation.
# Get prediction variance for "incl_S2" and "excl_S2".
# - Use bootstrap repetitions for this (February 2024 for consistency).
# Calculate prediction standard error (PSE) for "incl_S2" and "excl_S2". Then
# merge.
# Calculate peat probability using cube root mean and PSE layers.
# Calculate prediction quantiles on cube root scale, then back-transform them to
# linear scale.
# - Remove quantiles below 0 and above 60.
# - Round off decimals.
# Use 2022 data to validate uncertainty layers.

# Get mean cube root SOC:

# soc_mean_cuberoot <- combined_soc_map^(1/3)

soc_mean_cuberoot_filename <- dir_new_maps %>%
  paste0(., "soc_mean_cuberoot.tif")

# writeRaster(
#   soc_mean_cuberoot,
#   filename = soc_mean_cuberoot_filename,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )

soc_mean_cuberoot <- rast(soc_mean_cuberoot_filename)

# Get MSE for the two overall models

mse_excl_S2 <- dir_peat_excl_S2_overall2 %>%
  paste0(., "/test_metrics/res_cub_cbrt_soc_train.csv") %>%
  read.table(., header = TRUE, sep = ",") %>%
  summarise(
    mse = mse(fit_cubist.pred.obs, fit_cubist.pred.pred)
  ) %>%
  unlist() %>%
  unname()
# 0.376681

mse_incl_S2 <- dir_peat_incl_S2_overall2 %>%
  paste0(., "/test_metrics/res_cub_cbrt_soc_train - Copy (2).csv") %>%
  read.table(., header = TRUE, sep = ";") %>%
  summarise(
    mse = mse(fit_cubist.pred.obs, fit_cubist.pred.pred)
  ) %>%
  unlist() %>%
  unname()
# 0.1916157

# Calculate prediction variance for both sets of bootstrap models

bootmaps_excl_S2 <- dir_peat_excl_S2 %>%
  paste0(., "/maps/") %>%
  list.files(., pattern = ".tif", full.names = TRUE) %>%
  rast()

bootmaps_incl_S2 <- dir_peat_incl_S2 %>%
  paste0(., "/maps/") %>%
  list.files(., pattern = ".tif", full.names = TRUE) %>%
  rast()

varmap_excl_S2_outfile <- dir_new_maps %>%
  paste0(., "soc_variance_crbt_excl_S2.tif")
varmap_incl_S2_outfile <- dir_new_maps %>%
  paste0(., "soc_variance_crbt_incl_S2.tif")

# varmap_excl_S2 <- stdev(bootmaps_excl_S2, pop = FALSE, na.rm = TRUE)^2
# names(varmap_excl_S2) <- "var"
# crs(varmap_excl_S2) <- mycrs
# varmap_incl_S2 <- stdev(bootmaps_incl_S2, pop = FALSE, na.rm = TRUE)^2
# names(varmap_incl_S2) <- "var"
# crs(varmap_incl_S2) <- mycrs

# writeRaster(
#   varmap_excl_S2,
#   filename = varmap_excl_S2_outfile,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )
# 
# writeRaster(
#   varmap_incl_S2,
#   filename = varmap_incl_S2_outfile,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )

varmap_excl_S2 <- varmap_excl_S2_outfile %>% rast()
varmap_incl_S2 <- varmap_incl_S2_outfile %>% rast()


# Calculate prediction standard error

# pse_excl_S2 <- sqrt(varmap_excl_S2 + mse_excl_S2)
# pse_incl_S2 <- sqrt(varmap_incl_S2 + mse_incl_S2)

pse_excl_S2_outfile <- dir_new_maps %>%
  paste0(., "soc_pse_crbt_excl_S2.tif")
pse_incl_S2_outfile <- dir_new_maps %>%
  paste0(., "soc_pse_crbt_incl_S2.tif")

# names(pse_excl_S2) <- "pse"
# names(pse_incl_S2) <- "pse"
# 
# writeRaster(
#   pse_excl_S2,
#   filename = pse_excl_S2_outfile,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )
# 
# writeRaster(
#   pse_incl_S2,
#   filename = pse_incl_S2_outfile,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )

pse_excl_S2 <- pse_excl_S2_outfile %>% rast()
pse_incl_S2 <- pse_incl_S2_outfile %>% rast()

# Combine layers of the prediction standard error

pse_combined_outfile <- dir_new_maps %>%
  paste0(., "/soc_pse_crbt_combined.tif")

# pse_combined <- ifel(
#   test = is.na(pse_incl_S2),
#   yes = pse_excl_S2,
#   no = pse_incl_S2,
#   filename = pse_combined_outfile,
#   names = "pse",
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S"
# )

pse_combined <- pse_combined_outfile %>% rast()


# Split mean and pse into tiles

mean_pse_cbrt <- c(soc_mean_cuberoot, pse_combined) 
names(mean_pse_cbrt) <- c("mean", "pse")

dir_tiles_in <- dir_new_maps %>%
  paste0(., "/tiles_in/") %T>%
  dir.create()

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

tile_shapes <- dir_tiles %>%
  paste0(., "/tiles.shp") %>%
  vect()

tile_numbers_chr <- tile_shapes %>%
  values() %>%
  unlist() %>%
  unname() %>%
  formatC(width = 3, flag = "0")

dir_tiles_mean_pse <- dir_tiles_in %>%
  paste0(., "/mean_pse_cbrt/") %T>%
  dir.create()

# for (i in 1:length(tile_shapes)) {
#   terra::crop(
#     mean_pse_cbrt,
#     tile_shapes[i],
#     filename = paste0(
#       dir_tiles_mean_pse, "/mean_pse_cbrt_tile_", tile_numbers_chr[i], ".tif"
#     )
#   )
# }

tiles_mean_pse_cbrt <- dir_tiles_mean_pse %>%
  list.files(full.names = TRUE)

dir_tiles_out <- dir_new_maps %>%
  paste0(., "/tiles_out/") %T>%
  dir.create()

# Calculate peat probability

library(parallel)
library(foreach)

dir_tiles_prob_peat <- dir_tiles_out %>%
  paste0(., "/peat_probability_t22/") %T>%
  dir.create()

calc_prob <- function (x) {
  out <- pnorm(
    q = soc6pct_cbrt,
    mean = x[1],
    sd = x[2],
    lower.tail = FALSE
  ) %>%
    round(., digits = 3) %>%
    multiply_by(100)
  return(out)
}

numCores <- detectCores() - 1
numCores

# showConnections()
# 
# cl <- makeCluster(numCores)
# 
# clusterEvalQ(
#   cl,
#   {
#     library(terra)
#     library(magrittr)
#     library(dplyr)
#     library(tools)
#   }
# )
# 
# clusterExport(
#   cl,
#   c("dir_tiles_prob_peat",
#     "tile_numbers_chr",
#     "soc6pct_cbrt",
#     "calc_prob",
#     "tiles_mean_pse_cbrt",
#     "dir_dat"
#   )
# )
# 
# parSapplyLB(
#   cl,
#   1:length(tiles_mean_pse_cbrt),
#   function(x) {
#     tmpfolder <- paste0(dir_dat, "/Temp/")
#     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
#     
#     mean_pse_cbrt_x <- tiles_mean_pse_cbrt[x] %>%
#       rast()
#   
#     prob_peat_outfile_x <- dir_tiles_prob_peat %>%
#       paste0(., "/peat_prob_tile", tile_numbers_chr[x], ".tif")
#     
#     app(
#       mean_pse_cbrt_x,
#       fun = calc_prob,
#       filename = prob_peat_outfile_x,
#       overwrite = TRUE,
#       wopt = list(
#         gdal = "TILED=YES",
#         datatype = "FLT4S"
#       )  
#     )
# 
#     return(NULL)
#   }
# )
# 
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)
# 
# tiles_prob_peat <- dir_tiles_prob_peat %>%
#   list.files(full.names = TRUE)
# 
# tile1_i <- tiles_prob_peat[1] %>% rast()
# tile1_i
# 
# dtyp_i <- datatype(tile1_i)
# naflag_i <- NAflag(tile1_i)
# 
# if (!is.finite(naflag_i)) {
#   naflag_i <- -1
# }
# 
# outtiles_sprc <- tiles_prob_peat %>% sprc()

outfile_basename <- "peat_probability"

outfile_peat_probability <- dir_new_maps %>%
  paste0(., "/", outfile_basename, ".tif")

# merge(
#   outtiles_sprc,
#   filename = outfile_peat_probability,
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = dtyp_i,
#   NAflag = naflag_i,
#   names = outfile_basename
# )

peat2022_probability <- dir_new_maps %>%
  paste0(., "/peat_probability.tif") %>%
  rast()


# Calculate quantiles

rev((100 - c(68, 90, 95, 99)) / 2)
(c(68, 90, 95, 99) + 100) / 2

prob_q_out <- c(0.5, 2.5, 5.0, 16.0, 84.0, 95.0, 97.5, 99.5)/100
prob_q_out_chr <- prob_q_out %>%
  multiply_by(1000) %>%
  formatC(width = 4, flag = "0") %>%
  paste0("p", .)
  
dir_tiles_peat_quantiles <- dir_tiles_out %>%
  paste0(., "/peat_quantiles_t22/") %T>%
  dir.create()

calc_q_cbrt <- function (x) {
  out <- qnorm(
    p = prob_q_out,
    mean = x[1],
    sd = x[2]
  ) %>%
    raise_to_power(3) %>%
    round(., digits = 1) %>%
    set_names(prob_q_out_chr)
  out[out < 0] <- 0
  out[out > 60] <- 60
  return(out)
}

numCores <- detectCores() - 1
numCores

# showConnections()
# 
# cl <- makeCluster(numCores)
# 
# clusterEvalQ(
#   cl,
#   {
#     library(terra)
#     library(magrittr)
#     library(dplyr)
#     library(tools)
#   }
# )
# 
# clusterExport(
#   cl,
#   c("dir_tiles_prob_peat",
#     "tile_numbers_chr",
#     "prob_q_out",
#     "calc_q_cbrt",
#     "tiles_mean_pse_cbrt",
#     "dir_dat",
#     "prob_q_out_chr",
#     "dir_tiles_peat_quantiles"
#   )
# )
# 
# parSapplyLB(
#   cl,
#   1:length(tiles_mean_pse_cbrt),
#   function(x) {
#     tmpfolder <- paste0(dir_dat, "/Temp/")
#     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
#     
#     mean_pse_cbrt_x <- tiles_mean_pse_cbrt[x] %>%
#       rast()
#     
#     outfile_tmp_x <- tmpfolder %>%
#       paste0(., "/peat_qs_tile", tile_numbers_chr[x], ".tif")
#     
#     app(
#       mean_pse_cbrt_x,
#       fun = calc_q_cbrt,
#       filename = outfile_tmp_x,
#       overwrite = TRUE,
#       wopt = list(
#         gdal = "TILED=YES",
#         datatype = "FLT4S"
#       )  
#     )
#     
#     my_qs <- outfile_tmp_x %>% rast()
#     
#     outdir_tile_x <- dir_tiles_peat_quantiles %>%
#       paste0(., "/tile_", tile_numbers_chr[x], "/") %T>%
#       dir.create()
#     
#     for (i in 1:nlyr(my_qs)) {
#       outname_base <- paste0("soc_", prob_q_out_chr[i])
#       
#       outfile_xi <- outdir_tile_x %>%
#         paste0(., outname_base, ".tif")
#       
#       writeRaster(
#         my_qs[[i]],
#         filename = outfile_xi,
#         names = outname_base,
#         overwrite = TRUE,
#         gdal = "TILED=YES",
#         datatype = "FLT4S"
#       )
#     }
#     
#     return(NULL)
#   }
# )
# 
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)

# Merge tiles

tilenames <- tile_numbers_chr %>%
  paste0("tile_", .)

outfiles_basenames <- paste0(
  dir_tiles_peat_quantiles, "/", tilenames[1], "/") %>%
  list.files() %>%
  basename() %>%
  tools::file_path_sans_ext()

# for (i in 1:length(outfiles_basenames)) {
#   summary_tiles_i <- paste0(
#     dir_tiles_peat_quantiles, "/", tilenames, "/",
#     outfiles_basenames[i], ".tif"
#   )
#   
#   tile1_i <- summary_tiles_i[1] %>% rast()
#   
#   dtyp_i <- datatype(tile1_i)
#   naflag_i <- NAflag(tile1_i)
#   
#   if (dtyp_i == "INT1U") { naflag_i <- 101 }
#   if (!is.finite(naflag_i)) { naflag_i <- -1}
#   
#   outtiles_sprc <- summary_tiles_i %>% sprc()
#   
#   merge(
#     outtiles_sprc,
#     filename = paste0(dir_new_maps, "/", outfiles_basenames[i], ".tif"),
#     overwrite = TRUE,
#     gdal = "TILED=YES",
#     datatype = dtyp_i,
#     NAflag = naflag_i,
#     names = outfiles_basenames[i]
#   )
# }

# Load combination raster showing parts covered by Tørv2022 and my own map.
# Combination == 1 shows peat2022 extent
dir_final_maps_soc <- dir_peat %>%
  paste0(., "/final_combined_soc_maps/") %T>%
  dir.create()

SOC_combination_map <- root %>%
  paste0(
    ., "/Soil_maps_10m_new/Kulstof2022/SOC_000_030_cm/",
    "SOC_combination_000_030_cm.tif") %>%
  rast() %>%
  crop(
    y = ext(peat2022_probability),
    extend = TRUE,
    filename = paste0(dir_final_maps_soc, "SOC_combination_000_030_cm.tif"),
    overwrite = TRUE,
    gdal = "TILED=YES",
    datatype = "INT1U",
    NAflag = 101,
    names = "SOC_combination_000_030_cm"
    )

# Combine with Digijord probabilities and quantiles

basenames_finalmaps <- c(
  "peat_probability", "soc_p0005", "soc_p0025", "soc_p0050", "soc_p0160",
  "soc_p0840", "soc_p0950", "soc_p0975", "soc_p0995"
  )

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")
dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")
dir_pred_boot <- dir_boot %>%
  paste0(., "/predictions/")
dir_merged_depth <- dir_pred_boot %>%
  paste0(., "/final_maps/depth_000_030_cm/")
    
outmaps_peat2022 <- list()
outmaps_digijord <- list()
outmaps_complete <- list()

tmpfolder <- paste0(dir_dat, "/Temp/")
terraOptions(memfrac = 0.02, tempdir = tmpfolder)

# for (i in 1:length(basenames_finalmaps)) {
#   outmaps_peat2022[[i]] <- paste0(
#     dir_new_maps, "/", basenames_finalmaps[i], ".tif"
#   ) %>%
#     rast()
#   
#   outmaps_digijord[[i]] <- paste0(
#     dir_merged_depth, "/", basenames_finalmaps[i], ".tif"
#     ) %>%
#     rast()
#   
#   outmaps_complete[[i]] <- ifel(
#     test = SOC_combination_map == 1,
#     yes = outmaps_peat2022[[i]],
#     no = outmaps_digijord[[i]]
#   ) %>%
#     mask(
#       mask = SOC_combination_map,
#       filename = paste0(dir_final_maps_soc, basenames_finalmaps[i], ".tif"),
#       overwrite = TRUE,
#       gdal = "TILED=YES",
#       datatype = datatype(outmaps_digijord[[i]]),
#       NAflag = -1,
#       names = basenames_finalmaps[i]
#     )
# }

# Re-mask previous Kulstof2022 maps (remove lakes)
 
dir_old_SOC_maps <- root %>%
  paste0(., "/Soil_maps_10m_new/Kulstof2022/SOC_000_030_cm/")

list.files(
  dir_old_SOC_maps,
  pattern = ".tif"
  )

maps_to_crop <- c(
  "SOC_class_000_030_cm_combined",  
  "SOC_mean_000_030_cm_combined"
)

for (i in 1:length(maps_to_crop)) {
  map_old <- dir_old_SOC_maps %>%
    paste0(., maps_to_crop[i], ".tif") %>%
    rast()
  
  dtyp_i <- datatype(map_old)
  naflag_i <- NAflag(map_old)
  
  if (dtyp_i == "INT1U") { naflag_i <- 101 }
  if (!is.finite(naflag_i)) { naflag_i <- -1}
  
  map_old %>%
    mask(
      mask = SOC_combination_map,
      filename = paste0(dir_final_maps_soc, maps_to_crop[i], ".tif"),
      overwrite = TRUE,
      gdal = "TILED=YES",
      datatype = dtyp_i,
      NAflag = naflag_i,
      names = maps_to_crop[i]
    )
}

# END