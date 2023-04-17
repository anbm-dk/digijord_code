# 08: Script for making maps

library(parallel)
library(caret)
library(terra)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 8
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")


fractions <- c("clay", "silt", "fine_sand", "coarse_sand", "logSOC", "logCaCO3")

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)


bounds_lower <- c(0, 0, 0, 0, NA, NA)
bounds_upper <- c(100, 100, 100, 100, log(100), log(100))

dir_cov <- dir_dat %>% paste0(., "/covariates")
cov_files <- dir_cov %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230323.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

cov_selected <- cov_cats %>%
  filter(anbm_use == 1) %>%
  select(name) %>%
  unlist() %>%
  unname()

source("f_predict_passna.R")

# Load models

models_loaded <- lapply(
  1:6,
  function(x) {
    out <- fractions[x] %>%
      paste0(dir_results, "/model_", ., ".rds") %>%
      readRDS()
    return(out)
  }
)

# Tiles for model prediction

numCores <- detectCores()
numCores

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

subdir_tiles <- dir_tiles %>% list.dirs() %>% .[-1]

dir_pred_all <- dir_results %>%
  paste0(., "/predictions/") %T>%
  dir.create()

dir_pred_tiles <- dir_pred_all  %>%
  paste0(., "/tiles/") %T>%
  dir.create()

for (i in 1:length(fractions)) {
  frac <- fractions[i]
  
  dir_pred_tiles_frac <- dir_pred_tiles %>%
    paste0(., "/", names(models)[i], "/") %T>%
    dir.create()
  
  model_i <- models[[i]]
  
  showConnections()
  
  cl <- makeCluster(numCores)
  
  clusterEvalQ(
    cl,
    {
      library(terra)
      library(caret)
      library(Cubist)
      library(magrittr)
      library(dplyr)
    }
  )
  
  clusterExport(
    cl,
    c("model_i",
      "subdir_tiles",
      "dir_pred_tiles_frac",
      "frac",
      "cov_names",
      "cov_selected",
      "predict_passna",
      "dir_dat"
    )
  )
  
  parSapplyLB(
    cl,
    1:length(subdir_tiles),
    function(x) {
      tmpfolder <- paste0(dir_dat, "/Temp/")
      
      terraOptions(memfrac = 0.02, tempdir = tmpfolder)
      
      cov_x <- subdir_tiles[x] %>%
        list.files(full.names = TRUE) %>%
        rast()
      
      names(cov_x) <- cov_names
      
      cov_x2 <- subset(cov_x, cov_selected)
      
      tilename_x <- basename(subdir_tiles[x])
      
      outname_x <- dir_pred_tiles_frac %>%
        paste0(., "/", frac, "_", tilename_x, ".tif")
      
      predict(
        cov_x2,
        model_i,
        fun = predict_passna,
        na.rm = FALSE,
        filename = outname_x,
        overwrite = TRUE,
        const = data.frame(
          SOM_removed = TRUE,
          year = 2010
        )
      )
    }
  )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  
  outtiles_frac <- dir_pred_tiles_frac %>%
    list.files(full.names = TRUE) %>%
    sprc()
  
  merge(
    outtiles_frac,
    filename = paste0(dir_pred_all, frac, "_merged.tif"),
    overwrite = TRUE
  )
}

# March 16, 2023: 60 tiles, 121 predictors, clay: 32 hours
# March 22, 2023: 591 tiles, 121 predictors, clay: 3 h 48 min

# outfiles_table <- dir_pred_tiles_frac %>%
#   list.files(full.names = TRUE) %>%
#   file.info() %>%
#   rownames_to_column()
# 
# outfiles_dims <- dir_pred_tiles_frac %>%
#   list.files(full.names = TRUE) %>%
#   lapply(
#     function(x) {
#       r <- rast(x)
#       out <- dim(r)
#       return(out)
#     }
#   ) %>%
#   unlist() %>%
#   matrix(nrow = 3) %>%
#   t()
# 
# cbind(outfiles_table, outfiles_dims) %>%
#   write.table(
#     file = "out_tiles.csv",
#     sep = ";",
#     row.names = FALSE
#     )

# Old code without tiles:
# Maps for all of Denmark
# 2023-03-09: Took 24 hours for less than 25%. Not feasible.
# Need to test:
# - Variable selection
# - Tiles

# cov2 <- subset(cov, cov_selected)
# 
# predfolder2 <- paste0(dir_dat, "/predictions_", testn, "/") %T>% dir.create()
# 
# tmpfolder <- paste0(dir_dat, "/Temp/")
# 
# terraOptions(memfrac = 0.15, tempdir = tmpfolder)
# 
# maps <- list()
# 
# for(i in 1:length(fractions))
# {
#   frac <- fractions[i]
#   
#   showConnections()
#   
#   maps[[i]] <- predict(
#     cov2,
#     models[[i]],
#     fun = predict_passna,
#     na.rm = FALSE,
#     cores = 19,
#     filename = paste0(predfolder2, frac,  ".tif"),
#     overwrite = TRUE,
#     const = data.frame(
#       SOM_removed = TRUE,
#       year = 2010
#     )
#   )
# }


# Post processing
# Sum texture to 100
# transform SOC and CaCO3
# round values

# END