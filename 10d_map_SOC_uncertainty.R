# Map class uncertainties for SOC
# Only topsoil

library(parallel)
library(caret)
library(terra)
library(magrittr)
library(dplyr)
library(xgboost)
library(foreach)
library(stringr)
library(tools) # file_path_sans_ext

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/")

fractions_alt <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")

fractions <- fractions_alt

frac_ind_mineral <- c(1:4)
frac_ind_predict <- c(1:length(fractions))[-3]  # Exclude fine sand

fraction_names <- c(
  "Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3"
)

fraction_names_underscore <- c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand", "SOC", "CaCO3"
)

dir_cov <- dir_dat %>% paste0(., "/covariates")
cov_files <- dir_cov %>% list.files()
cov_names <- cov_files %>% tools::file_path_sans_ext()

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20231110.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

# cov_selected <- cov_cats %>%
#   filter(anbm_use == 1) %>%
#   dplyr::select(., name) %>%
#   unlist() %>%
#   unname()

source("f_predict_passna.R")
source("f_classify_SOC.R")

# Path for loading bootstrap models

dir_boot <- dir_results %>%
  paste0(., "/bootstrap/")

dir_boot_models <- dir_boot %>%
  paste0(., "/models/")

dir_boot_models_fractions <- dir_boot_models %>%
  paste0(., "/", fractions, "/")

models_boot_files <- lapply(
  1:length(fractions),
  function(x) {
    out <- fractions[x] %>%
      paste0(dir_boot_models, "/", ., "/") %>%
      list.files(full.names = TRUE)
    return(out)
  }
)

# Tiles for model prediction

numCores <- detectCores() - 1
numCores

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

subdir_tiles <- dir_tiles %>%
  list.dirs() %>%
  .[-1]

tilenames <- basename(subdir_tiles)

# Directory for saving bootstrap predictions

dir_pred_boot <- dir_boot %>%
  paste0(., "/predictions/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

dir_pred_tiles <- dir_pred_boot %>%
  paste0(., "/tiles/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Directory for temporary texture predictions, for standardization
# Overwrite for each depth + bootstrap repetition

dir_pred_tiles_tmp <- dir_pred_tiles %>%
  paste0(., "/tmp_fractions/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

for (x in 1:length(tilenames)) {
  dir_pred_tiles_tmp %>%
    paste0(., "/", tilenames[x], "/") %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)
}

dir_boot_sum_tiles <- dir_pred_boot %>%
  paste0(., "/summary/") %T>%
  dir.create(showWarnings = FALSE)

# Set number of bootstrap repetitions to use

nboot <- lapply(
  models_boot_files,
  function(x) {
    out <- length(x)
    return(out)
  }
) %>%
  unlist() %>%
  min()

nboot_final <- 100

nboot_max <- 100
# nboot_max <- 100

nboot <- min(c(nboot, nboot_max))

boot_all_chr <- c(1:nboot) %>%
  str_pad(
    .,
    nchar(nboot_final),
    pad = "0"
  ) %>%
  paste0("boot_", .)

# Set up loop for predicting each soil depth

n_digits <- 1

breaks <- c(0, 30, 60, 100, 200)

max_char <- breaks %>%
  as.character() %>%
  nchar() %>%
  max()

breaks_chr <- breaks %>%
  str_pad(
    .,
    max_char,
    pad = "0"
  )

j_all_depths <- 1:(length(breaks) - 1)
j_only_top <- 1

# only_top <- TRUE
# # only_top <- FALSE
# 
# if (only_top) {
#   j_depth <- j_only_top
# } else {
#   j_depth <- j_all_depths
# }

# j_depth <- 2 # Predict second layer
# j_depth <- 3 # Predict third layer
# j_depth <- 4 # Predict fourth layer

# Delete texture class predictions (force recalculation)

# dir_pred_tiles %>%
#   list.files(
#     pattern = "JB",
#     recursive = TRUE,
#     full.names = TRUE
#     ) %>%
#   file.remove()

# Execute loop for predicting each soil depth
for (j in j_depth) {
  breaks_j <- breaks[j:(j + 1)]
  breaks_j_chr <- breaks_chr[j:(j + 1)]
  print(paste0("Depth ", paste(breaks_j_chr, collapse = " - "), " cm"))
  dir_pred_tiles_depth <- dir_pred_tiles %>%
    paste0(
      ., "/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
    ) %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)
  
  # Loop for predicting maps from bootstrap repetitions
  for (bootr in 1:nboot) {
    bootr_chr <- bootr %>%
      str_pad(
        .,
        nchar(nboot_final),
        pad = "0"
      )
    print(paste0("Bootstrap repetition ", bootr_chr))
    
    dir_pred_tiles_bootr <- dir_pred_tiles_depth %>%
      paste0(., "/boot_", bootr_chr, "/") %T>%
      dir.create(showWarnings = FALSE, recursive = TRUE)
    
    for (x in 1:length(tilenames)) {
      dir_pred_tiles_bootr %>%
        paste0(., "/", tilenames[x], "/") %T>%
        dir.create(showWarnings = FALSE, recursive = TRUE)
    }
    
    n_outfiles_missing <- dir_pred_tiles_bootr %>%
      paste0(., "/", tilenames, "/") %>%
      rep(., each = length(fraction_names_underscore)) %>%
      paste0(., fraction_names_underscore, ".tif") %>%
      file.exists() %>%
      not() %>%
      sum()
    
    if (n_outfiles_missing > 0) {
      # for (i in frac_ind_predict) {
      for (i in 5) {  # Only SOC
        frac <- fraction_names_underscore[i]
        print(paste0("Mapping ", frac))
        
        model_i <- models_boot_files[[i]][bootr] %>% readRDS()
        
        cov_selected <- (varImp(model_i)$importance %>% row.names()) %>%
          .[. %in% cov_cats$name]
        
        showConnections()
        
        cl <- makeCluster(numCores)
        
        clusterEvalQ(
          cl,
          {
            library(terra)
            library(caret)
            library(xgboost)
            library(magrittr)
            library(dplyr)
            library(tools)
          }
        )
        
        clusterExport(
          cl,
          c(
            "i",
            "model_i",
            "subdir_tiles",
            "dir_pred_tiles_bootr",
            "dir_pred_tiles_tmp",
            "frac",
            "cov_selected",
            "predict_passna",
            "dir_dat",
            "n_digits",
            "breaks_j",
            "breaks_j_chr",
            "frac_ind_mineral",
            "classify_SOC"
          )
        )
        
        parSapplyLB(
          cl,
          1:length(subdir_tiles),
          function(x) {
            tilename_x <- basename(subdir_tiles[x])
            
            outname_x_final <- dir_pred_tiles_bootr %>%
              paste0(., "/", tilename_x, "/", frac, ".tif")
            
            const_i <- data.frame(
              upper = breaks_j[1],
              lower = breaks_j[2]
            )
            
            # Write temporary files for the mineral fractions
            if (i %in% frac_ind_mineral) {
              outname_x_tmp <- dir_pred_tiles_tmp %>%
                paste0(., "/", tilename_x, "/", frac, ".tif")

              outname_x <- outname_x_tmp

              const_i$SOM_removed <- 1
            } else {
              outname_x <- outname_x_final
            }
            
            # Omit predictions if the output file already exists
            if (!file.exists(outname_x_final)) {
              tmpfolder <- paste0(dir_dat, "/Temp/")
              
              terraOptions(memfrac = 0.02, tempdir = tmpfolder)
              
              cov_x_files <- subdir_tiles[x] %>%
                list.files(full.names = TRUE)
              
              cov_x_names <- cov_x_files %>%
                basename() %>%
                tools::file_path_sans_ext()
              
              cov_x <- cov_x_files %>% rast()
              
              names(cov_x) <- cov_x_names
              
              cov_x %<>% terra::subset(., cov_selected)
              
              predict(
                cov_x,
                model_i,
                fun = predict_passna,
                na.rm = FALSE,
                const = const_i,
                n_const = ncol(const_i),
                n_digits = 1,
                filename = outname_x,
                overwrite = TRUE
              )
              
              # Calculate SOC class for each bootstrap repetition
              
              rs_s2 <- outname_x %>% rast()
              
              names(rs_s2) <- c("SOC")
              
              outname_x2_SOC_class <- dir_pred_tiles_bootr %>%
                paste0(., "/", tilename_x, "/SOC_class.tif")
              
              lapp(
                rs_s2,
                classify_SOC,
                filename = outname_x2,
                overwrite = TRUE,
                wopt = list(
                  datatype = "INT1U",
                  NAflag = 13
                )
              )
            }
          
            return(NULL)
          }
        )
        
        stopCluster(cl)
        foreach::registerDoSEQ()
        rm(cl)
      }
      
      
      
      # Calculate the probability of each SOC class across repetitions
      
      # Merge SOC class probability rasters
      
      # # Standardize mineral texture predictions
      # print("Standardizing mineral fractions")
      # dir_pred_tiles_100 <- paste0(dir_pred_tiles, "/tmp_tex100") %T>%
      #   dir.create(showWarnings = FALSE, recursive = TRUE)
      # for (x in 1:length(tilenames)) {
      #   dir_pred_tiles_100 %>%
      #     paste0(., "/", tilenames[x], "/") %T>%
      #     dir.create(showWarnings = FALSE, recursive = TRUE)
      # }
      # showConnections()
      # cl <- makeCluster(numCores)
      # clusterEvalQ(
      #   cl,
      #   {
      #     library(terra)
      #     library(magrittr)
      #     library(dplyr)
      #     library(tools)
      #   }
      # )
      # clusterExport(
      #   cl,
      #   c(
      #     "subdir_tiles",
      #     "dir_dat",
      #     "n_digits",
      #     "dir_pred_tiles_100",
      #     "fraction_names_underscore",
      #     "breaks_j_chr",
      #     "dir_pred_tiles_tmp",
      #     "dir_pred_tiles_bootr"
      #   )
      # )
      # parSapplyLB(
      #   cl,
      #   1:length(subdir_tiles),
      #   function(x) {
      #     tilename_x <- basename(subdir_tiles[x])
      #     n_rs_tex <- dir_pred_tiles_bootr %>%
      #       paste0(., "/", tilename_x, "/") %>%
      #       list.files(pattern = ".tif", full.names = TRUE) %>%
      #       length()
      #     
      #     if (n_rs_tex < 6) {
      #       tmpfolder <- paste0(dir_dat, "/Temp/")
      #       
      #       terraOptions(memfrac = 0.02, tempdir = tmpfolder)
      #       
      #       tmp_files_tile_x <- dir_pred_tiles_tmp %>%
      #         paste0(., "/", tilename_x, "/") %>% 
      #         list.files(full.names = TRUE)
      #       
      #       rs_pred_tile_x <- tmp_files_tile_x %>% rast()
      #       
      #       layernames <- tmp_files_tile_x %>%
      #         tools::file_path_sans_ext() %>%
      #         basename() %>%
      #         c(., "Fine_sand")
      #       
      #       outname_x <- dir_pred_tiles_100 %>%
      #         paste0(., "/", tilename_x, ".tif")
      #       
      #       app(
      #         rs_pred_tile_x,
      #         function(r) {
      #           res <- max(c(100 - sum(r), 0))
      #           out <- c(r, res)
      #           out <- out * 100 / sum(out)
      #           out %<>% round(digits = n_digits)
      #           names(out) <- layernames
      #           return(out)
      #         },
      #         overwrite = TRUE,
      #         filename = outname_x
      #       )
      #       
      #       tex100_x <- outname_x %>% rast()
      #       
      #       for (i in 1:length(layernames)) {
      #         outname_x_final <- dir_pred_tiles_bootr %>%
      #           paste0(., "/", tilename_x, "/", layernames[i], ".tif")
      #         
      #         writeRaster(
      #           tex100_x[[i]],
      #           filename = outname_x_final,
      #           overwrite = TRUE
      #         )
      #       }
      #     }
      #     return(NULL)
      #   }
      # )
      # stopCluster(cl)
      # foreach::registerDoSEQ()
      # rm(cl)
    }
    
    # n_JB_outfiles_missing <- dir_pred_tiles_bootr %>%
    #   paste0(., "/", tilenames, "/") %>%
    #   paste0(., "JB.tif") %>%
    #   file.exists() %>%
    #   not() %>%
    #   sum()
    # 
    # if (n_JB_outfiles_missing > 0) {
    #   # Classify JB
    #   print("Calculating soil texture class")
    #   source("f_classify_soil_JB.R")
    #   
    #   showConnections()
    #   
    #   cl <- makeCluster(numCores)
    #   
    #   clusterEvalQ(
    #     cl,
    #     {
    #       library(terra)
    #       library(magrittr)
    #       library(dplyr)
    #       library(caret)
    #       library(xgboost)
    #       library(tools)
    #     }
    #   )
    #   
    #   clusterExport(
    #     cl,
    #     c(
    #       "breaks_j_chr",
    #       "subdir_tiles",
    #       "dir_dat",
    #       "dir_pred_tiles_bootr",
    #       "fraction_names_underscore",
    #       "classify_soil_JB"
    #     )
    #   )
    #   
    #   parSapplyLB(
    #     cl,
    #     1:length(subdir_tiles),
    #     function(x) {
    #       tilename_x <- basename(subdir_tiles[x])
    #       
    #       outname_x <- dir_pred_tiles_bootr %>%
    #         paste0(., "/", tilename_x, "/JB.tif")
    #       
    #       if (!file.exists(outname_x)) {
    #         tmpfolder <- paste0(dir_dat, "/Temp/")
    #         
    #         terraOptions(memfrac = 0.02, tempdir = tmpfolder)
    #         
    #         rs_tex <- dir_pred_tiles_bootr %>%
    #           paste0(
    #             ., "/", tilename_x, "/",
    #             fraction_names_underscore,
    #             ".tif"
    #           ) %>%
    #           rast()
    #         
    #         rs_s2 <- terra::subset(rs_tex, c(1, 2, 3, 5, 6))
    #         
    #         names(rs_s2) <- c("clay", "silt", "sand_f", "SOM", "CaCO3")
    #         
    #         lapp(
    #           rs_s2,
    #           classify_soil_JB,
    #           # SOM_factor = 1 / 0.587,
    #           SOM_factor = 1 / 0.6,
    #           filename = outname_x,
    #           overwrite = TRUE,
    #           wopt = list(
    #             datatype = "INT1U",
    #             NAflag = 13
    #           )
    #         )
    #         return(NULL)
    #       }
    #     }
    #   )
    #   stopCluster(cl)
    #   foreach::registerDoSEQ()
    #   rm(cl)
    # }
  }
  
  # Summarize results across bootstrap repetitions
  dir_sum_depth_tiles <- dir_boot_sum_tiles %>%
    paste0(
      ., "/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
    ) %T>%
    dir.create(showWarnings = FALSE)
  
  for (x in 1:length(tilenames)) {
    dir_sum_depth_tiles %>%
      paste0(., "/", tilenames[x], "/") %T>%
      dir.create(showWarnings = FALSE, recursive = TRUE)
  }
  
  # for (i in 1:length(fraction_names_underscore)) {
  #   frac <- fraction_names_underscore[i]
  #   print(paste0("Summarizing ", frac))
  #   
  #   showConnections()
  #   
  #   cl <- makeCluster(numCores)
  #   
  #   clusterEvalQ(
  #     cl,
  #     {
  #       library(terra)
  #       library(magrittr)
  #       library(dplyr)
  #       library(tools)
  #     }
  #   )
  #   
  #   clusterExport(
  #     cl,
  #     c(
  #       "boot_all_chr",
  #       "subdir_tiles",
  #       "frac",
  #       "dir_dat",
  #       "n_digits",
  #       "dir_sum_depth_tiles",
  #       "dir_pred_tiles_depth"
  #     )
  #   )
  #   
  #   parSapplyLB(
  #     cl,
  #     1:length(subdir_tiles),
  #     function(x) {
  #       tilename_x <- basename(subdir_tiles[x])
  #       
  #       r_frac_tile <- paste0(
  #         dir_pred_tiles_depth, "/", boot_all_chr, "/", tilename_x, "/",
  #         frac, ".tif") %>%
  #         rast()
  #       
  #       outname_x_mean <- dir_sum_depth_tiles %>%
  #         paste0(., "/", tilename_x, "/", frac, "_mean.tif")
  #       outname_x_sd <- dir_sum_depth_tiles %>%
  #         paste0(., "/", tilename_x, "/", frac, "_sd.tif")
  #       outname_x_q05 <- dir_sum_depth_tiles %>%
  #         paste0(., "/", tilename_x, "/", frac, "_q05.tif")
  #       outname_x_q95 <- dir_sum_depth_tiles %>%
  #         paste0(., "/", tilename_x, "/", frac, "_q95.tif")
  #       
  #       tmpfolder <- paste0(dir_dat, "/Temp/")
  #       
  #       terraOptions(memfrac = 0.02, tempdir = tmpfolder)
  #       
  #       # Calculate mean
  #       app(
  #         r_frac_tile,
  #         fun = mean
  #       ) %>%
  #         math(
  #           fun = "round",
  #           digits = n_digits,
  #           filename = outname_x_mean,
  #           overwrite = TRUE
  #         )
  #       # Calculate standard deviation
  #       app(
  #         r_frac_tile,
  #         fun = sd
  #       ) %>%
  #         math(
  #           fun = "round",
  #           digits = n_digits,
  #           filename = outname_x_sd,
  #           overwrite = TRUE
  #         )
  #       # Calculate 5% quantile
  #       terra::quantile(
  #         r_frac_tile,
  #         probs = 0.05,
  #         na.rm = TRUE,
  #         filename = outname_x_q05,
  #         overwrite = TRUE
  #       )
  #       # Calculate 95% quantile
  #       terra::quantile(
  #         r_frac_tile,
  #         probs = 0.95,
  #         na.rm = TRUE,
  #         filename = outname_x_q95,
  #         overwrite = TRUE
  #       )
  #       return(NULL)
  #     }
  #   )
  #   
  #   stopCluster(cl)
  #   foreach::registerDoSEQ()
  #   rm(cl)
  # }
  
  # Calculate mean soil class and soil class uncertainties
  
  # print("Summarizing texture classes")
  # source("f_classify_soil_JB.R")
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
  #   c(
  #     "boot_all_chr",
  #     "subdir_tiles",
  #     "dir_dat",
  #     "n_digits",
  #     "breaks_j",
  #     "breaks_j_chr",
  #     "dir_sum_depth_tiles",
  #     "fraction_names_underscore",
  #     "nboot",
  #     "dir_pred_tiles_depth",
  #     "classify_soil_JB"
  #   )
  # )
  # 
  # parSapplyLB(
  #   cl,
  #   1:length(subdir_tiles),
  #   function(x) {
  #     tilename_x <- basename(subdir_tiles[x])
  #     
  #     r_frac_tile_means <- paste0(
  #       dir_sum_depth_tiles, "/", tilename_x, "/",
  #       fraction_names_underscore, "_mean.tif") %>%
  #       rast()
  #     
  #     r_JB_tile_all <- paste0(
  #       dir_pred_tiles_depth, "/", boot_all_chr, "/", tilename_x, "/JB.tif") %>%
  #       rast()
  #     
  #     n_JB_lyr <- nlyr(r_JB_tile_all)
  #     
  #     outname_x_JB_mean <- dir_sum_depth_tiles %>%
  #       paste0(., "/", tilename_x, "/JB_mean.tif")
  #     outname_x_JB_unc <- dir_sum_depth_tiles %>%
  #       paste0(., "/", tilename_x, "/JB_unc.tif")
  #     
  #     tmpfolder <- paste0(dir_dat, "/Temp/")
  #     
  #     terraOptions(memfrac = 0.02, tempdir = tmpfolder)
  #     
  #     r_frac_tile_means <- terra::subset(r_frac_tile_means, c(1, 2, 3, 5, 6))
  #     
  #     names(r_frac_tile_means) <- c("clay", "silt", "sand_f", "SOM", "CaCO3")
  #     
  #     # JB texture class for mean texture
  #     lapp(
  #       r_frac_tile_means,
  #       classify_soil_JB,
  #       # SOM_factor = 1 / 0.587,
  #       SOM_factor = 1 / 0.6,
  #       filename = outname_x_JB_mean,
  #       overwrite = TRUE,
  #       wopt = list(
  #         datatype = "INT1U",
  #         NAflag = 13
  #       )
  #     )
  #     # Calculate JB uncertainty
  #     mode_JB <- terra::modal(r_JB_tile_all, na.rm = TRUE)
  #     terra::compare(r_JB_tile_all, mode_JB, "!=") %>%
  #       app(
  #         fun = sum,
  #         na.rm = TRUE,
  #         filename = outname_x_JB_unc,
  #         overwrite = TRUE,
  #         wopt = list(
  #           datatype = "INT1U",
  #           NAflag = 101
  #         )
  #       )
  #     
  #     return(NULL)
  #   }
  # )
  # 
  # stopCluster(cl)
  # foreach::registerDoSEQ()
  # rm(cl)
  
  # Merge summary rasters
  
  outfiles_basenames <- paste0(
    dir_sum_depth_tiles, "/", tilenames[1], "/") %>%
    list.files() %>%
    basename() %>%
    tools::file_path_sans_ext()
  
  dir_merged_depth <- dir_pred_boot %>%
    paste0(
      ., "/final_maps/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
    ) %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)
  
  for (i in 1:length(outfiles_basenames)) {
    summary_tiles_i <- paste0(
      dir_sum_depth_tiles, "/", tilenames, "/", outfiles_basenames[i], ".tif"
    )
    tile1_i <- summary_tiles_i[1] %>% rast()
    
    dtyp_i <- datatype(tile1_i)
    naflag_i <- NAflag(tile1_i)
    
    outtiles_sprc <- summary_tiles_i %>% sprc()
    
    merge(
      outtiles_sprc,
      filename = paste0(dir_merged_depth, "/", outfiles_basenames[i], ".tif"),
      overwrite = TRUE,
      gdal = "TILED=YES",
      datatype = dtyp_i,
      NAflag = naflag_i,
      names = outfiles_basenames[i]
    )
  }
}


# END