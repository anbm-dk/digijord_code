# 10c: Map the effect of SOM removal on texture predictions

# Procedure:
# For depth 1:
## For each tile:
### For each bootstrap iteration:
#### Predict texture fractions with and without som removal (tmp)
#### Standardize texture predictions (tmp)
#### Calculate differences for mineral texture fractions (tmp)
#### Predict SOC and CaCO3 only once (tmp)
#### Map soil texture class with and without som removal (tmp)
#### Calculate differences for texture classes (tmp)
### Summarize texture without SOM removal (save)
### Summarize only mean value for SOC and CaCO3 (save)
### Summarize texture classes without SOM removal (save)
### Summarize differences for texture with and without som removal (save)
### Summarize differences for texture classes (save)
## Merge maps of summarized mineral fractions classes (save)


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
frac_ind_mineral_predict <- frac_ind_mineral[frac_ind_mineral %in% frac_ind_predict]
frac_ind_nonmineral <- frac_ind_predict[!frac_ind_predict %in% frac_ind_mineral]

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

# Names for SOM removal

names_somremoval <- c("norem", "somrem")
numbers_somremoval <- 0:1

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

nboot_max <- 19
# nboot_max <- 100

nboot <- min(c(nboot, nboot_max))

boot_all_chr <- c(1:nboot) %>%
  str_pad(
    .,
    nchar(nboot_final),
    pad = "0"
  ) %>%
  paste0("boot_", .)


# Tiles for model prediction

numCores <- detectCores() - 1
numCores

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_591/")

subdir_tiles <- dir_tiles %>%
  list.dirs() %>%
  .[-1]

tilenames <- basename(subdir_tiles)

n_tiles_use <- 3
# n_tiles_use <- length(subdir_tiles)

if (n_tiles_use < length(subdir_tiles)) {
  tilenames <- tilenames[1:n_tiles_use]
}

# Directory for saving bootstrap predictions
dir_pred_diff <- dir_boot %>%
  paste0(., "/predictions_difference/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Folder for temporary bootstrap predictions for the current tile
dir_pred_tileboot_tmp <- dir_pred_diff %>%
  paste0(., "/tile_x_tmp/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Folder for non-standardized mineral fraction predictions:
# Subfolders (3 files each):
# - with SOM removed
# - without SOM removal

dir_pred_tileboot_raw <- dir_pred_tileboot_tmp %>%
  paste0(., "/mineral_raw/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

dir_pred_tileboot_raw_subdirs <- dir_pred_tileboot_raw %>%
  paste0(., "/", names_somremoval, "/") %T>%
  sapply(
    function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    }
  )

# Subfolders for bootstrap predictions
dir_pred_tileboot_raw_subdirs %>%
  sapply(., function(x2) {
    paste0(x2, "/", boot_all_chr) %>%
      sapply(
        function(x) {
          dir.create(x, showWarnings = FALSE, recursive = TRUE)
        }
      )
  }
  )

dir_pred_tex100 <- dir_pred_tileboot_tmp %>%
  paste0(., "/tex100/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Folder for standardized predictions:
# Files:
# - mineral texture fractions (2x4 files):
# -- with SOM removed
# -- without SOM removal
# - mineral fraction differences (4 files)
# - SOC and CaCO3 (2 files)
# - texture class (2 files)
# -- with SOM removed
# -- without SOM removal
# - texture class difference (1 file)

dir_pred_tileboot_standard <- dir_pred_tileboot_tmp %>%
  paste0(., "/standardized/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Subfolders for bootstrap predictions
dir_pred_tileboot_standard %>%
  paste0(., "/", boot_all_chr) %>%
  sapply(
    function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    }
  )

# Folder for tile summaries
dir_pred_tiles_sum <- dir_pred_diff %>%
  paste0(., "/summarized/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

# Subfolder for each tile
dir_pred_tiles_sum %>%
  paste0(., "/", tilenames) %>%
  sapply(
    function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    }
  )

# Files:
# - Mean and SD for each mineral fraction, without SOM removal (2x4 = 8 files)
# - Mean and SD for differences (2 x 4 files = 8 files)
# - Mean for SOC and CaCO3 (2 files)
# - Texture class for mean texture without SOM removal (1 file)
# - Uncertainty for texture class without SOM removal (1 file)
# - Sum of texture class differences (1 file)

names_summary <- c(
  paste0(fraction_names_underscore[frac_ind_mineral], "_norem_mean.tif"),
  paste0(fraction_names_underscore[frac_ind_mineral], "_norem_SD.tif"),
  paste0(fraction_names_underscore[frac_ind_mineral], "_diff_mean.tif"),
  paste0(fraction_names_underscore[frac_ind_mineral], "_diff_SD.tif"),
  paste0(fraction_names_underscore[5:6], "_mean.tif"),
  "JB_norem_mean.tif",
  "JB_norem_unc.tif",
  "JB_diff_mean.tif"
)
names_summary

# Merge tiles:
# Files:
# - Mean and SD for each mineral fraction, without SOM removal (2x4 = 8 files)
# - Mean and SD for differences (2 x 4 files = 8 files)
# - Texture class for mean texture without SOM removal (1 file)
# - Uncertainty for texture class without SOM removal (1 file)
# - Sum of texture class differences (1 file)

dir_pred_merged <- dir_pred_diff %>%
  paste0(., "/final/") %T>%
  dir.create(showWarnings = FALSE, recursive = TRUE)

names_merge <- grep("SOC|CaCO3", names_summary, value = TRUE, invert = TRUE)
names_merge

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

only_top <- TRUE
# # only_top <- FALSE

if (only_top) {
  j_depth <- j_only_top
} else {
  j_depth <- j_all_depths
}

# Execute loop for predicting each soil depth
for (j in j_depth) {
  breaks_j <- breaks[j:(j + 1)]
  breaks_j_chr <- breaks_chr[j:(j + 1)]
  print(paste0("Depth ", paste(breaks_j_chr, collapse = " - "), " cm"))
  
  # Loop for predicting tiles
  for (x in 1:length(tilenames)) {
    tilename_x <- basename(subdir_tiles[x])
    
    print(paste0("Mapping ", tilename_x))
    
    n_outfiles_missing <- dir_pred_tiles_sum %>%
      paste0(., "/", tilename_x, "/") %>%
      rep(., each = length(names_summary)) %>%
      paste0(., names_summary, ".tif") %>%
      file.exists() %>%
      not() %>%
      sum()
    
    if (n_outfiles_missing > 0) {
# Function for parallel bootstrap prediction
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
          "boot_all_chr",
          "x",
          "subdir_tiles",
          "tilename_x",
          "dir_pred_diff",
          "dir_pred_tileboot_raw",
          "dir_pred_tileboot_raw_subdirs",
          "dir_pred_tileboot_standard",
          "dir_pred_tileboot_tmp",
          "dir_pred_tex100",
          "predict_passna",
          "dir_dat",
          "n_digits",
          "breaks_j",
          "breaks_j_chr",
          "frac_ind_mineral",
          "frac_ind_mineral_predict",
          "frac_ind_nonmineral",
          "fraction_names_underscore",
          "models_boot_files",
          "cov_cats",
          "names_somremoval",
          "numbers_somremoval"
        )
      )
      
      parSapplyLB(
        cl,
        1:length(boot_all_chr),
        function(bootr) {
          name_bootr <- boot_all_chr[bootr]
          
          tmpfolder <- paste0(dir_dat, "/Temp/")
          
          terraOptions(memfrac = 0.02, tempdir = tmpfolder)
          
          # Load covariates
          cov_x_files <- subdir_tiles[x] %>%
            list.files(full.names = TRUE)
          
          cov_x_names <- cov_x_files %>%
            basename() %>%
            tools::file_path_sans_ext()
          
          cov_x <- cov_x_files %>% rast()
          
          names(cov_x) <- cov_x_names
          
          # Predict mineral fractions
          # - norem
          # - somrem
          for (i in frac_ind_mineral_predict) {
            frac <- fraction_names_underscore[i]
            
            model_i <- models_boot_files[[i]][bootr] %>% readRDS()
            
            cov_selected <- (varImp(model_i)$importance %>% row.names()) %>%
              .[. %in% cov_cats$name]
            
            for (m in 1:length(names_somremoval)) {
              outname_x <- dir_pred_tileboot_raw_subdirs[m] %>%
                paste0(., "/", name_bootr, "/", frac, ".tif")
              
              const_i <- data.frame(
                upper = breaks_j[1],
                lower = breaks_j[2],
                SOM_removed = numbers_somremoval[m]
              )
              
              cov_i <- cov_x %>% terra::subset(., cov_selected)
              
              predict(
                cov_i,
                model_i,
                fun = predict_passna,
                na.rm = FALSE,
                const = const_i,
                n_const = ncol(const_i),
                n_digits = 1,
                filename = outname_x,
                overwrite = TRUE
              )
            }
          }
          # standardize mineral fractions
          # - norem
          # - somrem
          for (m in 1:length(names_somremoval)) {
            tmp_files_bootr <- dir_pred_tileboot_raw_subdirs[m] %>%
              paste0(., "/", name_bootr, "/") %>% 
              list.files(full.names = TRUE)
            
            rs_pred_bootr <- tmp_files_bootr %>% rast()
            
            layernames <- tmp_files_bootr %>%
              tools::file_path_sans_ext() %>%
              basename() %>%
              c(., "Fine_sand")
            
            outname_x <- dir_pred_tex100 %>%
              paste0(., "/", name_bootr, "_", names_somremoval[m], ".tif")
            
            app(
              rs_pred_bootr,
              function(r) {
                res <- max(c(100 - sum(r), 0))
                out <- c(r, res)
                out <- out * 100 / sum(out)
                out %<>% round(digits = n_digits)
                names(out) <- layernames
                return(out)
              },
              overwrite = TRUE,
              filename = outname_x
            )
            
            tex100_x <- outname_x %>% rast()
            
            for (i in 1:length(layernames)) {
              outname_x_final <- dir_pred_tileboot_standard %>%
                paste0(
                  ., "/", name_bootr, "/", layernames[i], "_",
                  names_somremoval[m], ".tif"
                )
              
              writeRaster(
                tex100_x[[i]],
                filename = outname_x_final,
                overwrite = TRUE
              )
            }
          }
          # Predict SOC and CaCO3
          for (i in frac_ind_nonmineral) {
            frac <- fraction_names_underscore[i]
            
            model_i <- models_boot_files[[i]][bootr] %>% readRDS()
            
            cov_selected <- (varImp(model_i)$importance %>% row.names()) %>%
              .[. %in% cov_cats$name]
            
            
            outname_x <- dir_pred_tileboot_standard %>%
              paste0(., "/", name_bootr, "/", frac, ".tif")
            
            const_i <- data.frame(
              upper = breaks_j[1],
              lower = breaks_j[2]
            )
            
            cov_i <- cov_x %>% terra::subset(., cov_selected)
            
            predict(
              cov_i,
              model_i,
              fun = predict_passna,
              na.rm = FALSE,
              const = const_i,
              n_const = ncol(const_i),
              n_digits = 1,
              filename = outname_x,
              overwrite = TRUE
            )
          }

          # Calculate JB
          # - somrem
          # - norem
          source("f_classify_soil_JB.R")
          
          for (m in 1:length(names_somremoval)) {
            outname_x <- dir_pred_tileboot_standard %>%
              paste0(., "/", name_bootr, "/JB_", names_somremoval[m], ".tif")
            
            rs_tex <- c(
              paste0(
                dir_pred_tileboot_standard, "/", name_bootr, "/",
                fraction_names_underscore[frac_ind_mineral],
                "_", names_somremoval[m],
                ".tif"
              ),
              paste0(
                dir_pred_tileboot_standard, "/", name_bootr, "/",
                fraction_names_underscore[frac_ind_nonmineral],
                ".tif"
              )
            ) %>%
              rast()
            
            rs_s2 <- terra::subset(rs_tex, c(1, 2, 3, 5, 6))
            names(rs_s2) <- c("clay", "silt", "sand_f", "SOM", "CaCO3")
            
            lapp(
              rs_s2,
              classify_soil_JB,
              # SOM_factor = 1 / 0.587,
              # SOM_factor = 1 / 0.6,
              SOM_factor = 0,  # NB! Disregard SOM as, JB11 depends on Peat2022
              filename = outname_x,
              overwrite = TRUE,
              wopt = list(
                datatype = "INT1U",
                NAflag = 13
              )
            )
          }
          
          # Calculate mineral fraction differences
          for (i in frac_ind_mineral) {
            r1 <- dir_pred_tileboot_standard %>%
              paste0(
                ., "/", name_bootr, "/", fraction_names_underscore[i], "_",
                names_somremoval[1], ".tif"
              ) %>% rast()
            r2 <- dir_pred_tileboot_standard %>%
              paste0(
                ., "/", name_bootr, "/", fraction_names_underscore[i], "_",
                names_somremoval[2], ".tif"
              ) %>% rast()
            
            r_diff <- r2 - r1
            
            outname <- dir_pred_tileboot_standard %>%
              paste0(
                ., "/", name_bootr, "/", fraction_names_underscore[i],
                "_diff.tif"
              )
            
            writeRaster(
              r_diff,
              filename = outname,
              overwrite = TRUE
            )
          }
          
          # Calculate JB differences
          r1 <- dir_pred_tileboot_standard %>%
            paste0(., "/", name_bootr, "/JB_", names_somremoval[1], ".tif") %>%
            rast()
          r2 <- dir_pred_tileboot_standard %>%
            paste0(., "/", name_bootr, "/JB_", names_somremoval[2], ".tif") %>%
            rast()
          
          r_diff <- r1 != r2
          
          writeRaster(
            r_diff,
            filename = paste0(
              dir_pred_tileboot_standard,
              "/", name_bootr, "/JB_diff.tif"
            ),
            overwrite = TRUE,
            datatype = "INT1U"
          )

          return(NULL)
        }
      )
      
      stopCluster(cl)
      foreach::registerDoSEQ()
      rm(cl)
      
      # Summarize results across bootstrap repetitions
      
      names_summary_num <- names_summary %>%
        grep("JB", .,  invert = TRUE, value = TRUE)

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
          library(stringr)
        }
      )
      
      clusterExport(
        cl,
        c(
          "boot_all_chr",
          "x",
          "tilename_x",
          "dir_pred_diff",
          "dir_pred_tileboot_raw",
          "dir_pred_tileboot_raw_subdirs",
          "dir_pred_tileboot_standard",
          "dir_pred_tileboot_tmp",
          "dir_dat",
          "n_digits",
          "breaks_j",
          "breaks_j_chr",
          "frac_ind_mineral",
          "frac_ind_mineral_predict",
          "frac_ind_nonmineral",
          "fraction_names_underscore",
          "names_summary_num",
          "dir_pred_tiles_sum"
        )
      )
      
      parSapplyLB(
        cl,
        1:length(names_summary_num),
        function(i) {
          tmpfolder <- paste0(dir_dat, "/Temp/")
          terraOptions(memfrac = 0.02, tempdir = tmpfolder)
          
          splitname <- names_summary_num[i] %>%
            file_path_sans_ext() %>%
            str_split(., "_") %>%
            unlist()
          
          stat_i <- splitname[length(splitname)]
          input_name <- splitname[-length(splitname)] %>% paste(collapse = "_")
          
          input_layers_i <- dir_pred_tileboot_standard %>%
            list.files(recursive = TRUE, full.names = TRUE) %>%
            grep(input_name, ., value = TRUE) %>%
            rast()
          
          if (stat_i == "mean") {
            outname <- dir_pred_tiles_sum %>%
              paste0(., "/", tilename_x, "/", input_name, "_mean.tif")
            
            app(
              input_layers_i,
              fun = mean
            ) %>%
              math(
                fun = "round",
                digits = n_digits,
                filename = outname,
                overwrite = TRUE
              )
          }
          
          if (stat_i == "SD") {
            outname <- dir_pred_tiles_sum %>%
              paste0(., "/", tilename_x, "/", input_name, "_sd.tif")
            
            app(
              input_layers_i,
              fun = sd
            ) %>%
              math(
                fun = "round",
                digits = n_digits,
                filename = outname,
                overwrite = TRUE
              )
          }
          
          return(NULL)
        }
      )
      
      stopCluster(cl)
      foreach::registerDoSEQ()
      rm(cl)
      
      # Texture class for mean texture predictions without SOM removal
      
      source("f_classify_soil_JB.R")
      
      # Summarize texture classes without SOM removal
      
      r_frac_tile_means <- c(
        paste0(
          dir_pred_tiles_sum, "/", tilename_x, "/",
          fraction_names_underscore[frac_ind_mineral],
          "_", names_somremoval[1], "_mean",
          ".tif"
        ),
        paste0(
          dir_pred_tiles_sum, "/", tilename_x, "/",
          fraction_names_underscore[frac_ind_nonmineral], "_mean",
          ".tif"
        )
      ) %>%
        rast()
      
      r_JB_tile_all <- paste0(
        dir_pred_tileboot_standard, "/", boot_all_chr, "/",
        "/JB_norem.tif") %>%
        rast()
      
      n_JB_lyr <- nlyr(r_JB_tile_all)
      
      outname_x_JB_mean <- dir_pred_tiles_sum %>%
        paste0(., "/", tilename_x, "/JB_norem_mean.tif")
      outname_x_JB_unc <- dir_pred_tiles_sum %>%
        paste0(., "/", tilename_x, "/JB_norem_unc.tif")
      
      tmpfolder <- paste0(dir_dat, "/Temp/")
      
      terraOptions(tempdir = tmpfolder)
      
      r_frac_tile_means <- terra::subset(r_frac_tile_means, c(1, 2, 3, 5, 6))
      
      names(r_frac_tile_means) <- c("clay", "silt", "sand_f", "SOM", "CaCO3")
      
      # JB texture class for mean texture
      lapp(
        r_frac_tile_means,
        classify_soil_JB,
        # SOM_factor = 1 / 0.587,
        SOM_factor = 0,
        filename = outname_x_JB_mean,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT1U",
          NAflag = 13
        )
      )
      # Calculate JB uncertainty
      mode_JB <- terra::modal(r_JB_tile_all, na.rm = TRUE)
      terra::compare(r_JB_tile_all, mode_JB, "!=") %>%
        app(
          fun = sum,
          na.rm = TRUE,
          filename = outname_x_JB_unc,
          overwrite = TRUE,
          wopt = list(
            datatype = "INT1U",
            NAflag = 101
          )
        )
      
      # Summarize differences for texture classes
      r_JB_diff_all <- paste0(
        dir_pred_tileboot_standard, "/", boot_all_chr, "/",
        "/JB_diff.tif") %>%
        rast()
      
      outname_x_JB_diff <- dir_pred_tiles_sum %>%
        paste0(., "/", tilename_x, "/JB_diff_mean.tif")
      
      app(
        r_JB_diff_all,
        fun = sum,
        na.rm = TRUE,
        filename = outname_x_JB_diff,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT1U",
          NAflag = 101
        )
      )
    }
  }
  # Merge maps of summarized mineral fractions classes
  dir_pred_merged_depth <- dir_pred_merged %>%
    paste0(
      ., "/depth_", breaks_j_chr[1], "_", breaks_j_chr[2], "_cm/"
    ) %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)
  
  for (i in 1:length(names_merge)) {
    summary_tiles_i <- paste0(
      dir_pred_tiles_sum, "/", tilenames, "/", names_merge[i]
    )
    tile1_i <- summary_tiles_i[1] %>% rast()
    
    dtyp_i <- datatype(tile1_i)
    naflag_i <- NAflag(tile1_i)
    
    outtiles_sprc <- summary_tiles_i %>% sprc()
    
    outfile_basename <- names_merge[i] %>% file_path_sans_ext()
    
    merge(
      outtiles_sprc,
      filename = paste0(dir_pred_merged_depth, "/", names_merge[i]),
      overwrite = TRUE,
      gdal = "TILED=YES",
      datatype = dtyp_i,
      NAflag = naflag_i,
      names = outfile_basename
    )
  }
}


# END