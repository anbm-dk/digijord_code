# First maps for the topsoil (0 - 20 cm)
# Clay and SOM
# Only bare soil area
# Only DSC points (so far)
# Boosted cubist
# k fold CV
# Test all covariates
# Make maps with the best stack
# Make new covariate extracts for the points that I use in the model
# Get proper texture values from the access database (not rounded off)

library(terra)  # NB!
library(tidyverse)
library(magrittr)
library(stringr)

testn <- 3

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()
cov_dir <- root %>% paste0(., '/covariates')

# 1: Move bare soil composite
# - remove zeroes
# - don't include count layer

# reclasser <- matrix(c(0, NA), ncol = 2)
# 
# bare_dir <- root %>% paste0(., '/layers/bare-soil')
# bare_files <- bare_dir %>% list.files(pattern = "tif$"
#                                       , full.names = TRUE) %>%
#   head(., -1)
# for(i in 1:length(bare_files))
# {
#   f_i <- bare_files[[i]]
#   outname <- f_i %>% basename %>% paste0(cov_dir, '/', .)
#   r <- f_i %>%
#     rast %>%
#     classify(rcl = reclasser
#              , filename = outname
#              )
# }

# 2: Move newest vegetated composites
# - remove zeroes
# - don't include count layers

# veg_dir <- root %>% paste0(., '/layers/sentinel_composites')
# veg_files <- veg_dir %>% list.files(pattern = "tif$"
#                                     , full.names = TRUE
#                                     , recursive = TRUE
# ) %>% .[!grepl('count', .)]
# for(i in 1:length(veg_files))
# {
#   f_i <- veg_files[[i]]
#   outname <- f_i %>% basename %>% paste0(cov_dir, '/', .)
#   r <- f_i %>%
#     rast %>%
#     classify(rcl = reclasser
#              , filename = outname
#     )
# }


# 3: List all files and their categories

cov_dir %>% list.files %>% write.table(
  file = paste0(root, '/cov_list.csv')
  , sep = ';'
  )

cov_cats <- root %>%
  paste0(., '/cov_categories.csv') %>%
  read.table(sep = ';'
             , header = TRUE
             )

# Anything missing?

cov_names <- cov_dir %>% list.files

cov_names[!cov_names %in% cov_cats$name]

# 4: Load DSC points

dsc <- root %>%
  paste0(
    .
    , '/observations/DanishSoilClassification/DLJ/DLJmDecimaler_DKBH.shp'
  ) %>%
  vect

# Fix comma decimals (not relevant for new file)

# fix_these <- c('LABORATORI', 'LER', 'SILT', 'GROVSILT', 'GROVFINSAN', 'FINSAND'
#                ,'GROVSAND', 'CACO3', 'MULDTLSTAN')
# 
# values(dsc) %<>%
#   mutate(across(all_of(fix_these)
#                        , str_replace
#                        , ','
#                        , '.'
#                        )) %>%
#   mutate(across(all_of(fix_these)
#                                             , as.numeric
#                        ))

# 5: check if the clay values give percent of the mineral fraction

mineral <- c('Ler', 'Silt', 'FinSD', 'GrovSD')

dsc %>%
  values %>%
  select(all_of(mineral)) %>%
  apply(., 1, FUN = function(x) sum(x, na.rm = TRUE))

# They do not

# 6: Test which values add up to 100

# test_these <- c('Ler', 'Silt', 'FinSD', 'GrovSD', 'Humus', 'Gsilt', 'GfinSD'
#                 , 'CaCO3')
# 
# test_list <- lapply(1:length(test_these)
#                     , FUN = function(x)
#                     {
#                       combn(test_these, x, simplify = FALSE)
#                     }
# ) %>%
#   unlist(recursive = FALSE)
# 
# tex_sums <- numeric()
# merged_names <- lapply(test_list, function(x){paste0(x, collapse = '+')}) %>%
#   unlist()
# 
# for(i in 1:length(test_list))
# {
#   tex_sums[i] <- dsc %>%
#     values %>%
#     select(all_of(test_list[[i]])) %>%
#     apply(., 1, FUN = function(x) sum(x, na.rm = TRUE)) %>%
#     median(na.rm = TRUE)
#   c(i, tex_sums[i], merged_names[i]) %>% print
# }
# 
# merged_names[tex_sums == 100]

# These add up to 100
# "Ler+Silt+FinSD+GrovSD+Humus"
# "Ler+Silt+FinSD+GrovSD+Humus+CaCO3"       
# "Ler+Silt+GrovSD+Humus+Gsilt+GfinSD"
# "Ler+Silt+GrovSD+Humus+Gsilt+GfinSD+CaCO3"

# These ones come close

# tibble(sums = tex_sums
#        , vars = merged_names
#        , row = 1:length(test_list)) %>%
#   mutate(diff = abs(sums - 100)) %>%
#   arrange(diff)

#       sums vars                                       row  diff
#      <dbl> <chr>                                    <int> <dbl>
#   1  100   Ler+Silt+FinSD+GrovSD+Humus                163 0    
#   2  100   Ler+Silt+FinSD+GrovSD+Humus+CaCO3          221 0    
#   3  100   Ler+Silt+GrovSD+Humus+Gsilt+GfinSD         229 0    
#   4  100   Ler+Silt+GrovSD+Humus+Gsilt+GfinSD+CaCO3   252 0    
#   5   99.7 Silt+FinSD+GrovSD+Humus+Gsilt+CaCO3        241 0.270
#   6   99.6 Silt+FinSD+GrovSD+Humus+Gsilt              198 0.370
#   7   98.6 Ler+FinSD+GrovSD+Humus+Gsilt+CaCO3         235 1.35 
#   8   98.6 Ler+FinSD+GrovSD+Humus+Gsilt               183 1.43 
#   9   97.2 Ler+Silt+FinSD+GrovSD+CaCO3                166 2.83 
#   10  97.2 Ler+Silt+FinSD+GrovSD                       93 2.84   

#   2  100   Ler+Silt+FinSD+GrovSD+Humus+CaCO3          221 0
# dsc %>%
#   values %>%
#   select(all_of(test_list[[221]])) %>%
#   apply(., 1, FUN = function(x) sum(x, na.rm = TRUE))
# 
# dsc %>%
#   values %>%
#   mutate(sum_ex = Ler+Silt+FinSD+GrovSD+Humus+CaCO3) %>%
#   filter(CaCO3 > 0)

# CACO3 is probably not part of the texture sum, but Humus is.

# Correct clay contents so they represent the percentage of the mineral fraction

mineral_cols <- c('Ler', 'Silt', 'FinSD', 'GrovSD')

# dsc$mineral <- values(dsc) %>%
#   mutate(mineral = across(all_of(mineral_cols), sum)) %>%
#   
dsc$clay <- values(dsc) %>% 
  select(all_of(mineral_cols)) %>%
  mutate(clay = Ler * 100 / rowSums(., na.rm = TRUE), .keep = 'none') %>%
  .$clay

# log transform HUMUS

dsc$ln_SOM <- log(dsc$Humus)

# 7: Extracting covariates

# List all covariates

cov_files <- cov_dir %>% list.files(pattern = "tif$"
                                    , full.names = TRUE)

cov <- cov_files %>% rast

# The code below fixed the extents and projections of the covariates

# # Get the extents of all covariates
# 
# extents <- lapply(cov_files, function(x)
# {
#   out <- x %>% rast %>% ext %>% .@ptr %>% .$vector
#   return(out)
# }
# ) %>%
#   bind_cols %>%
#   as.matrix %>%
#   t

base_covnames <- basename(cov_files)

# cbind(extents, base_covnames)

# Crop all covariates to the same extent as the bare soil composite

cov_crop_dir <- root %>% paste0(., '/covariates_cropped/') %T>% dir.create()

terraOptions(tempdir = paste0(root, '/Temp'))

bare_ext <- root %>%
  paste0(., '/layers/bare-soil/s2-geomedian-B2.tif') %>%
  rast %>%
  ext

# for(i in 2:length(cov_files))
# {
#   r_i <- cov_files[[i]] %>% rast
#   
#   if(ext(r_i) != bare_ext)
#   {
#     if(nlyr(r_i) > 1)
#     {
#       fname <- paste0(cov_crop_dir, '/', names(r_i), '.tif')
#     } else {
#       fname <- paste0(cov_crop_dir, '/', base_covnames[i])
#     }
#     r_i %>% math('signif'
#                  , digits = 4
#                  , filename = paste0(root, '/Temp/tmp1.tif')
#                  , overwrite = TRUE) %>%
#       crop(bare_ext
#            , filename = paste0(root, '/Temp/tmp2.tif')
#            , overwrite = TRUE
#       ) %>%
#       writeRaster(filename = fname
#                   , overwrite = TRUE
#       )
#   }
#   print(cov_files[i])
# }

# Not working for categorical variables
# - geology (wrong numbers format)
# - georeg (missing projection, wrong numbers format)
# - IMK layers (wrong numbers format)
# - landscape (missing projection, wrong numbers format)
# - wetlands_10m (wrong numbers format)
# - LU

# Missing projection
# - S1_2020_1_VH.tif (etc)

# These look wierd:
# - General_Curvature: Nearly no variation
# - Flow_Line_Curvature: Nearly no variation
# - MRRTF: Rectangular artefacts
# - MRVBF: Rectangular artefacts
# - Negative_Openness: Angular artefacts
# - Flow_Accumulation_ArcGISPro: Dotted patterns in filled sinks
# - Flow_Accumulation_top-down: Angular flow lines in filled sinks
# - Ln_Flow_Accumulation_ArcGISPro_Planchon_Darboux: Dotted patterns in filled sinks

# Not useful:
# - Flow_Direction_ArcGISPro_Planchon_Darboux.tif
# - Flow_Direction_ArcGISPro.tif
# - LnGte9TimesDEM0_Flow_Accumulation_ArcGISPro_Planchon_Darboux.tif
# - LnGte8TimesDEM1_Flow_Accumulation_ArcGISPro.tif
# - LnGte9TimesDEM1Nibble_Flow_Accumulation_ArcGISPro_Planchon_Darboux.tif
# - LnGte8TimesDEM1Nibble_Flow_Accumulation_ArcGISPro.tif
# - LnGte9TimesDEM1NibbleGte0_Flow_Accumulation_ArcGISPro_Planchon_Darboux.tif
# - LnGte8TimesDEM1NibbleGte0_Flow_Accumulation_ArcGISPro.tif

# Fix the remaining variables

cov_files <- cov_dir %>%
  list.files(pattern = "tif$"
             , full.names = TRUE)

proj_ETRS89 <- root %>%
  paste0(., '/layers/bare-soil/s2-geomedian-B2.tif') %>%
  rast %>%
  crs

# for(i in 1:length(cov_files))
# {
#   r_i <- cov_files[[i]] %>% rast
# 
#   if(ext(r_i) != bare_ext)
#   {
#     if(nlyr(r_i) > 1)
#     {
#       fname <- paste0(cov_crop_dir, '/', names(r_i), '.tif')
#     } else {
#       fname <- paste0(cov_crop_dir, '/', base_covnames[i])
#     }
#     dtyp <- raster::dataType(raster::raster(cov_files[[i]]))
#     r_out <- r_i
#     crs(r_out) <- proj_ETRS89
#     r_out %>%
#       # math('signif'
#       #            , digits = 4
#       #            , filename = paste0(root, '/Temp/tmp1.tif')
#       #            , overwrite = TRUE) %>%
#       crop(bare_ext
#            , filename = paste0(root, '/Temp/tmp2.tif')
#            , overwrite = TRUE
#       ) %>%
#       writeRaster(filename = fname
#                   , overwrite = TRUE
#                   , datatype = dtyp
#       )
#     print(cov_files[i])
#   }
# }

cov_files <- cov_dir %>% list.files(pattern = "tif$"
                                    , full.names = TRUE)

cov_names_sans <- cov_files %>%
  basename %>%
  tools::file_path_sans_ext(.) %>%
  make.names()

cov <- cov_files %>% rast

names(cov) <- cov_names_sans

# dsc_extr <- terra::extract(
#   cov
#   , dsc
# )

extr_file <- root %>% paste0(., '/extracts/dsc_extr_all.csv')

# Drop these (Bornholm missing)
# Area_Solar_RadiationSpecialDays4hour_ArcGISPro_1
# Area_Solar_RadiationSpecialDays4hour_ArcGISPro_2
# Area_Solar_RadiationSpecialDays4hour_ArcGISPro_3

# dsc_extr %<>% select(-c(2:4))

# dsc_extr %>%
#   write.table(
#     file = extr_file
#     , sep = ';'
#     , row.names = FALSE
# )

dsc_extr <- extr_file %>%
  read.table(
    header = TRUE
    , sep = ';'
  )
  
# 8: Create folds

library(caret)

# set.seed(1)
# 
# dsc %>%
#   nrow %>%
#   sample(1:10, ., replace = TRUE) %>%
#   data.frame(fold = .) %>%
#   write.table(
#     file = 'dsc_folds_all.csv'
#     , sep = ';'
#     , row.names = FALSE
# )

folds_10 <- read.table(
  file = 'dsc_folds_all.csv'
  , header = TRUE
  , sep = ';'
)

# 10: Compile data

dsc_cov_folds <- bind_cols(
  values(dsc)
  , select(dsc_extr, -1)
  , folds_10
)

tr_dat <- dsc_cov_folds %>%
  filter(
    DybFra == 0
    , DybTil == 20
    , s2.geomedian.B2 > 0
  ) %>%
  filter_at(vars(all_of(names(cov)))
            , all_vars(!is.na(.)))

tr_dat_clay <- tr_dat %>%
  filter(!is.na(clay))

tr_dat_SOM <- tr_dat %>%
  filter(!is.infinite(ln_SOM))

# 11: Set up models
# Making sure the names match and are legitimate

f_clay <- names(cov) %>%
  paste(collapse = ' + ') %>%
  paste0('clay ~ ', .) %>%
  as.formula

f_SOM <- names(cov) %>%
  paste(collapse = ' + ') %>%
  paste0('ln_SOM ~ ', .) %>%
  as.formula

folds_clay <- sapply(
  1:10,
  function(x)
  {
    ind <- tr_dat_clay$fold != x
    out <- c(1:length(ind))[ind]
    return(out)
  }
)

folds_SOM <- sapply(
  1:10,
  function(x)
  {
    ind <- tr_dat_SOM$fold != x
    out <- c(1:length(ind))[ind]
    return(out)
  }
)

tgrid <- data.frame(
  committees = c(1, 10, 20)
  , neighbors = c(0, 0, 0))

# Clay model

library(doParallel)

showConnections()

cl <- makePSOCKcluster(12)
registerDoParallel(cl)

model_clay <- caret::train(
  f_clay
  , data = tr_dat_clay
  , method = 'cubist'
  # , trControl = trainControl(method = 'cv'
  #                            , number = 12)
  , tuneGrid = tgrid
  , trControl = trainControl(
    index = folds_clay
    , savePredictions = 'final'
    , predictionBounds = c(0, 100)
    )
)

registerDoSEQ()
rm(cl)

model_clay
model_clay %>% varImp
model_clay %>%
  varImp %>%
  .$importance %>%
  rownames_to_column(var = 'covariate') %>%
  arrange(-Overall)

# SOM model

set.seed(1)

showConnections()

cl <- makePSOCKcluster(12)
registerDoParallel(cl)

model_SOM <- caret::train(
  f_SOM
  , data = tr_dat_SOM
  , method = 'cubist'
  , tuneGrid = tgrid
  , trControl = trainControl(
    index = folds_SOM
    , savePredictions = 'final'
    , predictionBounds = c(NA, log(100))
    )
)

registerDoSEQ()
rm(cl)

showConnections()

model_SOM
model_SOM %>% varImp
model_SOM %>%
  varImp %>%
  .$importance %>%
  rownames_to_column(var = 'covariate') %>%
  arrange(-Overall)

# Make covariates for small test area

squareshape <- root %>%
  paste0(., '/testarea_10km/square10km.shp') %>%
  vect

square_ext <- squareshape %>%
  ext %>%
  round(-1)

outfolder <- root %>%
  paste0(., '/testarea_10km/covariates/')

# outfolder %>% dir.create
# 
# cropstack <- function(x # list of files
#                       , y # extent
#                       , folder # target folder
#                       )
# {
#   for(i in 1:length(x))
#   {
#     r <- x[i] %>% rast
#     dtype <- r %>%
#       sources %>%
#       raster::raster(.) %>%
#       raster::dataType(.)
#     outname <- r %>%
#       sources %>%
#       basename %>%
#       tools::file_path_sans_ext(.) %>%
#       make.names %>%
#       paste0(folder, ., '.tif')
#     r %>%
#       crop(y = y) %>%
#       writeRaster(datatype = dtype
#                   , filename = outname
#                   , overwrite = TRUE
#       )
#   }
# }

# cov_files %>%
#   cropstack(y = square_ext
#             , folder = outfolder
#             )

# Test for small area

cov_10km <- outfolder %>%
  list.files(full.names = TRUE) %>%
  rast

names(cov_10km) <- names(cov)

predfolder <- root %>%
  paste0(., '/testarea_10km/predictions_', testn, '/') %T>%
  dir.create()

t1 <- Sys.time()
 
predict(cov_10km
        , model_clay
        , na.rm = TRUE
        , cores = 12
        , filename = paste0(predfolder, 'clay.tif')
        , overwrite = TRUE
        )

Sys.time() - t1
# Time difference of 2.222116 mins

t1 <- Sys.time()

predict(cov_10km
        , model_SOM
        , na.rm = TRUE
        , cores = 12
        , filename = paste0(predfolder, 'SOM.tif')
        , overwrite = TRUE
)

Sys.time() - t1

# How come the triangular artefacts?
# OGC, rounding

# Import Vindum data for comparison

# exp transformation for SOM outputs

# Move everything to ssd?

# Make maps for all of Denmark

# END