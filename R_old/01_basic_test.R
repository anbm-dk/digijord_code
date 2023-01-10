# Script 1: First, basic test

library(raster)
library(magrittr)
library(dplyr)
library(caret)
library(readr)
library(tidyr)

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

source('loadandstack.R')

pts <- paste0(root, '/DanishSoilClassification/DSC_ETRS89/dsc_pts.shp') %>%
  shapefile

locale_comma <- readr::locale(decimal_mark = ",")

pts$LER %<>% readr::parse_number(locale = locale_comma)
pts$SILT %<>% readr::parse_number(locale = locale_comma)
pts$GROVSILT %<>% readr::parse_number(locale = locale_comma)
pts$GROVFINSAN %<>% readr::parse_number(locale = locale_comma)
pts$FINSAND %<>% readr::parse_number(locale = locale_comma)
pts$GROVSAND %<>% readr::parse_number(locale = locale_comma)
pts$HUMUS %<>% readr::parse_number(locale = locale_comma)
pts$CACO3 %<>% readr::parse_number(locale = locale_comma)

pts_top <- pts[pts$DYBDE_FRA == 0, ]

s <- paste0(root, '/cov_basic') %>% loadandstack()

# extr <- raster::extract(s, pts_top)
# 
# extr %>% write.table(file = 'dsc_top_extr_basic.csv'
#                      , sep = ';'
#                      )

extr <- read.table(file = 'dsc_top_extr_basic.csv'
                                        , sep = ';'
                                        )

extr$geology_10m %<>% as.factor
extr$georeg_10m %<>% as.factor
extr$landscapes_10m %<>% as.factor

pts_extr <- pts_top

pts_extr@data <- cbind(pts_top@data, extr)

texclas <- c('LER', 'SILT', 'GROVSILT', 'GROVFINSAN', 'FINSAND', 'GROVSAND'
             , 'HUMUS', 'CACO3')

# names(pts_extr@data)



# BARE SOIL COMPOSITE

# Load data

bare <- paste0(root, '/digijord_download/bare-soil') %>% loadandstack()

# extr_bare <- raster::extract(bare, pts_top)
# 
# extr_bare %>% write.table(file = 'dsc_top_extr_bare.csv'
#                           , sep = ';'
# )

extr_bare <- read.table(file = 'dsc_top_extr_bare.csv'
                        , sep = ';'
)

pts_bare <- pts_top

pts_bare@data <- cbind(pts_top@data, extr_bare)


# Basic models

models_basic <- list()

for(i in 1:length(texclas))
# for(i in c(1, 7))
{
  dat <- cbind(pts_top@data[, names(pts_top) == texclas[i]]
               , extr) #!
  
  dat %<>% as.data.frame
  colnames(dat)[1] <- texclas[i]
  dat %<>% drop_na
  
  fm_i <- names(s) %>% #!
    paste(collapse = ' + ') %>%
    paste0(texclas[i], ' ~ ', .) %>%
    as.formula
  
  models_basic[[i]] <- caret::train(fm_i
                              , data = dat
                              , method = 'cubist'
                              , trControl = trainControl(method = 'cv'
                                                         , number = 3)
  )
}


# Bare soil models

models_bare <- list()

for(i in 1:length(texclas))
# for(i in c(1, 7))
{
  dat <- cbind(pts_top@data[, names(pts_top) == texclas[i]]
               , extr_bare) #!
  
  dat %<>% as.data.frame
  colnames(dat)[1] <- texclas[i]
  dat %<>% drop_na
  
  fm_i <- names(bare) %>% #!
    paste(collapse = ' + ') %>%
    paste0(texclas[i], ' ~ ', .) %>%
    as.formula
  
  models_bare[[i]] <- caret::train(fm_i
                                    , data = dat
                                    , method = 'cubist'
                                    , trControl = trainControl(method = 'cv'
                                                               , number = 3)
  )
}

# Combined models

pts_comb <- pts_top

pts_comb@data <- cbind(pts_top@data, extr, extr_bare)

models_comb <- list()

for(i in 1:length(texclas))
# for(i in c(1, 7))
{
  dat <- cbind(pts_top@data[, names(pts_top) == texclas[i]]
               , extr
               , extr_bare
               ) #!
  
  dat %<>% as.data.frame
  colnames(dat)[1] <- texclas[i]
  dat %<>% drop_na
  
  
  fm_i <- c(names(s), names(bare)) %>% #!
    paste(collapse = ' + ') %>%
    paste0(texclas[i], ' ~ ', .) %>%
    as.formula
  
  models_comb[[i]] <- caret::train(fm_i
                                   , data = dat
                                   , method = 'cubist'
                                   , trControl = trainControl(method = 'cv'
                                                              , number = 3)
  )
}

# END