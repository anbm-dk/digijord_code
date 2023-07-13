# Tests for 10 km square area

# A: Get ready

library(raster)
library(tidyr)
library(magrittr)
library(rgdal)
library(dplyr)
library(caret)
library(readr)
library(lubridate)
library(foreach)
library(tibble)
library(stringr)
library(doParallel)
library(ggupset) # for combination plots
library(forcats) # reordering factors

# Test number

testn <- 2

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

source("loadandstack.R")

squareshape <- root %>%
  paste0(., "/testarea_10km/square10km.shp") %>%
  shapefile()

square_ext <- squareshape %>%
  extent() %>%
  round(-1)

outfolder <- root %>%
  paste0(., "/testarea_10km/covariates/")

outfolder %>% dir.create()

cropstack <- function(x, y, folder, replaceNA = FALSE) {
  dtypes <- x %>%
    dataType()

  rnames <- x %>%
    names()

  x_cropped <- x %>%
    crop(y)

  if (replaceNA == TRUE) {
    x_cropped[is.na(x_cropped[])] <- 0
  }

  for (i in 1:nlayers(x_cropped))
  {
    outname <- folder %>%
      paste0(., rnames[i], ".tif")

    raster::writeRaster(x_cropped[[i]],
      datatype = dtypes[i],
      filename = outname,
      overwrite = TRUE
    )
  }

  return(x_cropped)
}

# B: Basic covariates

# cov_basic <- root %>%
#   paste0(., '/cov_basic/') %>%
#   loadandstack()
#
# cov_basic %>%
#   cropstack(y = square_ext
#             , folder = outfolder
#             )

# C: Bare soil composite

# cov_bare <- root %>%
#   paste0(., '/digijord_download/bare-soil/') %>%
#   loadandstack()
#
# cov_bare %>%
#   cropstack(y = square_ext
#             , folder = outfolder
#   )

# D: Radar

# cov_radar <- root %>%
#   paste0(., '/Radar/processed/') %>%
#   loadandstack()
#
# cov_radar %>%
#   cropstack(y = square_ext
#             , folder = outfolder
#   )

# E: Time series

# rdirs <- paste0(root, '/digijord_download') %>%
#   list.dirs %>%
#   grep("Bornholm", .
#        , invert = TRUE
#        , value = TRUE) %>%
#   grep("bare", .
#        , invert = TRUE
#        , value = TRUE)
#
# for(j in 1:length(rdirs))
# {
#   dir_j <- rdirs[j]
#
#   ntifs <- dir_j %>%
#     list.files %>%
#     grep('.tif', .) %>%
#     length
#
#   if(ntifs > 0)
#   {
#     s <- loadandstack(dir_j)
#
#     s %>%
#       cropstack(y = square_ext
#                 , folder = outfolder
#                 , replaceNA = TRUE
#       )
#   }
#
# }

rm(j, i, ntifs, dir_j, rdirs, outname, s, s2, square_r, folders, rnames)

# F: Redo big test, but also make some maps

# 1 Skipped here

# 2 Load DSC points

pts <- paste0(root, "/DanishSoilClassification/DSC_ETRS89/dsc_pts.shp") %>%
  shapefile()

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

texclas <- c(
  "LER", "SILT", "GROVSILT", "GROVFINSAN", "FINSAND", "GROVSAND",
  "HUMUS", "CACO3"
)

# We skip calcium carbonates, as we couldn't predict them last time
texclas %<>% setdiff("CACO3")

pts_top_tex <- pts_top@data[, names(pts_top) %in% texclas]

texclas <- names(pts_top_tex)

# 3 Load covariates

## 3.1 Basic covariates (redo extraction if necessary)

# cov_basic <- paste0(root, '/cov_basic') %>% loadandstack()
#
# extr_basic <- raster::extract(cov_basic, pts_top)
#
# extr_basic <- extr_basic[, -c(1:2)]
#
# extr_basic %>% write.table(file = 'dsc_top_extr_basic.csv'
#                      , sep = ';'
#                      , row.names = FALSE
#                      )

extr_basic <- read.table(
  file = "dsc_top_extr_basic.csv",
  sep = ";",
  header = TRUE
)

nacolumns <- lapply(extr_basic, function(x) sum(is.na(x))) %>% unlist()

extr_basic <- extr_basic[, nacolumns < 100]

# extr_basic$geology_10m %<>% as.factor
# extr_basic$georeg_10m %<>% as.factor
# extr_basic$landscapes_10m %<>% as.factor

## 3.2 Bare soil composite

extr_bare <- read.table(
  file = "dsc_top_extr_bare.csv",
  sep = ";"
)

## 3.3 Vegetated time series

# extrs_merged <- paste0(root, '/sat_extracts_merged')
#
# extrs_merged_files <- extrs_merged %>% list.files(full.names = TRUE)
#
# f_nas <- function(x)
# {
#   d <- read.table(x
#                   , sep = ';'
#                   , header = TRUE)
#   out <- is.na(d[, 1]) %>% mean
#   return(out)
# }
#
# n_nas <- extrs_merged_files %>% sapply(f_nas) %>% unname()
#
# extrs_formodel <- extrs_merged_files[n_nas < 0.5]
#
# f_read2 <- function(x)
# {
#   out <- read.table(x
#                     , sep = ';'
#                     , header = TRUE)
#   return(out)
# }
#
# extrs_formodel %<>%
#   lapply(f_read2) %>%
#   bind_cols()
#
# extrs_formodel[is.na(extrs_formodel)] <- 0
#
# extrs_formodel %>%
#   write.table(file = 'dsc_top_extr_veg.csv'
#               , sep = ';'
#               , row.names = FALSE
#   )
#
# rm(extrs_formodel, f_read2, n_nas, f_nas, extrs_merged_files, extrs_merged)

extr_veg <- read.table(
  file = "dsc_top_extr_veg.csv",
  sep = ";",
  header = TRUE
)

## 3.4 Radar images

extr_radar <- read.table(
  file = "dsc_top_extr_radar.csv",
  sep = ";",
  header = TRUE
)

# 4 Combined models

# Prepare for the test

getimps <- function(x) {
  d <- x$importance
  d$image <- rownames(d)
  d %<>% arrange(image)
  return(d)
}

getimps_list <- function(x, ynames) {
  getimps <- function(x) {
    d <- x$importance
    d$name <- rownames(d)
    d %<>% arrange(name)
    return(d)
  }

  x %<>% lapply(getimps)
  varnames <- x[[1]]$name
  scores <- x %>%
    lapply("[", "Overall") %>%
    bind_cols(varnames, .) %>%
    set_names(., c("Variable", ynames))

  return(scores)
}

# Tuning grid for simple model

tgrid_simp <- data.frame(
  committees = 1,
  neighbors = 0
)

# Tuning grid for optimized model

tgrid_opt <- data.frame(
  committees = c(1, 10, 20),
  neighbors = c(0, 0, 0)
)

# Grid for covariate types

cov_grid <- expand_grid(
  basic = 0:1,
  bare = 0:1,
  series = 0:1,
  radar = 0:1
) %>%
  filter(!row_number() %in% 1)

# pts_top_tex$CACO3 %<>% replace_na(0)  # No longer necessary

all_imps <- list()

all_imps_opt <- list()

all_RMSEs <- matrix(numeric(),
  ncol = ncol(pts_top_tex),
  nrow = nrow(cov_grid),
  dimnames = list(
    1:nrow(cov_grid),
    colnames(pts_top_tex)
  )
)

all_R2s <- all_RMSEs
all_times <- all_RMSEs


# createDataPartition(1:10, times = 10, p = 0.5) # Create data partition (later)

trc <- trainControl(
  method = "LGOCV",
  number = 10,
  p = 0.5
)

# Load covariates for test area

covs_test <- root %>%
  paste0(., "/testarea_10km/covariates") %>%
  loadandstack()

# Make folder for predictions

pred_dir <- root %>%
  paste0(., "/testarea_10km/predictions_", testn) %T>%
  dir.create()

# Names of the prediction covariates

pred_covs <- cov_grid %>%
  apply(1, function(x) {
    out <- x[x == 1] %>%
      names() %>%
      paste(collapse = "_")
    return(out)
  })

# List for prediction models

modelsformaps <- list()

# Files for accuracies

rmse_file <- paste0("bigtest", testn, "_rmses.csv")
r2_file <- paste0("bigtest", testn, "_r2s.csv")

# Run test

set.seed(4914)

for (i in 15)
# for(i in 1:nrow(cov_grid))  # i is the row in the covariate grid
{
  all_imps[[i]] <- list()

  covs <- list()

  j <- 1

  if (cov_grid$basic[i] == 1) {
    covs[[j]] <- extr_basic
    j %<>% +(1)
  }
  if (cov_grid$bare[i] == 1) {
    covs[[j]] <- extr_bare
    j %<>% +(1)
  }
  if (cov_grid$series[i] == 1) {
    covs[[j]] <- extr_veg
    j %<>% +(1)
  }
  if (cov_grid$radar[i] == 1) {
    covs[[j]] <- extr_radar
    j %<>% +(1) # Not really necessary, but I'll keep it for now
  }

  covs %<>% bind_cols()

  for (k in 1:ncol(pts_top_tex)) # k is is texture class in the vector
  {
    dat <- pts_top_tex %>%
      select(all_of(k)) %>%
      bind_cols(., covs) %>%
      drop_na()

    fm_k <- covs %>%
      colnames() %>%
      paste(collapse = " + ") %>%
      paste0(colnames(dat)[1], " ~ ", .) %>%
      as.formula()

    cl <- makePSOCKcluster(10)
    registerDoParallel(cl)

    model_k <- caret::train(fm_k,
      data = dat,
      method = "cubist",
      tuneGrid = tgrid_simp,
      trControl = trc
    )

    registerDoSEQ()
    rm(cl)

    all_imps[[i]][[k]] <- model_k %>% varImp()

    all_RMSEs[i, k] <- model_k$results$RMSE
    all_R2s[i, k] <- model_k$results$Rsquared



    all_RMSEs %>%
      write.table(
        file = rmse_file,
        sep = ";",
        row.names = FALSE
      )

    all_R2s %>%
      write.table(
        file = r2_file,
        sep = ";",
        row.names = FALSE
      )

    # Write maps for the last set of covariates

    if (i == nrow(cov_grid)) {
      top100 <- model_k %>%
        varImp() %>%
        .$importance %>%
        as_tibble(rownames = "covariate") %>%
        arrange(-Overall) %>%
        slice_head(n = 100)

      predvars <- top100$covariate

      covs_test_k <- covs_test %>%
        subset(subset = predvars)

      pred_name <- pred_dir %>%
        paste0(
          ., "/P", k, "_", texclas[k],
          "_M", formatC(i,
            width = cov_grid %>%
              nrow() %>%
              nchar(),
            flag = "0"
          ),
          "_", pred_covs[i],
          ".tif"
        )

      f_top100 <- covs_test_k %>%
        names() %>%
        paste(collapse = " + ") %>%
        paste0(colnames(dat)[1], " ~ ", .) %>%
        as.formula()

      cl <- makePSOCKcluster(10)
      registerDoParallel(cl)

      modelsformaps[[k]] <- caret::train(
        f_top100,
        data = dat,
        method = "cubist",
        tuneGrid = tgrid_opt,
        trControl = trc
      )

      registerDoSEQ()
      rm(cl)

      beginCluster(12)

      predict(
        covs_test_k,
        modelsformaps[[k]],
        filename = pred_name
      )

      endCluster()
    }
  }

  outname <- paste0(
    "bigtest", testn, "_imps_",
    formatC(i,
      width = cov_grid %>%
        nrow() %>%
        nchar(),
      flag = "0"
    ),
    ".csv"
  )

  all_imps[[i]] %>%
    getimps_list(ynames = colnames(pts_top_tex)) %>%
    write.table(
      file = outname,
      sep = ";",
      row.names = FALSE
    )
}

# # Test rfe (doesn't seem to work very well)
#
# rfe_control <- rfeControl(functions = caretFuncs   #caretFuncs here
#                           , method = "cv"
#                           , number = 12)
#
# siz_k <- covs %>%
#   ncol %>%
#   `/` (10) %>%
#   floor %>%
#   `*` (10) %>%
#   seq(10, ., 10)
#
# subset_k <- sample(1:nrow(dat), 1000)
#
# cl <- makePSOCKcluster(12)
# registerDoParallel(cl)
#
# rfe_fit <- rfe(fm_k
#                , data = dat
#                , sizes = siz_k
#                , rfeControl = rfe_control
#                , method = 'cubist'
#                , tuneGrid = tr_grid
#                , subset = subset_k
#                , metric = "Rsquared"
#                )
#
# registerDoSEQ()
# rm(cl)
#
# rfe_fit$optVariables
#
# plot(rfe_fit$results$RMSE, rfe_fit$results$Rsquared)
#
# plot(rfe_fit$results$Rsquared)
#
# f_rfe <- rfe_fit$optVariables %>%
#   paste(collapse = ' + ') %>%
#   paste0(colnames(dat)[1], ' ~ ', .) %>%
#   as.formula
#
# f_rfe <- covs_rfe %>%
#   names %>%
#   paste(collapse = ' + ') %>%
#   paste0(colnames(dat)[1], ' ~ ', .) %>%
#   as.formula
#
# cl <- makePSOCKcluster(12)
# registerDoParallel(cl)
#
# model_rfe <- caret::train(
#   # fm_k  # standard formula
#   f_rfe
#   , data = dat
#   , method = 'cubist'
#   # , tuneGrid = tr_grid
#   , trControl = trainControl(method = 'cv'
#                              , number = 12)
# )
#
# registerDoSEQ()
# rm(cl)
#
# covs_rfe <- covs_test %>%
#   subset(subset = rfe_fit$optVariables)

# 5 Present results

plots <- root %>%
  paste0(., "bigtest", testn, "_plots/")

plots %>% dir.create()

cov_tib <- (cov_grid == 1) %>%
  as_tibble(rownames = "row") %>%
  pivot_longer(-row) %>%
  filter(value) %>%
  select(-value) %>%
  group_by(row) %>%
  summarise(Covariates = list(name)) %>%
  mutate(row = as.numeric(row)) %>%
  arrange(row)

# Plot for R2

for (i in 1:ncol(all_R2s))
{
  cov_tib_i <- cov_tib %>%
    add_column(R2 = all_R2s[, i]) %>%
    arrange(R2) %>%
    mutate(cov_collapsed = sapply(
      Covariates,
      function(x) {
        paste0(sort(x),
          collapse = "-"
        )
      }
    )) %>%
    mutate(cov_collapsed = fct_reorder(cov_collapsed, R2, .desc = TRUE))

  tiff(
    filename = paste0(plots, "R2_", i, "_", colnames(all_R2s)[i], ".tif"),
    units = "cm",
    width = 16 * 1.5,
    height = 10 * 1.5,
    res = 600
  )

  p <- cov_tib_i %>%
    ggplot(
      aes(
        x = cov_collapsed,
        y = R2
      )
    ) +
    geom_col() +
    axis_combmatrix(sep = "-") +
    ggtitle(colnames(all_R2s)[i]) +
    xlab("Covariates")

  print(p)

  dev.off()
  # dev.off()
}


# Plot for RMSE

for (i in 1:ncol(all_RMSEs))
{
  cov_tib_i <- cov_tib %>%
    add_column(RMSE = all_RMSEs[, i]) %>%
    arrange(RMSE) %>%
    mutate(cov_collapsed = sapply(
      Covariates,
      function(x) {
        paste0(sort(x),
          collapse = "-"
        )
      }
    )) %>%
    mutate(cov_collapsed = fct_reorder(cov_collapsed, RMSE, .desc = FALSE))

  tiff(
    filename = paste0(plots, "RMSE_", i, "_", colnames(all_RMSEs)[i], ".tif"),
    units = "cm",
    width = 16 * 1.5,
    height = 10 * 1.5,
    res = 600
  )

  p <- cov_tib_i %>%
    ggplot(
      aes(
        x = cov_collapsed,
        y = RMSE
      )
    ) +
    geom_col() +
    axis_combmatrix(sep = "-") +
    ggtitle(colnames(all_RMSEs)[i]) +
    xlab("Covariates")

  print(p)

  dev.off()
  # dev.off()
}

# Covariate importance

best_R2s_ind <- all_R2s %>% apply(2, which.max)
all_RMSEs %>% apply(2, which.min)

all_R2s %>%
  apply(2, function(x) max(x, na.rm = TRUE)) %>%
  plot()
all_RMSEs %>%
  apply(2, function(x) min(x, na.rm = TRUE)) %>%
  plot()

plot(all_imps[[15]][[1]], top = 20)

all_imps[[15]][[1]]$importance %>%
  as_tibble(rownames = "covariate") %>%
  arrange(-Overall) %>%
  slice_head(n = 20)

l <- list()

ntop <- 20

for (i in 1:length(best_R2s_ind))
{
  best_R2 <- which.max(all_R2s[, i])
  l[[i]] <- all_imps[[best_R2]][[i]]$importance %>%
    as_tibble(rownames = "covariate") %>%
    drop_na() %>%
    arrange(-Overall) %>%
    slice_head(n = ntop) %>%
    mutate(target = colnames(all_RMSEs)[i]) %>%
    mutate(rank = 1:ntop)
}

l %<>% bind_rows() %>%
  mutate(target = factor(target, levels = texclas))

basic_ind <- which(cov_grid$basic == 1)[1]

l_cat <- list()

i <- 1
l_cat[[i]] <- extr_bare %>%
  colnames() %>%
  tibble(
    covariate = .,
    category = "bare"
  )
i %<>% +(1)
l_cat[[i]] <- all_imps[[basic_ind]][[1]]$importance %>% # ugly solution due to factors
  rownames() %>%
  tibble(
    covariate = .,
    category = "basic"
  )
i %<>% +(1)
l_cat[[i]] <- extr_radar %>%
  colnames() %>%
  tibble(
    covariate = .,
    category = "radar"
  )
i %<>% +(1)
l_cat[[i]] <- extr_veg %>%
  colnames() %>%
  tibble(
    covariate = .,
    category = "series"
  )

l_cat %<>% bind_rows()

l %<>%
  left_join(l_cat)

l %<>%
  ungroup() %>%
  arrange(target, Overall) %>%
  mutate(order = row_number())

l %>%
  ggplot(aes(x = order, y = Overall, bg = category)) +
  geom_col() +
  facet_wrap(~target,
    ncol = 3,
    scales = "free"
  ) +
  xlim(1, ntop) +
  ylim(0, 100) +
  coord_flip() +
  scale_x_continuous(
    breaks = l$order,
    labels = l$covariate,
    expand = c(0, 0)
  )


pdf(
  file = "bigtest", testn, "_importance.pdf",
  paper = "a4",
  width = 0,
  height = 0
)

l %>%
  ggplot(aes(x = order, y = Overall, bg = category)) +
  geom_col(width = 1) +
  facet_wrap(~target,
    ncol = 2,
    scales = "free_y"
  ) +
  xlim(1, ntop) +
  ylim(0, 100) +
  coord_flip() +
  scale_x_continuous(
    breaks = l$order,
    labels = l$covariate,
    expand = c(0, 0)
  )

dev.off()
dev.off()

# END
