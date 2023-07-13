# Big test

# 1 Get ready

library(raster)
library(magrittr)
library(dplyr)
library(caret)
library(readr)
library(tidyr)
library(lubridate)
library(foreach)
library(tibble)
library(stringr)
library(doParallel)
library(ggupset) # for combination plots
library(forcats) # reordering factors

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

source("loadandstack.R")

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

pts_top_tex <- pts_top@data[, names(pts_top) %in% texclas]

texclas <- names(pts_top_tex)

# 3 Load covariates

## 3.1 Basic covariates (redo extraction if necessary)

cov_basic <- paste0(root, "/cov_basic") %>% loadandstack()

extr_basic <- raster::extract(cov_basic, pts_top)

# extr_basic <- extr_basic[, -c(1:2)]

# extr_basic %>% write.table(file = 'dsc_top_extr_basic.csv'
#                      , sep = ';'
#                      , row.names = FALSE
#                      )

extr_basic <- read.table(
  file = "dsc_top_extr_basic.csv",
  sep = ";",
  header = TRUE
)

extr_basic$geology_10m %<>% as.factor
extr_basic$georeg_10m %<>% as.factor
extr_basic$landscapes_10m %<>% as.factor

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

tr_grid <- data.frame(
  committees = 1,
  neighbors = 0
)

cov_grid <- expand_grid(
  basic = 0:1,
  bare = 0:1,
  series = 0:1,
  radar = 0:1
) %>%
  filter(!row_number() %in% 1)

pts_top_tex$CACO3 %<>% replace_na(0)

all_imps <- list()
all_RMSEs <- matrix(numeric(),
  ncol = ncol(pts_top_tex),
  nrow = nrow(cov_grid),
  dimnames = list(
    1:nrow(cov_grid),
    colnames(pts_top_tex)
  )
)

all_R2s <- matrix(numeric(),
  ncol = ncol(pts_top_tex),
  nrow = nrow(cov_grid),
  dimnames = list(
    1:nrow(cov_grid),
    colnames(pts_top_tex)
  )
)

for (i in 1:nrow(cov_grid))
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

  for (k in 1:ncol(pts_top_tex))
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

    cl <- makePSOCKcluster(12)
    registerDoParallel(cl)

    model_k <- caret::train(fm_k,
      data = dat,
      method = "cubist",
      tuneGrid = tr_grid,
      trControl = trainControl(
        method = "cv",
        number = 12
      )
    )

    registerDoSEQ()
    rm(cl)

    all_imps[[i]][[k]] <- model_k %>% varImp()

    all_RMSEs[i, k] <- model_k$results$RMSE
    all_R2s[i, k] <- model_k$results$Rsquared

    all_RMSEs %>%
      write.table(
        file = "bigtest1_rmses.csv",
        sep = ";",
        row.names = FALSE
      )

    all_R2s %>%
      write.table(
        file = "bigtest1_r2s.csv",
        sep = ";",
        row.names = FALSE
      )
  }

  outname <- paste0(
    "bigtest1_imps_",
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

# 5 Present results

plots <- root %>%
  paste0(., "/bigtest1_plots/")

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
  file = "bigtest1_importance.pdf",
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
