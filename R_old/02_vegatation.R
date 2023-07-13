# Test correlation between soil and vegation

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

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

source("loadandstack.R")

# Load DSC points

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

# Load vegetation rasters

## Store extracted values here

extr_files <- paste0(root, "/sat_extracts")

extr_files %>% dir.create()

## Find the structure of the satellite images

folders <- paste0(root, "/digijord_download") %>%
  list.dirs(recursive = FALSE) %>%
  head(-1)

parts <- c("DK", "BH", "DK", "BH")

rdirs <- paste0(root, "/digijord_download") %>%
  list.dirs()

## All dates

all_dates <- rdirs %>%
  stringr::str_extract("[0-9]{4}[0-9]{2}[0-9]{2}") %>%
  lubridate::ymd() %>%
  unique() %>%
  sort()

## Routine for loop

for (j in 1:length(folders))
{
  rdirs_j <- folders[j] %>%
    list.dirs()

  for (i in 1:length(rdirs_j))
  {
    rfiles <- rdirs_j[i] %>%
      list.files(
        pattern = "tif$",
        full.names = TRUE
      )

    if (length(rfiles) > 0) {
      rs <- loadandstack(rdirs_j[i])

      extr <- raster::extract(rs,
        pts_top,
        sp = FALSE,
        df = TRUE
      ) %>%
        dplyr::select(-1)

      rdate <- names(rs)[1] %>%
        stringr::str_extract("[0-9]{4}[0-9]{2}[0-9]{2}") %>%
        lubridate::ymd()

      outname <- paste0(extr_files, "/extr_", rdate, "_", parts[j], ".csv")

      extr %>% write.table(
        file = outname,
        sep = ";",
        row.names = FALSE
      )
    }
  }
}
rm(outname, extr, rdate, i, j, rdirs_j, rfiles, rs)

# Analysis

## Count NAs

all_extrs <- extr_files %>% list.files(full.names = TRUE)

not_nas <- numeric()

for (i in 1:length(all_extrs))
{
  t_i <- read.table(all_extrs[i],
    sep = ";",
    header = TRUE
  )

  if (ncol(t_i) == 10 &
    sum(t_i, na.rm = TRUE) > 0) {
    not_nas[i] <- apply(t_i, 1, FUN = function(x) {
      sum(x) %>%
        is.na() %>%
        isFALSE()
    }) %>% sum()
  } else {
    not_nas[i] <- 0
  }
}

## NAs per point

pts_nas <- matrix(0,
  ncol = length(all_extrs),
  nrow = nrow(extr)
)

for (i in 1:length(all_extrs))
{
  t_i <- read.table(all_extrs[i],
    sep = ";",
    header = TRUE
  )

  if (ncol(t_i) == 10 &
    sum(t_i, na.rm = TRUE) > 0) {
    pts_nas[, i] <- apply(t_i, 1, FUN = function(x) {
      sum(x) %>%
        is.na() %>%
        isFALSE()
    })
  }
}

na_pts <- pts_nas %>% apply(1, FUN = function(x) {
  out <- sum(x) == 0
  return(out)
})

na_pts_sp <- pts_top[na_pts, ]

plot(pts_top)
plot(na_pts_sp, col = "red", add = TRUE)

all_extrs_dates <- all_extrs %>%
  stringr::str_extract("[0-9]{4}\\-[0-9]{2}\\-[0-9]{2}") %>%
  lubridate::ymd()

weeks_df <- data.frame(
  year = year(all_extrs_dates),
  week = isoweek(all_extrs_dates)
)

weeks_simp <- paste0(weeks_df$year, weeks_df$week)

weeks_uniq <- unique(weeks_simp)

f <- function(i) {
  dat <- pts_nas[, weeks_simp == i]
  if (is.matrix(dat)) {
    out_i <- apply(dat, 1, FUN = function(x) {
      out <- sum(x, na.rm = TRUE) > 0
      out %<>% as.numeric()
      return(out)
    })
  } else {
    out_i <- dat %>% as.numeric()
  }
  return(out_i)
}

nas_week <- foreach(i = weeks_uniq) %do% f(i)

lapply(nas_week, sum)

length(nas_week)

nas_biweek <- numeric()

for (i in 2:length(nas_week))
{
  x <- nas_week[[i]] + nas_week[[i - 1]]
  x <- x == 0
  nas_biweek[i] <- sum(x)
}

nas_week %>%
  lapply(sum) %>%
  unlist() %>%
  plot()
nas_biweek %>% points(col = "red")

nas_triweek <- numeric()

for (i in 3:length(nas_week))
{
  x <- nas_week[[i]] + nas_week[[i - 1]] + nas_week[[i - 2]]
  x <- x == 0
  nas_triweek[i] <- sum(x)
}


nas_week %>%
  lapply(sum) %>%
  unlist() %>%
  plot()
nas_biweek %>% points(col = "red")
nas_triweek %>% points(col = "blue")

cbind(weeks_uniq, nas_biweek)

nas_week %>%
  lapply(sum) %>%
  unlist()


# Merge extracts

extrs_merged <- paste0(root, "/sat_extracts_merged")

extrs_merged %>% dir.create()

f_read <- function(x) {
  out <- read.table(x,
    sep = ";",
    header = TRUE
  )
  return(out)
}

f_null <- function(x) {
  if (is.data.frame(x)) {
    if (ncol(x) == 10 &
      sum(x, na.rm = TRUE) > 0) {
      out <- TRUE
    } else {
      out <- FALSE
    }
  } else {
    out <- FALSE
  }
  return(out)
}

for (i in 1:length(all_dates))
{
  xs <- all_extrs[all_extrs_dates == all_dates[i]]

  l <- lapply(xs, f_read)

  l_ok <- lapply(l, f_null) %>% unlist()

  l <- l[l_ok]

  if (length(l) > 0) {
    extr_i <- l[[1]]
    if (length(l) == 2) {
      ind <- extr_i %>% is.na()
      extr_i[ind] <- l[[2]][ind]
    }

    outname <- paste0(extrs_merged, "/extr_", all_dates[i], ".csv")

    extr_i %>% write.table(
      file = outname,
      sep = ";",
      row.names = FALSE
    )
  }
}


# Check correlation

extrs_merged_files <- extrs_merged %>% list.files(full.names = TRUE)

extrs_merged_dates <- extrs_merged_files %>%
  stringr::str_extract("[0-9]{4}\\-[0-9]{2}\\-[0-9]{2}") %>%
  lubridate::ymd()

pts_top_tex$CACO3[is.na(pts_top_tex$CACO3)] <- 0

l_cors <- list()

for (i in 1:length(extrs_merged_files))
{
  extr_i <- read.table(extrs_merged_files[i],
    sep = ";",
    header = TRUE
  )

  l_cors[[i]] <- cor(
    x = pts_top_tex,
    y = extr_i,
    use = "pairwise.complete.obs"
  )
}

l_cors_na <- l_cors %>%
  lapply(FUN = function(x) is.na(x) %>% sum()) %>%
  unlist()

l_cors <- l_cors[l_cors_na == 0]
l_cors_dates <- extrs_merged_dates[l_cors_na == 0]

maxcors <- l_cors %>%
  lapply(max) %>%
  unlist()
mincors <- l_cors %>%
  lapply(min) %>%
  unlist()

plot(x = l_cors_dates, y = maxcors)
plot(x = l_cors_dates, y = mincors)

cors_all <- l_cors %>%
  lapply(function(x) as_tibble(x)) %>%
  lapply(function(x) {
    add_column(x,
      texclass = texclas
    )
  }) %>%
  lapply(function(x) {
    pivot_longer(x,
      cols = !texclass,
      names_to = "date_band"
    )
  }) %>%
  bind_rows() %>%
  add_column(date = .$date_band %>%
    stringr::str_extract("[0-9]{4}[0-9]{2}[0-9]{2}") %>%
    lubridate::ymd()) %>%
  add_column(band = .$date_band %>%
    str_sub(-3, -1))

library(ggplot2)

pdf(
  file = "bands_texture.pdf",
  paper = "a4",
  width = 0,
  height = 0
)

cors_all %>% ggplot(aes(x = band, y = value)) +
  geom_violin(fill = "green") +
  geom_boxplot(width = 0.1) +
  facet_grid(rows = vars(texclass))

dev.off()
dev.off()

cors_all %>%
  group_by(texclass, band) %>%
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = mean,
    names_from = band
  )

cors_all %>%
  group_by(texclass, band) %>%
  summarise(mean = mean(value^2, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = mean,
    names_from = band
  )

tclasses <- c("LER", "SILT", "GROVSAND", "HUMUS")
bands <- c("B11", "B11", "B11", "B04")

cors_all$year <- cors_all$date %>% year()
cors_all$doy <- cors_all$date %>% yday()

cor_plots <- list()

for (i in 1:length(tclasses))
{
  cor_plots[[i]] <- cors_all %>%
    dplyr::filter(texclass == tclasses[i] & band == bands[i]) %>%
    ggplot(aes(x = doy, y = value)) +
    geom_hline(yintercept = 0) +
    geom_point(
      shape = 21,
      fill = "white"
    ) +
    geom_smooth(method = "loess") +
    facet_wrap(vars(year), ncol = 2) +
    ggtitle(paste0(LETTERS[i], ": Correlation between ", tclasses[i], " and ", bands[i]))
}

library(patchwork)

tiff(
  filename = "correlation_bands.tiff",
  units = "cm",
  width = 16 * 1.5,
  height = 10 * 1.5,
  res = 600
)

(cor_plots[[1]] | cor_plots[[2]]) / (cor_plots[[3]] | cor_plots[[4]])

dev.off()
dev.off()


# Test model

f_nas <- function(x) {
  d <- read.table(x,
    sep = ";",
    header = TRUE
  )
  out <- is.na(d[, 1]) %>% mean()
  return(out)
}

n_nas <- extrs_merged_files %>%
  sapply(f_nas) %>%
  unname()

plot(x = extrs_merged_dates, y = n_nas)

extrs_formodel <- extrs_merged_files[n_nas < 0.5]
dates_formodel <- extrs_merged_dates[n_nas < 0.5]

f_read2 <- function(x) {
  out <- read.table(x,
    sep = ";",
    header = TRUE
  )
  return(out)
}

extrs_formodel %<>%
  lapply(f_read2) %>%
  bind_cols()

extrs_formodel[is.na(extrs_formodel)] <- 0

models_veg <- list()

for (i in 1:length(texclas))
{
  dat <- cbind(
    pts_top@data[, names(pts_top) == texclas[i]],
    extrs_formodel
  ) # !

  dat %<>% as.data.frame
  colnames(dat)[1] <- texclas[i]
  dat %<>% drop_na

  fm_i <- names(extrs_formodel) %>% # !
    paste(collapse = " + ") %>%
    paste0(texclas[i], " ~ ", .) %>%
    as.formula()

  library(doParallel)
  cl <- makePSOCKcluster(12)
  registerDoParallel(cl)

  models_veg[[i]] <- caret::train(fm_i,
    data = dat,
    method = "cubist",
    trControl = trainControl(
      method = "cv",
      number = 3
    )
  )

  registerDoSEQ()
  rm(cl)
}

get_r2 <- function(x) {
  out <- x$results$Rsquared %>% max()
  return(out)
}

r2s <- models_veg %>%
  lapply(get_r2) %>%
  unlist()

get_RMSE <- function(x) {
  out <- x$results$RMSE %>% min()
  return(out)
}

RMSEs <- models_veg %>%
  lapply(get_RMSE) %>%
  unlist()

getimps <- function(x) {
  d <- x$importance
  d$image <- rownames(d)
  d %<>% arrange(image)
  return(d)
}

imps <- models_veg %>%
  lapply(varImp) %>%
  lapply(getimps) %>%
  bind_rows() %>%
  add_column(tclass = texclas %>% rep(each = 200)) %>%
  add_column(date = .$image %>%
    stringr::str_extract("[0-9]{4}[0-9]{2}[0-9]{2}") %>%
    lubridate::ymd()) %>%
  add_column(band = .$image %>%
    str_sub(-3, -1))

imps %>%
  ggplot(aes(x = band, y = Overall)) +
  geom_violin(fill = "green") +
  geom_boxplot(width = 0.1) +
  facet_grid(rows = vars(tclass))


imps %>%
  group_by(tclass, band) %>%
  summarise(mean = mean(Overall, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = mean,
    names_from = band
  )

imps %>%
  group_by(tclass, band) %>%
  summarise(max = max(Overall, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = max,
    names_from = band
  )

imps %>%
  group_by(tclass, date) %>%
  summarise(mean = mean(Overall, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = mean,
    names_from = date
  ) %>%
  write.table(
    file = "imps_veg_date.csv",
    sep = ";"
  )


imps %>%
  group_by(tclass, date) %>%
  summarise(max = max(Overall, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = max,
    names_from = date
  ) %>%
  write.table(
    file = "imps_veg_date_max.csv",
    sep = ";"
  )

imps %>%
  group_by(band, date) %>%
  summarise(mean = mean(Overall, na.rm = TRUE)) %>%
  pivot_wider(
    values_from = mean,
    names_from = date
  ) %>%
  write.table(
    file = "imps_veg_date_band.csv",
    sep = ";"
  )

imps_wide <- imps %>% pivot_wider(
  values_from = Overall,
  names_from = tclass
)

imps_wide %>% write.table(
  file = "imps_veg_all.csv",
  sep = ";",
  row.names = FALSE
)

# END
