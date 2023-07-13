# Process radar information
# Anders BM

library(raster)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)
library(tibble)
library(stringr)

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

originals <- "D:/anbm/digijord/Radar/Mosaics"

orig_files <- originals %>% list.files(
  pattern = ".tif",
  full.names = TRUE
)

orig_names <- originals %>%
  list.files(pattern = ".tif") %>%
  strsplit("[.]") %>%
  unlist() %>%
  extract(c(TRUE, FALSE))

b <- orig_files[1] %>% brick()

coast_shp <- shapefile("D:/anbm/digijord/Coast/Kyst.shp")

coast_shp$OBJECTID %<>% as.numeric()

dem <- "D:/anbm/digijord/cov_basic/DHM2015_terraen_10m.tif" %>% raster()

# coast_r <- raster::rasterize(x = coast_shp
#                              , y = dem
#                              , filename = paste0(root, '/Coast/Kyst.tif')
#                              , datatype = 'INT2U'
#                              , overwrite = TRUE)

coast_r <- paste0(root, "/Coast/Kyst.tif") %>% raster()

pol_names <- c("VV", "VH")

originals_all <- orig_files %>% stack()

# 1 Crop to extent

extended <- extend(originals_all,
  y = coast_r
)

cropped <- crop(extended,
  y = coast_r,
  filename = paste0(root, "/Temp/cropped.tif")
)

radar_out <- paste0(root, "/Radar/processed/")

radar_out %>% dir.create()

outnames <- lapply(
  orig_names,
  function(x) {
    out <- paste0(x, "_", pol_names)
  }
) %>%
  unlist() %>%
  paste0(radar_out, ., ".tif")

# 2 Mask and round

s <- stack(coast_r, cropped)

f <- function(x) {
  if (sum(x) %>% is.na()) {
    out <- base::rep(NA, length(x) - 1)
  } else {
    out <- x[-1] %>% round(digits = 1)
  }
  return(out)
}

beginCluster(12)

masked <- clusterR(s,
  calc,
  args = list(fun = f),
  filename = paste0(root, "/Temp/masked.tif"),
  datatype = "FLT4S",
  overwrite = TRUE,
  format = "GTiff"
)

endCluster()

# 3 Write as separate layers

# raster::writeRaster(masked
#                     , filename = outnames
#                     , bylayer = TRUE
#                     , overwrite = TRUE
#                     , format = 'GTiff'
#                     )

# 4 Extract to points

source("loadandstack.R")

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

radar <- radar_out %>% loadandstack()

# extr_radar <- raster::extract(radar, pts_top)
#
# extr_radar %>% write.table(file = 'dsc_top_extr_radar.csv'
#                           , sep = ';'
#                           , row.names = FALSE
# )

extr_radar <- read.table(
  file = "dsc_top_extr_radar.csv",
  sep = ";",
  header = TRUE
)

cor(pts_top_tex,
  extr_radar,
  use = "pairwise.complete.obs"
)

pts_radar <- pts_top

pts_radar@data <- cbind(pts_top@data, extr_radar)

# END
