# 11: Sampling points for validation fields

library(terra)
library(magrittr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

dir_sampling <- paste0(dir_dat, "/sampling/")

dir_sampling_2023 <- paste0(dir_sampling, "/digijord_marker_2023/")

fields_2023 <- dir_sampling_2023 %>%
  paste0(., "/digijord_marker_2023.shp") %>%
  vect()

# Load covariates

dir_cov <- dir_dat %>% paste0(., "/covariates")

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20230712.csv") %>%
  read.table(
    sep = ";",
    header = TRUE
  )

cov_files <- dir_cov %>% list.files(full.names = TRUE)

cov_names <- cov_files %>% basename() %>% tools::file_path_sans_ext()

covs_use <- cov_cats$name[cov_cats$anbm_use == 1]

cov <- cov_files %>% rast() %>% subset(covs_use)

covs_for_sampling <- c(
  "s2_geomedian_b11", "s2_geomedian_b12", "s2_geomedian_b2", "s2_geomedian_b3",
  "s2_geomedian_b4", "s2_geomedian_b5", "s2_geomedian_b6", "s2_geomedian_b7",
  "s2_geomedian_b8", "s2_geomedian_b8a", "s1_baresoil_composite_vh_8_days",  
  "s1_baresoil_composite_vv_8_days", "dhm2015_terraen_10m"
)

cov_sampling <- subset(cov, covs_for_sampling)

safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)

source("f_sample_kmeans.R")

# Field processing

i <- 1

pts_list <- list()

for(i in 1:length(fields_2023)) {
  covs_i <- crop(
    cov_sampling,
    fields_2023[i],
    mask = TRUE
  )
  
  field_i_rast <- fields_2023[i] %>%
    as.lines() %>%
    rasterize(
      covs_i,
      touches = TRUE
    )
  
  covs_i %<>% mask(field_i_rast, inverse = TRUE)
  pca <- covs_i %>% as.data.frame(na.rm = TRUE) %>% prcomp(scale. = TRUE)
  pcs_i <- predict(covs_i, pca) %>% subset(1:4)
  
  set.seed(1)
  
  pts_list[[i]] <- sample_kmeans(
    input = pcs_i,
    use_xy = TRUE,
    clusters = 9,
    scale = TRUE,
    sp_pts = TRUE,
    filename_cl = paste0(
      dir_sampling_2023, "/", fields_2023$CVR[i] ,"_zoner.tif"
      ),
    args_cl = list(
      overwrite = TRUE,
      datatype = "INT2U"
    ),
    filename_pts = paste0(
      dir_sampling_2023, "/", fields_2023$CVR[i] ,"_punkter.shp"
    ),
    args_pts = list(overwrite = TRUE),
  )
}

pdf(paste0(dir_sampling_2023, "/punkter_2023.pdf"))

for(i in 1:length(pts_list)) {
  
  pts_list[[i]]$clusters %>%
    as.factor() %>%
    plot(
      main = fields_2023$CVR[i],
      col = safe_colorblind_palette
    )
  points(pts_list[[i]]$points, pch = 21, bg = "white")
}

try(dev.off())
try(dev.off())

# END