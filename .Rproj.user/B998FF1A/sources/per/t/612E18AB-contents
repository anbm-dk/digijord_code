# 19: Partial dependency plots

library(terra)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(sf)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

testn <- 14
mycrs <- "EPSG:25832"

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") 

fractions_alt  <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")
fractions      <- fractions_alt
frac_ind_predict <- c(1:length(fractions))[-3]  # Exclude fine sand

obs <- paste0(dir_results, "/observations_texture.rds") %>%
  readRDS()

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

# 1 Load importance

imp_all <- read.table(
  paste0(dir_results, "/varImp_all.csv"),
  sep = ";",
  header = TRUE
)

# 2 Analyze correlation between covariates based on random points

# 2.1 Extract points from tiles in parallel

# cells_per_pt <- 2000
# 
# dir_tiles <- dir_dat %>%
#   paste0(., "/tiles_591/")
# 
# subdir_tiles <- dir_tiles %>%
#   list.dirs() %>%
#   .[-1]
# 
# library(parallel)
# 
# numCores <- detectCores()
# numCores
# 
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
#     "subdir_tiles",
#     "dir_dat",
#     "cells_per_pt"
#   )
# )
# 
# pts_tiles <- parSapplyLB(
#   cl,
#   1:length(subdir_tiles),
#   function(j) {
#     set.seed(1)
#     
#     cov_x_files <- subdir_tiles[j] %>%
#       list.files(full.names = TRUE)
#     
#     cov_x_names <- cov_x_files %>%
#       basename() %>%
#       file_path_sans_ext()
#     
#     cov_x <- cov_x_files %>% rast()
#     
#     names(cov_x) <- cov_x_names
#     
#     cov_pts_x <- terra::spatSample(
#       x = cov_x,
#       size = ncell(cov_x) / cells_per_pt,
#       na.rm = FALSE,
#       # as.points = TRUE,
#       xy = TRUE,
#       exp = 1
#     )
#     
#     return(cov_pts_x)
#   },
#   simplify = FALSE,
#   USE.NAMES = FALSE
# )
# 
# stopCluster(cl)
# foreach::registerDoSEQ()
# rm(cl)
# 
# cov_pts <- bind_rows(pts_tiles) %>%
#   na.omit() %>%
#   sample_n(
#     min(
#       10^5,
#       nrow(.)
#     )
#   )

# cov_pts %>%
#   saveRDS(
#     file = paste0(dir_results, "/cov_pts_100k.rds")
#   )

cov_pts <- paste0(dir_results, "/cov_pts_100k.rds") %>% readRDS()

# 2.2 Analyze correlation for top 20

top20_names <- imp_all %>%
  filter(
    !(Covariate %in% grep('fuzzy', Covariate, value = TRUE)),
    !(Covariate %in% c("upper", "lower", "SOM_removed")),
    !(Covariate %in% grep('ogc_pi', Covariate, value = TRUE)),
  ) %>%
  select(Covariate) %>%
  slice_head(n = 20) %>%
  unlist() %>%
  unname()

cor_top20 <- cov_pts %>%
  select(any_of(top20_names)) %>%
  cor()

# 2.3 Find important uncorrelated covariates

# imp_selected <- imp_all %>%
#   filter(
#     !(Covariate %in% grep('fuzzy', Covariate, value = TRUE)),
#     !(Covariate %in% c("upper", "lower", "SOM_removed")),
#     !(Covariate %in% grep('ogc_pi', Covariate, value = TRUE)),
#     is.finite(mean)
#   )

imp_selected <- imp_all %>%
  filter(
    !(Covariate %in% grep('fuzzy', Covariate, value = TRUE)),
    !(Covariate %in% c("upper", "lower", "SOM_removed")),
    # !(Covariate %in% grep('ogc_pi', Covariate, value = TRUE)),
    is.finite(mean)
  )

top_i <- imp_selected %>%
  slice_head(n = 1) %>%
  select(Covariate) %>%
  unlist() %>%
  unname()

for (i in 1:19) {
  imp_i <- imp_selected %>%
    filter(!(Covariate %in% top_i)) %>%
    mutate(imp100 = mean/max(mean, na.rm = TRUE))
  
  pts_top_i <- cov_pts %>% select(any_of(top_i))
  pts_imp_i <- cov_pts %>% select(any_of(imp_i$Covariate))
  
  max_cor <- cor(
    pts_imp_i,
    pts_top_i
  ) %>%
    # "^"(2) %>%
    abs() %>%
    apply(., 1, max) %>%
    unname()
  
  imp_i %<>%
    mutate(
      max_cor = max_cor,
      imp_new = imp100*(1 - max_cor)
    ) %>%
    arrange(-imp_new)
  
  top_i <- imp_i %>%
    slice_head(n = 1) %>%
    select(Covariate) %>%
    unlist() %>%
    unname() %>%
    c(top_i, .)
}

top_i

# Everything:
# "cwl_10m_fuzzy"                          "fuzzy_geology_7"                       
# "filled_s1_baresoil_composite_vh_8_days" "fuzzy_georeg_3"                        
# "vdtochn"                                "ogc_pi906"                             
# "wetlands_10m_fuzzy"                     "fuzzy_landscape_11"                    
# "dhm2015_terraen_10m"                    "filled_s2_geomedian_b11"               
# "fuzzy_geology_5"                        "fuzzy_landscape_9"                     
# "chelsa_bio17_1981_2010_10m"             "filled_s2_geomedian_b3"                
# "flooded_depth_10m_mean"                 "s2_count_max10_fuzzy"                  
# "filled_s1_baresoil_composite_vv_8_days" "ogc_pi438"                             
# "chelsa_bio15_1981_2010_10m"             "slope" 

# No categorical covariates
# "filled_s1_baresoil_composite_vh_8_days" "vdtochn"                               
# "chelsa_bio17_1981_2010_10m"             "dhm2015_terraen_10m"                   
# "filled_s2_geomedian_b11"                "flooded_depth_10m_mean"                
# "chelsa_bio15_1981_2010_10m"             "filled_s2_geomedian_b3"                
# "ogc_pi906"                              "filled_s1_baresoil_composite_vv_8_days"
# "ogc_pi375"                              "valley_depth"                          
# "slope"                                  "filled_s2_geomedian_b12"               
# "sin_aspect_radians"                     "s2_geomedian_20180701_20180731_b4"     
# "ogc_pi125"                              "cos_aspect_radians"                    
# "s2_geomedian_20180408_20180515_b3"      "standardized_height" 

# No categorical covariates or OGCs:
#  "filled_s1_baresoil_composite_vh_8_days" "vdtochn"                               
#  "chelsa_bio17_1981_2010_10m"             "dhm2015_terraen_10m"                   
#  "filled_s2_geomedian_b11"                "flooded_depth_10m_mean"                
#  "chelsa_bio15_1981_2010_10m"             "filled_s2_geomedian_b3"                
#  "filled_s1_baresoil_composite_vv_8_days" "valley_depth"                          
#  "slope"                                  "filled_s2_geomedian_b12"               
#  "sin_aspect_radians"                     "s2_geomedian_20180701_20180731_b4"     
#  "cos_aspect_radians"                     "s2_geomedian_20180408_20180515_b3"     
#  "standardized_height"                    "rvb_bios"                              
#  "chelsa_bio16_1981_2010_10m"             "s2_geomedian_20180701_20180731_b8" 

# 2.4 Good candidates:

candidates_top <- c(
  "filled_s1_baresoil_composite_vh_8_days",
  "vdtochn",
  "ogc_pi906",
  "dhm2015_terraen_10m",
  "filled_s2_geomedian_b11",
  "chelsa_bio17_1981_2010_10m",
  "filled_s2_geomedian_b3",
  "flooded_depth_10m_mean",
  "s2_count_max10_fuzzy",
  "filled_s1_baresoil_composite_vv_8_days"
)

# Other variables to check:

candidates_other <- c(
  "ogc_pi438",
  "chelsa_bio15_1981_2010_10m",
  "slope",
  "valley_depth",
  "filled_s2_geomedian_b12",
  "sin_aspect_radians",
  "s2_geomedian_20180701_20180731_b4",
  "cos_aspect_radians",
  "s2_geomedian_20180408_20180515_b3",
  "standardized_height"
)

# Probably not worth checking:
# "rvb_bios" (categorical)
# "chelsa_bio16_1981_2010_10m" (probably just a proxy for geology)
# "s2_geomedian_20180701_20180731_b8" (ok effect for coarse sand)

# 3 Make PDPs for all selected covariates for each fraction

# 3.0 Trial plots using observations

var_i <- 10

var1 <- "s2_geomedian_20180701_20180731_b8"
# var1 <- candidates_top[var_i]
# var1 <- candidates_other[var_i]


# obs %>%
#   select(any_of(c(fractions, var1))) %>%
#   tidyr::pivot_longer(
#     cols = any_of(fractions),
#     names_to = "fraction",
#     values_to = "value",
#     values_drop_na = TRUE
#     ) %>%
#   mutate(fraction = factor(fraction, levels = fractions)) %>%
#   ggplot(
#     aes(x = .data[[var1]], y = value)
#   ) +
#   geom_smooth() +
#   facet_wrap(vars(fraction))


# Not good enough
# "vdtochn" (big effect close to channels)
# "chelsa_bio17_1981_2010_10m" (probably just a proxy for geology)
# "flooded_depth_10m_mean" (skewed variable)
# "s2_count_max10_fuzzy" (mainly important for SOC)
# "chelsa_bio15_1981_2010_10m" (probably just a proxy for geology)
# "slope" (skewed variable)
# "valley_depth" (big effect at 0)
# "sin_aspect_radians" (no clear effect)
# "cos_aspect_radians" (no clear effect)
# "standardized_height" (ok, but the effect is not strong)
# "s2_geomedian_20180701_20180731_b8" (ok, but the effect is not strong)

# Good for plotting:
goodforplot <- c(
  "filled_s1_baresoil_composite_vh_8_days",
  "ogc_pi906",
  "dhm2015_terraen_10m",
  "filled_s2_geomedian_b11",
  "filled_s2_geomedian_b3",
  "filled_s1_baresoil_composite_vv_8_days",
  "ogc_pi438",
  "filled_s2_geomedian_b12",
  "s2_geomedian_20180701_20180731_b4",
  "s2_geomedian_20180408_20180515_b3"
)

for_plot_bare_soil <- c(
  "filled_s1_baresoil_composite_vh_8_days",
  "filled_s1_baresoil_composite_vv_8_days",
  "filled_s2_geomedian_b3",
  "filled_s2_geomedian_b11",
  "filled_s2_geomedian_b12"
)

for_plot_other <- c(
  "dhm2015_terraen_10m",
  "ogc_pi438",
  "ogc_pi906",
  "s2_geomedian_20180408_20180515_b3",
  "s2_geomedian_20180701_20180731_b4"
)


# Bare soil plot

# obs %>%
#   select(any_of(c(fractions, for_plot_bare_soil))) %>%
#   select(-fine_sand) %>%
#   tidyr::pivot_longer(
#     cols = any_of(fractions),
#     names_to = "fraction",
#     values_to = "value_y",
#     values_drop_na = TRUE
#   ) %>%
#   pivot_longer(
#     cols = any_of(for_plot_bare_soil),
#     names_to = "Covariate",
#     values_to = "value_x",
#     values_drop_na = TRUE
#   ) %>%
#   mutate(fraction = factor(fraction, levels = fractions)) %>%
#   ggplot(
#     aes(x = value_x, y = value_y)
#   ) +
#   geom_smooth() +
#   facet_grid(fraction ~ Covariate, scales = "free")

# Other variables

# obs %>%
#   select(any_of(c(fractions, for_plot_other))) %>%
#   select(-fine_sand) %>%
#   tidyr::pivot_longer(
#     cols = any_of(fractions),
#     names_to = "fraction",
#     values_to = "value_y",
#     values_drop_na = TRUE
#   ) %>%
#   pivot_longer(
#     cols = any_of(for_plot_other),
#     names_to = "Covariate",
#     values_to = "value_x",
#     values_drop_na = TRUE
#   ) %>%
#   mutate(fraction = factor(fraction, levels = fractions)) %>%
#   ggplot(
#     aes(x = value_x, y = value_y)
#   ) +
#   geom_smooth() +
#   facet_grid(fraction ~ Covariate, scales = "free")

# Get quantiles for covariates

cov_pts_q <- cov_pts %>%
  lapply(
    function(x) {
      out <- stats::quantile(
        x,
        probs = seq(0.01, 0.99, 0.02),
        # seq(0.05, 0.95, 0.1),
        na.rm = TRUE
      )
      return(out)
    }
  ) %>%
  bind_cols()

# 3.1 Plots using models

# Load models

library(stringr)
library(pdp)
library(caret)

nboot_final <- 100

nboot_max <- 100
# nboot_max <- 100

nboot <- min(c(nboot_final, nboot_max))

boot_number_chr <- c(1:nboot) %>%
  str_pad(
    .,
    nchar(nboot_final),
    pad = "0"
  )

boot_all_chr <- boot_number_chr %>%
  paste0("boot_", .)

# pdp_outlist <- list()
# 
# for (i in frac_ind_predict) {
#   print(fractions[i])
#   for (bootr in 1:nboot) {
#     print(bootr)
#     model_i <- models_boot_files[[i]][bootr] %>% readRDS()
#     names_model_i <- varImp(model_i)$importance %>% rownames()
#     for (k in 1:length(goodforplot)) {
#       if (goodforplot[k] %in% names_model_i) {
#         print(goodforplot[k])
#         pg_k <- cov_pts_q %>%
#           select(any_of(goodforplot[k]))
#         outgrid <- partial(
#           model_i,
#           pred.var = goodforplot[k],
#           pred.grid = pg_k
#         )
#         names(outgrid) <- c("x", "yhat")
#         outgrid$Fraction <- fractions[i]
#         outgrid$Covariate <- goodforplot[k]
#         pdp_outlist[[length(pdp_outlist) + 1]] <- outgrid
#       }
#     }
#   }
# }
# 
# pdp_df_raw <- pdp_outlist %>%
#   bind_rows()
# 
# saveRDS(
#   pdp_df_raw,
#   file = paste0(dir_results, "/pdp_df_raw.rds")
#   )

pdp_df_raw <- readRDS(paste0(dir_results, "/pdp_df_raw.rds"))

pdp_df <- pdp_df_raw %>%
  group_by(Fraction, Covariate, x) %>%
  summarise(
    mean = mean(yhat),
    low = stats::quantile(yhat, probs = 0.25),
    high = stats::quantile(yhat, probs = 0.75)
    )

frac_labels <- c(
  "Clay",
  "Silt",
  "Fine sand",
  "Coarse<br>sand",
  "SOC",
  "CaCO<sub>*3*</sub>"
)

bare_soil_labels <- c(
  "S1 VH - Bare (dB)",  # "filled_s1_baresoil_composite_vh_8_days"
  "S1 VV - Bare (dB)",  # "filled_s1_baresoil_composite_vv_8_days"
  "S2 B3 - Bare",  # "filled_s2_geomedian_b3"                 
  "S2 B11 - Bare",  # "filled_s2_geomedian_b11"               
  "S2 B12 - Bare"  # "filled_s2_geomedian_b12"  
)

equal_breaks <- function(n = 3, s = 0.1, ...){
  function(x){
    d <- s * diff(range(x)) / (1 + 2*s)
    seq = seq(min(x) + d, max(x) - d, length = n)
    round(seq, -floor(log10(abs(seq[2] - seq[1]))))
  }
}

equal_breaks2 <- function(n = 3, s = 0.1, ...) {
  function(x){
    diff(range(x))
    lowest <- min(x) + diff(range(x))*s
    highest <- max(x) - diff(range(x))*s
    round1 <- -floor(log10(abs(highest - lowest)))
    lowest2 <- ceiling(lowest*10^round1) / 10^round1
    highest2 <- floor(highest*10^round1) / 10^round1
    if (highest2 == lowest2) {
      round1 <- round1 + 1
      lowest2 <- ceiling(lowest*10^round1) / 10^round1
      highest2 <- floor(highest*10^round1) / 10^round1
    }
    if (n > 2) {
      seq1 <- seq(lowest2, highest2, length = n)
      intermediates <- seq1[-c(1, length(seq1))]
      intermediates <- round(intermediates, round1 + 1)
      seq1 <- c(lowest2, intermediates, highest2)
    } else {
      seq1 <- c(lowest2, highest2)
    }
    seq1
  }
}

library(ggtext)
library(extrafont)
font_import()
loadfonts(device = "win")

tiff(
  paste0(dir_results, "/pdp_texture_bare_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 600
)

pdp_df %>%
  filter(Covariate %in% for_plot_bare_soil) %>%
  mutate(
    Fraction = factor(
      Fraction,
      levels = fractions,
      labels = frac_labels
      ),
    Covariate = factor(
      Covariate,
      levels = for_plot_bare_soil,
      labels = bare_soil_labels
    ),
    x = case_when(
      Covariate %in% bare_soil_labels[1:2] ~ x / 100,
      .default = x / 10^4
      ) 
    ) %>%
  ggplot(
    aes(x = x, y = mean)
  ) +
  geom_line(aes(x = x, y = low), colour = "cadetblue") +
  geom_line(aes(x = x, y = high), colour = "cadetblue") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "cadetblue") +
  geom_line() +
  facet_grid(
    Fraction ~ Covariate, scales = "free"
    ) +
  theme(
    text = element_text(family = "serif"),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 9),
    axis.text.y = element_markdown(size = 9),
    strip.text.y = element_markdown(size = 9, lineheight = 1),
    strip.text.x = element_markdown(size = 9, lineheight = 1)
  ) +
  ylab("ŷ (%)") +
  scale_x_continuous(
    expand = c(0, 0),
    # n.breaks = 4,
    breaks = equal_breaks2(n = 2, s = 0.1)
  ) +
  scale_y_continuous(breaks = equal_breaks2(n = 2, s = 0.1))

try(dev.off())

other_plot_labels <- c(
  "Elevation (m)",   # "dhm2015_terraen_10m"              
  "OGC 0.438π (km)",    # "ogc_pi438"                        
  "OGC 0.906π (km)",     # "ogc_pi906"                        
  "S2 B3 - Spring",  # "s2_geomedian_20180408_20180515_b3"
  "S2 B4 - July"    # "s2_geomedian_20180701_20180731_b4"
)

tiff(
  paste0(dir_results, "/pdp_texture_other_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 600
)

pdp_df %>%
  filter(Covariate %in% for_plot_other) %>%
  mutate(Fraction = factor(
    Fraction,
    levels = fractions,
    labels = frac_labels
  ),
  Covariate = factor(
    Covariate,
    levels = for_plot_other,
    labels = other_plot_labels
  ),
  x = case_when(
    Covariate == "OGC 0.438π (km)" ~ x / 10^3,
    Covariate == "OGC 0.906π (km)" ~ x / 10^3,
    Covariate == "S2 B3 - Spring" ~ x / 10^4,
    Covariate == "S2 B4 - July" ~ x / 10^4,
    .default = x
  )
  ) %>%
  ggplot(
    aes(x = x, y = mean)
  ) +
  geom_line(aes(x = x, y = low), colour = "cadetblue") +
  geom_line(aes(x = x, y = high), colour = "cadetblue") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "cadetblue") +
  geom_line() +
  facet_grid(Fraction ~ Covariate, scales = "free") +
  theme(
    text = element_text(family = "serif"),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 9),
    # axis.text.x = element_markdown(angle = 270, size = 9, vjust = 0.5, hjust = 1),
    axis.text.y = element_markdown(size = 9),
    strip.text.y = element_markdown(size = 9, lineheight = 1),
    strip.text.x = element_markdown(size = 9, lineheight = 1)
  ) +
  ylab("ŷ (%)") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = equal_breaks2(n = 2, s = 0.10)
  ) +
  scale_y_continuous(breaks = equal_breaks2(n = 2, s = 0.10))

try(dev.off())

# Save models

# library(xgboost)
# 
# dir_boot_models_raw <- dir_boot %>%
#   paste0(., "/models_raw/") %T>%
#   dir.create()
# 
# for (i in 1:length(fractions)) {
#   print(fractions[i])
#   dir_boot_models_raw_i <- dir_boot_models_raw %>%
#     paste0(., fractions[i], "/") %T>%
#     dir.create()
#   for (bootr in 1:nboot) {
#     print(bootr)
#     
#     model_i <- models_boot_files[[i]][bootr] %>% readRDS()
#     
#     xgb_i <- xgb.Booster.complete(model_i$finalModel)
#     
#     outname_i <- paste0(
#       dir_boot_models_raw_i,
#       "/xgbmodel_", fractions[i], "_", boot_number_chr[bootr], ".model"
#     )
#     
#     xgb.save(xgb_i, outname_i)
#   }
# }

# Plots with depth

library(xgboost)
library(pdp)
library(tibble)

dir_boot_models_raw <- dir_boot %>%
  paste0(., "/models_raw/")

xvalues <- seq(
  from = quantile(
    obs$filled_s1_baresoil_composite_vh_8_days,
    probs = 0.01, na.rm = TRUE, names = FALSE),
  to = quantile(
    obs$filled_s1_baresoil_composite_vh_8_days,
    probs = 0.99, na.rm = TRUE, names = FALSE),
  length.out = 32
)

pgrid <- expand.grid(
  upper = seq(0, 190, 10),
  filled_s1_baresoil_composite_vh_8_days = xvalues
) %>%
  mutate(lower = upper + 10)

predvars_full <- c("upper", "lower", "filled_s1_baresoil_composite_vh_8_days")

plotratio <- ((max(pgrid[, 2]) - min(pgrid[, 2])) / 200) * (20 / 32)

# pdp_outlist_depth <- list()
# 
# for (i in frac_ind_predict) {
#   for (bootr in 1:nboot) {
#     
#     model_ib <- models_boot_files[[i]][bootr] %>% readRDS()
#     
#     names_model_i <- varImp(model_ib)$importance %>% rownames()
#     
#     print(paste(fractions[i], "model", bootr))
#     
#     ok1 <- c("upper", "lower") %>%
#       is_in(names_model_i) %>%
#       sum() %>%
#       is_greater_than(0)
#     ok2 <- "filled_s1_baresoil_composite_vh_8_days" %>%
#       is_in(names_model_i) %>%
#       add(ok1) %>%
#       is_greater_than(1)
#     
#     if (ok2) {    
#       predvars_ib <- predvars_full %>%
#         magrittr::extract(. %in% names_model_i)
#       
#       pgrid_ib <- pgrid %>% select(any_of(predvars_ib))
#     
#       p1xv <- pdp::partial(
#         model_ib,
#         pred.var = predvars_ib,
#         pred.grid = pgrid_ib,
#         type = "regression"
#       )
#       
#       p1xv %<>%
#         add_column(
#           !!!pgrid[setdiff(names(pgrid), names(.))]
#         ) %>%
#         mutate(
#           Depth = (upper + lower)/2,
#           Fraction = fractions[i]
#           )
#         
#       pdp_outlist_depth[[length(pdp_outlist_depth) + 1]] <- p1xv
#     }
#   }
# }
# 
# pdp_depth_raw <- pdp_outlist_depth %>%
#   bind_rows()
# 
# saveRDS(
#   pdp_depth_raw,
#   file = paste0(dir_results, "/pdp_depth_raw.rds")
#   )

pdp_depth_raw <- readRDS(file = paste0(dir_results, "/pdp_depth_raw.rds"))

pdp_depth_mean <- pdp_depth_raw %>%
  group_by(Depth, filled_s1_baresoil_composite_vh_8_days, Fraction) %>%
  summarise(yhat = mean(yhat, na.rm = TRUE)) %>%
  ungroup()

pdp_depth_mean %<>%
  arrange(filled_s1_baresoil_composite_vh_8_days)

S1_VH_spring <- rep(
  pgrid$filled_s1_baresoil_composite_vh_8_days,
  each = nrow(pdp_depth_mean) / nrow(pgrid)
)

pdp_depth_mean$filled_s1_baresoil_composite_vh_8_days <- S1_VH_spring

frac_labels2 <- c(
  "Clay (%)",
  "Silt (%)",
  "Fine sand (%)",
  "Coarse sand (%)",
  "SOC (%)",
  "CaCO<sub>*3*</sub> (%)"
)

pdp_plotlist1 <- list()

for (i in 1:length(frac_ind_predict)) {
  pdp_plotlist1[[i]] <- pdp_depth_mean %>%
    mutate(
      S1_VH_spring = filled_s1_baresoil_composite_vh_8_days/100
    ) %>%
    filter(Fraction == fractions[frac_ind_predict[i]]) %>%
    ggplot() +
    geom_raster(aes(x = S1_VH_spring, y = Depth, fill = yhat)) +
    scale_fill_viridis_c(name = "ŷ", breaks = equal_breaks2(n = 2, s = 0.10)) +
    coord_fixed(ratio = plotratio/100, expand = FALSE) +
    scale_y_reverse(breaks = c(0, 100, 200)) +
    scale_x_continuous(breaks = equal_breaks2(n = 2, s = 0.10)) +
    labs(title = paste0(LETTERS[i], ": ", frac_labels2[frac_ind_predict[i]])) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_markdown(margin = unit(c(0, 0, -0.05, 0), "cm")),
      legend.title = ggtext::element_markdown(),
      legend.key.size = unit(0.5, "line"),
      legend.box.spacing = unit(1, "line"),
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, -10, -10, -10),
      plot.background = element_rect(fill = 'transparent', color = NA),
      plot.margin = margin(0, 0, 0, 0, "cm")
      ) +
    guides(col = guide_colourbar(title = "ŷ"))
}

library(grid)
library(ggpubr)

figure <- ggarrange(
    pdp_plotlist1[[1]] + rremove("ylab") + rremove("xlab") + rremove("x.text"),
    pdp_plotlist1[[2]] +
      rremove("ylab") +
      rremove("xlab") + 
      rremove("x.text") + 
      rremove("y.text"),
    pdp_plotlist1[[3]] + rremove("ylab") + rremove("xlab") + rremove("y.text"),
    pdp_plotlist1[[4]] + rremove("ylab") + rremove("xlab"),
    pdp_plotlist1[[5]] + rremove("ylab") + rremove("xlab") + rremove("y.text"),
    align = "hv",
    nrow = 2,
    ncol = 3
    # ,
    # labels = fractions[frac_ind_predict]
) +
  theme(plot.margin = margin(0.1, 0.5, 0, 0.25, "cm")) 

figure

tiff(
  paste0(dir_results, "/pdp_texture_depth_VH_test_", testn, ".tiff"),
  width = 16,
  height = 6.7,
  units = "cm",
  res = 600
)

annotate_figure(
  figure,
  left = text_grob(
    "Depth (cm)",
    rot = 90, vjust = 1,
    # gp = gpar(cex = 1.3),
    family = "serif"
    ),
  bottom = text_grob(
    "S1 VH - Bare (dB)",
    # gp = gpar(cex = 1.3),
    family = "serif"
    )
)

try(dev.off())
try(dev.off())

# Additional depth plots - test these:

# Use these:
# clay / filled_s2_geomedian_b12 [OK]
# clay / dhm2015_terraen_10m [OK]
# silt / dhm2015_terraen_10m [OK]
# silt / s2_geomedian_20180408_20180515_b3 [OK]
# coarse_sand / dhm2015_terraen_10m [OK]
# coarse_sand / s2_geomedian_20180601_20180630_b8 [OK]
# SOC / s2_count_max10_fuzzy [OK] - good for plot
# SOC / filled_s2_geomedian_b3 [OK] - effect in the topsoil
# CaCO3 / vdtochn [OK] - best candidate

# Dropped
# clay / s2_geomedian_20180408_20180515_b8a [ - ]    # Use dhm instead
# SOC / vdtochn [OK] - strong effect close to channels
# SOC / s2_geomedian_20180701_20180731_b4 [ - ] - same as bare soil effect
# CaCO3 / standardized_height [OK]
# (CaCO3 / dhm2015_terraen_10m) [ - ]

combinations_pdp <- data.frame(
  fraction = c(
    "clay",
    "clay",
    "silt",
    "silt",
    "coarse_sand",
    "coarse_sand",
    "SOC",
    "SOC",
    "CaCO3"
  ),
  covaritate = c(
    "filled_s2_geomedian_b12",
    "dhm2015_terraen_10m",
    "dhm2015_terraen_10m",
    "s2_geomedian_20180408_20180515_b3",
    "dhm2015_terraen_10m",
    "s2_geomedian_20180601_20180630_b8",
    "s2_count_max10_fuzzy",
    "filled_s2_geomedian_b3",
    "vdtochn"
  )
)

cov_pts_q2 <- cov_pts %>%
  lapply(
    function(x) {
      out <- seq(
        from = quantile(
          x,
          probs = 0.01, na.rm = TRUE, names = FALSE),
        to = quantile(
          x,
          probs = 0.99, na.rm = TRUE, names = FALSE),
        length.out = 32
      )
      return(out)
    }
  ) %>%
  bind_cols()

# pdp_outlist_combo <- list()
# 
# # loop
# for (j in 1:nrow(combinations_pdp)) {
#   i <- which(fractions %in% combinations_pdp$fraction[j])
#   
#   predvars_full_j <- c("upper", "lower", combinations_pdp$covaritate[j])
#   
#   pgrid_j <- expand.grid(
#     upper = seq(0, 190, 10),
#     x = cov_pts_q2 %>%
#       select(any_of(combinations_pdp$covaritate[j])) %>%
#       unlist() %>%
#       unname()
#   ) %>%
#     mutate(lower = upper + 10) %>%
#     rename_with(~ combinations_pdp$covaritate[j], x)
#   
#   plotratio_j <- ((max(pgrid_j[, 2]) - min(pgrid_j[, 2])) / 200) * (20 / 32)
#   
#   pdp_outlist_j <- list()
#   
#   for (bootr in 1:nboot) {
#     
#     model_ib <- models_boot_files[[i]][bootr] %>% readRDS()
#     
#     names_model_i <- varImp(model_ib)$importance %>% rownames()
#     
#     ok1 <- c("upper", "lower") %>%
#       is_in(names_model_i) %>%
#       sum() %>%
#       is_greater_than(0)
#     ok2 <- combinations_pdp$covaritate[j] %>%
#       is_in(names_model_i) %>%
#       add(ok1) %>%
#       is_greater_than(1)
#     
#     if (ok2) {
#       print(
#         paste(
#           fractions[i], combinations_pdp$covaritate[j], "model", bootr,
#           "processing")
#       )
#       
#       predvars_ib <- predvars_full_j %>%
#         magrittr::extract(. %in% names_model_i)
#       
#       pgrid_ib <- pgrid_j %>% select(any_of(predvars_ib))
#       
#       p1xv <- pdp::partial(
#         model_ib,
#         pred.var = predvars_ib,
#         pred.grid = pgrid_ib,
#         type = "regression"
#       )
#       
#       p1xv %<>%
#         add_column(
#           !!!pgrid[setdiff(names(pgrid), names(.))]
#         ) %>%
#         mutate(
#           Depth = (upper + lower)/2
#         )
#       
#       pdp_outlist_j[[length(pdp_outlist_j) + 1]] <- p1xv
#     } else {
#       print(
#         paste(
#           fractions[i], combinations_pdp$covaritate[j], "model", bootr,
#           "dropped")
#       )
#     }
#   }
#   
#   if (length(pdp_outlist_j) > 0) {
#     pdp_depth_mean_j <- pdp_outlist_j %>%
#       bind_rows() %>%
#       group_by(Depth, .data[[combinations_pdp$covaritate[j]]]) %>%
#       summarise(yhat = mean(yhat, na.rm = TRUE)) %>%
#       ungroup() %>%
#       arrange(.data[[combinations_pdp$covaritate[j]]])
#     
#     pdp_depth_mean_j$x <- pgrid_j[, 2]
#     
#     pdp_depth_mean_j %<>%
#       mutate(
#         Fraction = fractions[i],
#         covariate = combinations_pdp$covaritate[j],
#         plotratio = plotratio_j
#       ) %>%
#       select(Fraction, covariate, Depth, x, plotratio, yhat)
#     
#     pdp_outlist_combo[[length(pdp_outlist_combo) + 1]] <- pdp_depth_mean_j
#   } else {
#     pdp_outlist_combo[[length(pdp_outlist_combo) + 1]] <- "no results"
#   }
# }

# saveRDS(
#   pdp_outlist_combo,
#   file = paste0(dir_results, "/pdp_outlist_combo.rds")
#   )

pdp_outlist_combo <- readRDS(paste0(dir_results, "/pdp_outlist_combo.rds"))

frac_labels3 <- c(
  "A: Clay (%)",
  "B: Clay (%)",
  "C: Silt (%)",
  "D: Silt (%)",
  "E: Coarse sand (%)",
  "F: Coarse sand (%)",
  "G: SOC (%)",
  "H: SOC (%)",
  "I: CaCO<sub>*3*</sub> (%)"
)

cov_labels <- c(
  "S2 B12 - Bare",
  "Elevation (m)",
  "Elevation (m)",
  "S2 B3 - Spring",
  "Elevation (m)",
  "S2 B8 - June",
  "Bare images (n)",
  "S2 B3 - Bare",
  "vdchn (m)"
)

s2_indices <- c(1, 0, 0, 1, 0, 1, 0, 1, 0)

xlab_digits <- c(2, 0, 0, 2, 0, 2, 1, 2, 0)

pdp_plotlist2 <- list()

for (i in 1:nrow(combinations_pdp)) {
  pdp_combo_i <- pdp_outlist_combo[[i]]
  if (s2_indices[i]) {
    pdp_combo_i %<>%
      mutate(
        x = x/10^4,
        plotratio = plotratio/10^4
      )
  }
  
  x_breaks <- c(
    pdp_combo_i$x %>%
      min() %>%
      multiply_by(10^xlab_digits[i]) %>%
      ceiling() %>%
      divide_by(10^xlab_digits[i]),
    pdp_combo_i$x %>%
      max() %>%
      multiply_by(10^xlab_digits[i]) %>%
      floor() %>%
      divide_by(10^xlab_digits[i])
  )
  
  pdp_plotlist2[[i]] <- pdp_combo_i %>%
    ggplot() +
    geom_raster(aes(x = x, y = Depth, fill = yhat)) +
    coord_fixed(ratio = pdp_combo_i$plotratio[1], expand = FALSE) +
    scale_y_reverse(breaks = c(0, 100, 200)) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = x_breaks
    ) +
    scale_fill_viridis_c(name = "ŷ", breaks = equal_breaks2(n = 2, s = 0.10)) +
    ggtitle(frac_labels3[i]) +
    xlab(cov_labels[i]) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_markdown(margin = unit(c(0, 0, 0, 0), "cm")),
      legend.title = ggtext::element_markdown(),
      legend.key.size = unit(0.5, "line"),
      legend.box.spacing = unit(1, "line"),
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, -10, -10, -10),
      plot.background = element_rect(fill = 'transparent', color = NA),
      plot.margin = margin(0, 0, 0, 0.1, "cm"),
      axis.title.x = element_text(vjust = +5)
    ) +
    guides(col = guide_colourbar(title = "ŷ"))
}

library(grid)
library(ggpubr)

figure2 <- ggarrange(
  pdp_plotlist2[[1]] + rremove("ylab"),
  pdp_plotlist2[[2]] + rremove("ylab") + rremove("y.text"),
  pdp_plotlist2[[3]] + rremove("ylab") + rremove("y.text"),
  pdp_plotlist2[[4]] + rremove("ylab"),
  pdp_plotlist2[[5]] + rremove("ylab") + rremove("y.text"),
  pdp_plotlist2[[6]] + rremove("ylab") + rremove("y.text"),
  pdp_plotlist2[[7]] + rremove("ylab"),
  pdp_plotlist2[[8]] + rremove("ylab") + rremove("y.text"),
  pdp_plotlist2[[9]] + rremove("ylab") + rremove("y.text"),
  align = "hv",
  nrow = 3,
  ncol = 3
) +
  theme(plot.margin = margin(0.1, 0.5, 0, 0.25, "cm")) 

# figure2

tiff(
  paste0(dir_results, "/pdp_texture_depth_combo_test_", testn, ".tiff"),
  width = 16,
  height = 11,
  units = "cm",
  res = 600
)

annotate_figure(
  figure2,
  left = text_grob(
    "Depth (cm)",
    rot = 90, vjust = 1,
    # gp = gpar(cex = 1.3),
    family = "serif"
  )
)

try(dev.off())
try(dev.off())


# Other potential covariates (for individual fractions)

# Selected:
# vdtochn
# (- silt)
# - SOC  # NB
# - CaCO3  # NB
# s2_count_max10_fuzzy
# - SOC  # NB
# s2_geomedian_20180601_20180630_b8
# - coarse sand  # NB
# s2_geomedian_20180601_20180630_b7
# - coarse sand  # NB
# s2_geomedian_20180601_20180630_b11
# - silt  # NB
 
# Not selected:
# slope (log transform?)
# - SOC  # NB
# - CaCO3  # NB
# flooded_depth_10m_mean
# - SOC  # NB
# - CaCO3  # NB
# valley_depth
# - clay  
# standardized_height
# - coarse sand
# - CaCO3

# - Save the xgboost models to their native format.

# END