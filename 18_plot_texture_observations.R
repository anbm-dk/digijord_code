# 18: Plot texture observations

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

obs <- paste0(dir_results, "/observations_texture.rds") %>%
  readRDS()

mask_LU <- paste0(dir_dat, "/layers/Mask_LU.tif") %>% rast()  # reference

raster3km <- terra::aggregate(mask_LU, fact = 300, fun = "max")

coastline <- paste0(dir_dat, "/layers/Kyst/Kyst.shp") %>%
  st_read() %>%
  st_set_crs(mycrs)

# base_gg <- ggplot() +
#   geom_sf(data = coastline, fill = "grey", linewidth = 0, color = NA)

base_gg <- ggplot() +
  geom_sf(data = coastline, fill = NA, linewidth = 0.5, color = "black")

obs_v <- obs %>%
  vect(geom = c("UTMX", "UTMY"), crs = mycrs, keepgeom = TRUE) %>%
  filter(imputed == FALSE)

autoplot(obs_v)

obs_dens_pts <- obs_v %>%
  rasterize(raster3km, fun = function(i) {length(i)}, background = 0) %>%
  as.points() %>%
  filter(V1 > 0) %>%
  mutate(
    dens = V1,
    cat = case_when(
      dens < 25 ~ "A",
      dens < 50 ~ "B",
      dens < 100 ~ "C",
      dens < 250 ~ "D",
      dens < 500 ~ "E",
      dens < 1500 ~ "F",
      TRUE ~ "G"
    )
  )

# First plot

tiff(
  paste0(dir_results, "/texture_obs_density_map_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 600
)

base_gg +
  # Layer, this object is a SpatVector instead of sf object
  geom_spatvector(
    obs_dens_pts,
    mapping = aes(size = cat),
    # Use a point with a border AND a fill
    pch = 21, color = "white", fill = "black", linewidth = 0.05
  ) +
  scale_size_manual(
    values = c(0.3, 1, 1.25, 1.5, 2, 2.5, 3.5),
    labels = c(
      "< 25", "[25, 50)", "[50, 100)", "[100, 250)",
      "[250, 500)", "[500, 1500)", "≥ 1500"
    ),
    guide = guide_legend(
      nrow = 1,
      title.position = "top",
      keywidth = 1,
      label.position = "bottom"
    )
  ) +
  labs(
    size = "Observations per 3x3 km"
  ) +
  theme_void() +
  theme(
    text = element_text(family = "serif"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.margin = margin(r = 20, unit = "pt"),
    legend.key.width = unit(50, "pt"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

try(dev.off())

# Hexagonal grid

target_area <- cellSize(raster3km, unit = "m") %>%
  pull() %>%
  mean()

target_area

area_ratio_hex <- target_area / (9*10^6)

# Infer diam hex
diam_hex <- sqrt(2 * target_area / sqrt(3))
# Create hexagonal grid
pop_agg_sf <- st_make_grid(coastline, cellsize = diam_hex, square = FALSE)

counts_hex <- st_intersects(pop_agg_sf, as_sf(obs_v)) %>% lengths()

pop_sf_points <- pop_agg_sf %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  vect() %>%
  geom() %>%
  as.data.frame() %>%
  select(x, y) %>%
  mutate(
    count = counts_hex,
    dens = count / area_ratio_hex
    ) %>%
  filter(dens > 0) %>%
  vect(geom = c("x", "y"), crs = mycrs, keepgeom = TRUE) %>%
  select(dens) %>%
  # Categorize
  mutate(cat = case_when(
    dens < 25 ~ "A",
    dens < 50 ~ "B",
    dens < 100 ~ "C",
    dens < 250 ~ "D",
    dens < 500 ~ "E",
    dens < 1500 ~ "F",
    TRUE ~ "G"
  ))

size_mult <- 0.5

# base_gg <- ggplot() +
#   geom_sf(data = coastline, fill = "grey", linewidth = 0, color = NA)

base_gg <- ggplot() +
  geom_sf(
    data = coastline,
    fill = NA,
    linewidth = 0.5*size_mult,
    color = "black"
    )

tiff(
  paste0(dir_results, "/texture_obs_density_hex_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 600
)

base_gg +
  # Layer, this object is a SpatVector instead of sf object
  geom_spatvector(
    pop_sf_points,
    mapping = aes(size = cat),
    # Use a point with a border AND a fill
    pch = 21, color = "white", fill = "black", linewidth = 0.01*size_mult
  ) +
  scale_size_manual(
    values = c(0.3, 1, 1.25, 1.5, 2, 2.5, 3.5)*size_mult,
    labels = c(
      "< 25", "[25, 50)", "[50, 100)", "[100, 250)",
      "[250, 500)", "[500, 1500)", "≥ 1500"
    ),
    guide = guide_legend(
      nrow = 1,
      title.position = "top",
      keywidth = 1,
      label.position = "bottom"
    )
  ) +
  labs(
    size = "Observations per 9 km^2"
  ) +
  theme_void() +
  theme(
    text = element_text(family = "serif"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.margin = margin(r = 20, unit = "pt"),
    legend.key.width = unit(50, "pt")
  )

try(dev.off())

# Process observations

breaks <- c(0, 30, 60, 100, 200)

obs_v_long <- obs_v %>%
  values() %>%
  rowwise() %>%
  mutate(
    Texture = is.finite(
      mean(
        c(clay, silt, fine_sand, coarse_sand),
        na.rm = TRUE
        )
      )*1,
    SOC = is.finite(SOC)*1,
    CaCO3 = is.finite(CaCO3)*1,
    Depth = cut((upper + lower)/2, breaks, include.lowest = TRUE)
  ) %>%
  ungroup() %>%
  select(Texture, SOC, CaCO3, Depth, UTMX, UTMY) %>%
  pivot_longer(
    cols = c(Texture, SOC, CaCO3),
    values_to = "value",
    names_to = "Type"
  ) %>%
  mutate(
    Type = factor(Type, c("Texture", "SOC", "CaCO3"))
  ) %>%
  filter(
    value == 1,
    is.finite(Depth)
    ) %>%
  select(-value) %>%
  vect(geom = c("UTMX", "UTMY"), crs = mycrs, keepgeom = TRUE)

# autoplot(obs_v_long) +
#   facet_grid(Depth ~ Type)

# Hex grid with specific size

target_area <- 100*10^6

area_ratio_hex <- target_area / 10^6

# Infer diam hex
diam_hex <- sqrt(2 * target_area / sqrt(3))
# Create hexagonal grid points
pop_agg_sf <- st_make_grid(
  coastline,
  cellsize = diam_hex,
  square = FALSE
)

pop_agg_sf_pts <- pop_agg_sf %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  vect() %>%
  geom() %>%
  as.data.frame() %>%
  select(x, y)

# Analyse each combination

combination_grid <- expand.grid(
  depth = levels(obs_v_long$Depth),
  type = levels(obs_v_long$Type)
)

obs_dens_combinations <- list()

for (i in 1:nrow(combination_grid)) {
  obs_v_long_i <- obs_v_long %>%
    filter(
      Depth == combination_grid$depth[i],
      Type == combination_grid$type[i]
    )
  
  counts_hex_i <- st_intersects(
    pop_agg_sf,
    as_sf(obs_v_long_i)
  ) %>%
    lengths()
  
  obs_dens_combinations[[i]] <- pop_agg_sf_pts %>%
    mutate(
      count = counts_hex_i,
      Depth = combination_grid$depth[i],
      Type = combination_grid$type[i]
      ) %>%
    filter(count > 0) %>%
    vect(geom = c("x", "y"), crs = mycrs, keepgeom = TRUE) %>%
    select(count, Type, Depth, x, y)
}

obs_dens_combinations %<>% vect()

obs_dens_combinations

obs_dens_combinations %<>% mutate(
  # dens = count / area_ratio_hex
  dens = count
)

n_categories <- 5

get_unique_quantiles <- function(
  x = NULL,
  n = n_categories + 1,
  max_n = 20
) {
  test_n <- c(3:max_n)
  calc_unique_q <- function(x2) {
    out2 <- seq(0, 1, length.out = x2) %>%
      quantile(x = x, probs = .) %>%
      unique() %>%
      length()
    return(out2)
  }
  n_quantiles <- sapply(test_n, calc_unique_q)
  out <- quantile(
    x = x,
    probs = seq(
      0,
      1,
      length.out = min(test_n[n_quantiles == n]))
  ) %>%
    unique()
  return(out)
}

myquantiles <- get_unique_quantiles(
  obs_dens_combinations$dens,
  n_categories + 1
)

obs_dens_combinations %<>% mutate(
  cat = cut(dens, myquantiles, include.lowest = TRUE)
) %>%
  arrange(x - y)

levels(obs_dens_combinations$Type) <- c(
  expression("Texture"),
  expression("SOC"),
  expression("CaCO"[3])
)

levels(obs_dens_combinations$Depth) %<>% sapply(deparse) %>% unname()

# Maps for fractions and depths

size_mult <- 0.4

point_sizes <- c(0.3, 1, 1.25, 1.5, 2, 2.5, 3.5) %>%
  rev() %>%
  magrittr::extract(1:n_categories) %>%
  rev() %>%
  "*"(size_mult)

# base_gg <- ggplot() +
#   geom_sf(data = coastline, fill = "grey", linewidth = 0, color = NA)

# base_gg <- ggplot() +
#   geom_sf(
#     data = coastline,
#     fill = NA,
#     linewidth = 0.5*size_mult,
#     color = "black"
#   )

base_gg <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "white",
    color = NA
  )

myplot <- base_gg +
  # Layer, this object is a SpatVector instead of sf object
  geom_spatvector(
    obs_dens_combinations,
    mapping = aes(size = cat),
    # Use a point with a border AND a fill
    pch = 21, color = "white", fill = "black", stroke = 1*size_mult
  ) +
  scale_size_manual(
    values = point_sizes
  ) +
  labs(
    size = "Observations (n)"
  ) +
  theme_void() +
  theme(
    text = element_text(family = "serif", size = 12),
    panel.spacing.x = unit(0.05, "lines"),
    panel.spacing.y = unit(0.05 ,"lines"),
    panel.background = element_rect(fill = "grey80", color = NA),
    strip.text.y.left = element_text(
      angle = 90,
      margin = unit(rep(3, 4), "pt"),
      size = 12
      ),
    strip.text.x.top = element_text(
      margin = unit(rep(3, 4), "pt"),
      size = 12
      ),
    axis.title.y = element_text(angle = 90)
  ) +
  facet_grid(Depth ~ Type, switch = "y", , labeller = label_parsed) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(
    expand = expansion(mult = 0.02)
  ) +
  labs(y = "Depth (cm)")

tiff(
  paste0(dir_results, "/texture_obs_density_hex_depth_test_", testn, ".tiff"),
  width = 16,
  height = 13,
  units = "cm",
  res = 600
)

print(myplot)

try(dev.off())

pdf(
  paste0(dir_results, "/texture_obs_density_hex_depth_test_", testn, ".pdf"),
    width = 16/2.54,
    height = 13/2.54
  )

print(myplot)

try(dev.off())
try(dev.off())

# END