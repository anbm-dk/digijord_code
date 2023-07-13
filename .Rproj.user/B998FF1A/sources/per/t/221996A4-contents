# 08: Clustering of raster cells

# 1: Startup

library(terra)
library(magrittr)
library(graphics)
library(cluster)
library(randomcoloR)
library(dplyr)
library(MultivariateRandomForest)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")
dir_cov <- dir_dat %>% paste0(., "/covariates")

mycrs <- "EPSG:25832"

dem <- dir_cov %>%
  paste0("/dhm2015_terraen_10m.tif") %>%
  rast()

# dem1km <- dem %>%
#   aggregate(
#     fact = 100
#     , na.rm = TRUE
#     , cores = 20
#     , filename = paste0(dir_dat, '/layers/dhm2015_terraen_1km.tif')
#     , overwrite = TRUE
#   )

dem1km <- dir_dat %>%
  paste0(., "/layers/dhm2015_terraen_1km.tif") %>%
  rast()

coords <- crds(dem1km)


# 2: Distance-based clustering for islands

# distances <- dist(coords, method = 'maximum')
#
# clustering <- hclust(distances, method = 'single')
#
# groups <- cutree(clustering, k = 5)
#
# groups_df <- cbind(coords, groups) %>%
#   as.data.frame()
#
# groups_r <- rast(groups_df)
#
# set.seed(1)
#
# myrandomcolors <- randomColor(length(unique(groups)))
#
# plot(groups_r, col = myrandomcolors)


# 3: 4-km raster for MRF

dem4km <- dem1km %>%
  is.na(.) %>%
  ifel(., NA, 1) %>%
  aggregate(
    fact = 4,
    na.rm = TRUE,
    overwrite = TRUE,
    fun = "sum"
  ) %T>%
  plot

dem4km_df <- as.data.frame(dem4km, xy = TRUE, cells = TRUE)

coords_4km <- dem4km_df %>%
  select(c(x, y)) %>%
  as.matrix()


# Multivariate decision tree

# nrows <- nrow(coords_4km)
# # my_minleaf <- nrows/40  # 60 tiles
# my_minleaf <- nrows/300  # 440 tiles
#
# # nrows <- 1000
#
# mytree <- build_single_tree(
#   coords_4km[1:nrows, ],
#   coords_4km[1:nrows, ],
#   min_leaf = my_minleaf,
#   m_feature = 2,
#   Command = 2,
#   Inv_Cov_Y = solve(cov(coords_4km[1:nrows, ]))
#   )
#
# mypred <- single_tree_prediction(
#   mytree,
#   coords,
#   2
#   )
#
# mypred <- paste(mypred[, 1], mypred[, 2])
#
# mydf <- as.data.frame(coords)
#
# mydf$pred <- as.factor(mypred) %>% as.numeric()
#
# set.seed(1)
#
# myrandomcolors2 <- randomColor(length(unique(mydf$pred)))
#
# mydf %>% rast %>% plot(col = myrandomcolors2)
#
# # Turn groups into polygons
#
# rast_mygroups <- rast(mydf)
#
# crs(rast_mygroups) <- mycrs
#
# group_polygons <- rast_mygroups %>%
#   as.polygons() %>%
#   disagg()
#
# exts <- list()
#
# for (i in 1:length(group_polygons)) {
#   exts[[i]] <- ext(group_polygons[i]) %>% as.polygons(crs = mycrs)
#   exts[[i]]$pred <- group_polygons$pred[i]
# }
#
# group_boxes <- vect(exts)
#
# crs(group_boxes) <- mycrs
#
# group_boxes_dissolved <- terra::aggregate(group_boxes, by = "pred")
#
# group_polygons2 <- group_boxes_dissolved %>%
#   disagg()
#
# exts2 <- list()
#
# for (i in 1:length(group_polygons2)) {
#   exts2[[i]] <- ext(group_polygons2[i]) %>% as.polygons(crs = mycrs)
#   exts2[[i]]$pred <- group_polygons2$pred[i]
# }
#
# group_boxes2 <- vect(exts2)
#
# group_boxes_dissolved2 <- terra::aggregate(group_boxes2, by = "pred")
#
# group_polygons3 <- group_boxes_dissolved2 %>%
#   disagg()
#
# exts3 <- list()
#
# for (i in 1:length(group_polygons3)) {
#   exts3[[i]] <- ext(group_polygons3[i]) %>% as.polygons(crs = mycrs)
#   exts3[[i]]$pred <- group_polygons3$pred[i]
#   exts3[[i]]$id <- i
# }
#
# group_boxes3 <- vect(exts3)
#
# terra::union(group_boxes3) # no more overlaps
#
# crs(group_boxes3) <- mycrs
#
# myrandomcolors_boxes <- randomColor(length(group_boxes3))
#
# plot(is.na(dem1km))
# plot(group_boxes3, add = TRUE)
# plot(group_boxes3, "id", col = myrandomcolors_boxes)
# plot(final_boxes, add = TRUE)

# Square tiles

lapply(1:30, function(x) {
  out <- aggregate(dem1km, x, na.rm = TRUE) %>%
    as.data.frame() %>%
    nrow()
  return(out)
}) %>%
  unlist() %>%
  data.frame()

# External calculation indicates optimum at 11 km side

dem11km <- aggregate(dem1km, 11, na.rm = TRUE)

dem11_polygons <- dem11km %>%
  as.data.frame(xy = TRUE, cells = TRUE) %>%
  select(c(x, y, cell)) %>%
  rast() %>%
  as.polygons() %T>% plot()

dem1km_pol <- dem1km %>%
  as.polygons() %>%
  aggregate()

terra::intersect(dem11_polygons, dem1km_pol) %>% plot()

group_polygons_loop <- terra::intersect(dem11_polygons, dem1km_pol) %>% disagg()

for (j in 1:3) {
  exts <- list()

  for (i in 1:length(group_polygons_loop)) {
    exts[[i]] <- ext(group_polygons_loop[i]) %>% as.polygons(crs = mycrs)
    exts[[i]]$cell <- group_polygons_loop$cell[i]
  }

  group_boxes <- vect(exts)

  crs(group_boxes) <- mycrs

  group_boxes_dissolved <- terra::aggregate(group_boxes, by = "cell")

  group_polygons_loop <- group_boxes_dissolved %>%
    disagg()
}

terra::union(group_polygons_loop) # no more overlaps

crs(group_polygons_loop) <- mycrs

groups_merged <- group_polygons_loop %>%
  aggregate(dissolve = TRUE) %>%
  disagg()

groups_merged$id <- 1:nrow(groups_merged)

exts <- list()

for (i in 1:length(groups_merged)) {
  exts[[i]] <- ext(groups_merged[i]) %>% as.polygons(crs = mycrs)
  exts[[i]]$cell <- groups_merged$id[i]
}

groups_merged_boxes <- vect(exts)

crs(groups_merged_boxes) <- mycrs

groups_merged_boxes$diagonals_m <- sqrt(
  perim(groups_merged_boxes)^2 - 8 * expanse(groups_merged_boxes)
) / 2

max_diagonal <- 1556 * 10

groups_merged_boxes <- subset(
  groups_merged_boxes,
  diagonals_m < max_diagonal,
  NSE = TRUE
)

plot(groups_merged_boxes, "diagonals_m")

groups_merged_boxes$id <- 1000 * (1:nrow(groups_merged_boxes))

group_polygons_loop$id <- 1:nrow(group_polygons_loop)

square_groups <- cover(group_polygons_loop, groups_merged_boxes)

square_groups$id <- 1:nrow(square_groups)

square_groups <- subset(square_groups, select = id, NSE = TRUE)

myrandomcolors_boxes <- randomColor(length(square_groups))

plot(is.na(dem1km))
plot(square_groups, add = TRUE)
plot(square_groups, col = myrandomcolors_boxes)
plot(square_groups, "id")

# Write groups to shapefile

dir_tiles <- dir_dat %>%
  paste0(., "/tiles_", length(square_groups), "/") %T>%
  dir.create()

writeVector(
  square_groups,
  filename = paste0(dir_tiles, "tiles.shp"),
  filetype = "ESRI Shapefile"
)

# # Combine groups
#
# df_all <- mydf %>%
#   select(pred) %>%
#   bind_cols(groups_df, .) %>%
#   mutate(
#     group_final = ifelse(
#       groups == 1,
#       pred,
#       groups*100
#     )
#   ) %>% select(
#     c(x, y, group_final)
#   ) %>%
#   mutate(
#     group_final = group_final %>% as.factor() %>% as.numeric()
#   )
#
# myrandomcolors3 <- randomColor(max(df_all$group_final))
#
# df_all %>% rast() %>% as.factor() %>% plot(col = myrandomcolors3)

# # Turn combined groups into polygons
#
# rast_all <- rast(df_all)
#
# crs(rast_all) <- mycrs
#
# polygons_all <- as.polygons(rast_all) %T>% plot("group_final")
#
# exts <- list()
#
# for (i in 1:length(unique(df_all$group_final))) {
#   exts[[i]] <- ext(polygons_all[i]) %>% as.polygons(crs = mycrs)
# }
#
# final_boxes <- vect(exts)
#
# plot(rast_all, col = myrandomcolors3)
# plot(final_boxes, add = TRUE)
#
# # Write combined groups to shapefile
#
# dir_tiles <- dir_dat %>%
#   paste0(., "/tiles_", length(myrandomcolors3), "/") %T>%
#   dir.create()
#
# writeVector(
#   final_boxes,
#   filename = paste0(dir_tiles, "tiles.shp"),
#   filetype = "ESRI Shapefile"
# )

# END
