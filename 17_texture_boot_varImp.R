# Covariate importance plots for bootstrap texture models

library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)

testn <- 14
mycrs <- "EPSG:25832"

# Fractions

fractions_alt  <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC", "CaCO3")
fractions      <- fractions_alt
fraction_names <- c("Clay", "Silt", "Fine sand", "Coarse sand", "SOC", "CaCO3")

# Results folder

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

dir_results <- dir_dat %>%
  paste0(., "/results_test_", testn, "/") %T>%
  dir.create()

dir_boot <- dir_results %>%
  paste0(., "/bootstrap/") %T>%
  dir.create()

# Covariate table

cov_cats <- dir_code %>%
  paste0(., "/cov_categories_20240305.csv") %>%
  read.table(
    sep = ",",
    header = TRUE
  )

l_cat <- cov_cats %>%
  mutate(
    Covariate = name,
    category = case_when(
      category == "basic" ~ scorpan,
      category == "WATEM" ~ "OR",
      category == "sentinel_composite" ~ "S2 time series",
      category == "bare_soil" ~ "Bare soil",
      desc_text %in% grep(
        "bare soil",
        cov_cats$desc_text,
        value = TRUE
      ) ~ "Bare soil", 
      .default = "Other"
    )
  ) %>%
  mutate(
    category = case_when(
      category == "RP" ~ "P",
      category == "CR" ~ "R",
      category == "OR" ~ "R",
      .default = category
    )
  )

l_cat %>% select(name, category)

# Covariate importance results

boot_imp_all <- dir_boot %>%
  paste0(., "/models_boot_importance_all.rds") %>%
  readRDS()

boot_imp_sand <- boot_imp_all %>%
  pivot_longer(
    cols = -c(fraction, Covariate),
    values_to = "Importance",
    names_to = "rep"
  ) %>%
  filter(
    fraction %in% c("clay", "silt", "coarse_sand")
  ) %>%
  group_by(Covariate, rep) %>%
  summarise(Importance = sum(Importance)/3) %>%
  ungroup() %>%
  arrange(-Importance) %>%
  mutate(fraction = "fine_sand")

l <- boot_imp_all %>%
  pivot_longer(
    cols = -c(fraction, Covariate),
    values_to = "Importance",
    names_to = "rep"
  ) %>%
  filter(fraction != "fine_sand") %>%
  bind_rows(boot_imp_sand)


# Join categories

l %<>%
  left_join(l_cat) %>%
  mutate(
    category = case_when(
      Covariate == "upper" ~ "Depth",
      Covariate == "lower" ~ "Depth",
      Covariate == "year" ~ "Time",
      Covariate == "SOM_removed" ~ "Method",
      category == "N" ~ "Spatial position",
      category == "R" ~ "Topography",
      category == "C" ~ "Climate",
      category == "C " ~ "Climate",
      category == "P" ~ "Parent materials",
      category == "S" ~ "Soil",
      category == "SO" ~ "Soil and organisms",
      category == "CR" ~ "Climate and topography",
      category == "OR" ~ "Organisms and topography",
      category == "O" ~ "Organisms",
      category == "RP" ~ "Parent materials",
      .default = category
    ),
    fraction = factor(
      fraction,
      levels = fractions,
      labels = fraction_names
    )
  )

cat_overall_order <- l %>%
  group_by(category) %>%
  summarise(
    mean = sum(Importance)
  ) %>%
  arrange(-mean) %>%
  select(category) %>%
  unlist() %>%
  as.character()

l %<>% mutate(
  category = factor(category, levels = rev(cat_overall_order))
)

levels(l$fraction) <- c(
  expression("Clay"),
  expression("Silt"),
  expression("Fine~sand"),
  expression("Coarse~sand"),
  expression("SOC"),
  expression("CaCO"[3])
)

# Colors for categories

library(colorRamps)
library(rcartocolor) # for colorblind palette

mycolors <- carto_pal(12, "Safe") %>% sort()

library(TSP)

catcolors <- l$category %>%
  levels() %>%
  length() %>%
  carto_pal(., "Safe")
names(catcolors) <- levels(l$category)
colScale <- scale_fill_manual(name = "category", values = catcolors)

# Make plot for categories

tiff(
  paste0(dir_results, "/boot_importance_categories_test_", testn, ".tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 600
)

l %>%
  filter(Importance > 0) %>%
  group_by(category, fraction, rep) %>%
  summarise(Importance = sum(Importance)) %>%
  ungroup() %>%
  group_by(category, fraction) %>%
  mutate(mean = mean(Importance)) %>%
  ggplot(aes(x = category, y = Importance, fill = category)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_violin(scale = "width") +
  facet_wrap(~ fraction, labeller = label_parsed) +
  coord_flip() +
  geom_point(aes(x = category, y = mean, fill = category)) +
  theme(legend.position = "none") +
  xlab("Category")

try(dev.off())

# Plot importance for OGCs

ndir_plot <- 64

imp_ogc <- varImp_boot_mean %>%
  select(-mean_imp) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "target",
    values_to = "Overall"
  ) %>%
  filter(covariate %in% grep('ogc_pi', colnames(obs), value = TRUE)) %>%
  group_by(target) %>%
  mutate(Overall = Overall*100/max(Overall)) %>%
  arrange(covariate) %>%
  mutate(
    dir = substr(
      covariate,
      nchar(covariate) - 2,
      nchar(covariate)
    ) %>%
      as.numeric() %>%
      "*" (ndir_plot) %>%
      "/" (1000) %>%
      "+" (1) %>%
      round(digits = 0)
  ) %>%
  filter(is.finite(Overall))

imp_ogc2 <- imp_ogc %>%
  mutate(dir = dir + ndir_plot)

imp_ogc %<>% bind_rows(., imp_ogc2)

# imp_ogc %<>% mutate(dir = row_number())

imp_ogc %<>% mutate(
  target = factor(
    target,
    levels = fractions,
    labels = fraction_names
  )
)

brks <- seq(
  from = min(imp_ogc$dir),
  by = (max(imp_ogc$dir) + 1 - min(imp_ogc$dir))/4,
  length.out = 4
)

tiff(
  paste0(dir_results, "/boot_importance_ogc_test", testn, ".tiff"),
  width = 40,
  height = 20,
  units = "cm",
  res = 300
)

ggplot(imp_ogc, aes(x = dir, y = Overall)) +
  coord_polar(
    start = - pi/2 - pi/(ndir_plot*2),
    direction = -1) +
  geom_polygon(colour = 'black', fill = rgb(0,2/3,2/3,1/2)) + 
  # geom_col(width = 1, colour = 'black', fill = rgb(0,2/3,2/3,1/2)) +
  facet_wrap(
    ~ target,
    ncol = 3
  ) +
  scale_x_continuous(
    breaks = brks,
    labels = c('E', 'N', 'W', 'S')
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ylab('Covariate importance') +
  theme_bw() +
  theme(axis.text.x = element_text(
    colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.y = element_text(colour = 'black'),
    panel.grid.major = element_line(color = 'grey'),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1)
  )

try(dev.off())

# END