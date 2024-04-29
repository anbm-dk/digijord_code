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

# END