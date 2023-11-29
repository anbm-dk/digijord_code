# Effect of SOM removal before texture measurements

library(terra)
library(magrittr)
library(tools)
library(dplyr)
library(caret)
library(tibble)
library(tidyr)
library(tidyterra)
library(ggplot2)
library(akima) # For ggplot2 interpolation
library(viridis)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")

dir_fig <- dir_dat %>% paste0(., "/figures/") %T>% dir.create()

mycrs <- "EPSG:25832"

# 2: Load observations

dir_obs_proc <- dir_dat %>%
  paste0(., "/observations/processed/")

dsc <- dir_obs_proc %>%
  paste0(., "dsc.csv") %>%
  read.table(
    header = TRUE,
    sep = ";"
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

SEGES <- dir_obs_proc %>%
  paste0(., "SEGES.csv") %>%
  read.table(
    header = TRUE,
    sep = ";"
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

profiles_texture <- dir_obs_proc %>%
  paste0(., "profiles_texture.csv") %>%
  read.table(
    header = TRUE,
    sep = ";"
  ) %>%
  vect(
    geom = c("UTMX", "UTMY"),
    crs = mycrs,
    keepgeom = TRUE
  )

# Test 1: Profiles

profiles_analysis <- profiles_texture %>%
  values() %>%
  mutate(
    HUMUS = SOC / 0.587,
    fine_sum = clay + silt,
    clay_fines = clay / (fine_sum)
  ) %>%
  filter(
    HUMUS < 10,
    is.finite(clay_fines)
  )

profiles_sum <- profiles_analysis %>%
  mutate(hum_bins = cut(HUMUS, breaks = c(0:10))) %>%
  group_by(hum_bins) %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
    ) %>%
  rownames_to_column(var = "bin_num") %>%
  tidyr::extract(hum_bins, c("min", "max"), "\\((.*),(.*)]") %>%
  mutate(
    across(where(is.character), as.numeric)
  ) %>%
  drop_na()

prof_clayplot <- ggplot() +
  geom_point(
    data = profiles_analysis,
    aes(
      x = HUMUS,
      y = clay
      ),
    col = rgb(1, 0, 0, 0.5)
    ) +
  geom_errorbarh(
    data = profiles_sum,
    aes(
      xmin = min,
      xmax = max,
      y = clay_q,
      col = q
    ), 
    linewidth = 1
  ) +
  xlab("SOM (%)") +
  ylab("Clay (%)")

tiff(
  paste0(dir_fig, "/profiles_clay.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(prof_clayplot)
dev.off()

prof_siltplot <- ggplot() +
  geom_point(
    data = profiles_analysis,
    aes(
      x = HUMUS,
      y = silt
    ),
    col = rgb(1, 0, 0, 0.5)
  ) +
  geom_errorbarh(
    data = profiles_sum,
    aes(
      xmin = min,
      xmax = max,
      y = silt_q,
      col = q
    ), 
    linewidth = 1
  ) +
  xlab("SOM (%)") +
  ylab("Silt (%)")

tiff(
  paste0(dir_fig, "/profiles_silt.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(prof_siltplot)
dev.off()

prof_finesplot <- ggplot() +
  geom_point(
    data = profiles_analysis,
    aes(
      x = HUMUS,
      y = fine_sum
    ),
    col = rgb(1, 0, 0, 0.5)
  ) +
  geom_errorbarh(
    data = profiles_sum,
    aes(
      xmin = min,
      xmax = max,
      y = fines_q,
      col = q
    ), 
    linewidth = 1
  ) +
  xlab("SOM (%)") +
  ylab("Clay (%) + Silt (%)")

tiff(
  paste0(dir_fig, "/profiles_fines.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(prof_finesplot)
dev.off()

prof_clayfinesplot <- ggplot() +
  geom_point(
    data = profiles_analysis,
    aes(
      x = HUMUS,
      y = clay_fines
    ),
    col = rgb(1, 0, 0, 0.5)
  ) +
  geom_errorbarh(
    data = profiles_sum,
    aes(
      xmin = min,
      xmax = max,
      y = clayfines_q,
      col = q
    ), 
    linewidth = 1
  ) +
  xlab("SOM (%)") +
  ylab("Clay [%] / (Clay [%] + Silt [%])")


tiff(
  paste0(dir_fig, "/profiles_clayfines.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(prof_clayfinesplot)
dev.off()

qs_45 <- profiles_sum %>%
  filter(min == 4) %>%
  select(c(clay_q, silt_q, fines_q, clayfines_q))
qs_56 <- profiles_sum %>%
  filter(min == 5) %>%
  select(c(clay_q, silt_q, fines_q, clayfines_q))

prof_dif_rel <- qs_56 / qs_45
prof_dif_abs <- qs_56 - qs_45

round(qs_45, 3) %>% as.data.frame()
round(qs_56, 3) %>% as.data.frame()
round(prof_dif_rel, 3) - 1
round(prof_dif_abs, 3)


# Test 2: SEGES vs dsc comparison

dsc_analysis <- dsc %>%
  mutate(
    HUMUS = SOC / 0.587,
    fine_sum = clay + silt,
    clay_fines = clay / (fine_sum)
  ) %>% 
  filter(
    upper == 0,
    HUMUS < 10,
    is.finite(clay_fines)
  )

seges_analysis <- SEGES %>%
  mutate(
    HUMUS = SOC / 0.587,
    fine_sum = clay + silt,
    clay_fines = clay / (fine_sum)
  ) %>% 
  filter(
    upper == 0,
    HUMUS < 10,
    is.finite(clay_fines)
  )

distances <- distance(seges_analysis, dsc_analysis)

mindist_dsc <- apply(distances, 2, min)
mindist_seges <- apply(distances, 1, min)

dsc_selected <- dsc_analysis %>%
  filter(mindist_dsc < 500)
seges_selected <- seges_analysis %>%
  filter(mindist_seges < 500)

dsc_summary <- dsc_selected %>%
  values %>%
  reframe(
  clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
  silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
  fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
  clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
  q = c(0.25, 0.5, 0.75)
)
seges_summary <- seges_selected %>%
  values %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
  )
topsoil_dif_rel <- seges_summary / dsc_summary
topsoil_dif_abs  <- seges_summary - dsc_summary
round(dsc_summary, 3)
round(seges_summary, 3)
round(topsoil_dif_rel, 3) - 1
round(topsoil_dif_abs, 3)

dsc_seges <- rbind(dsc_selected, seges_selected)

autoplot(dsc_seges, aes(col = db))

top_clay_plot <- ggplot(data = dsc_seges, aes(x = HUMUS, y = clay, col = db)) +
  geom_point() +
  geom_smooth() +
  xlab("SOM (%)") +
  ylab("Clay (%)")

tiff(
  paste0(dir_fig, "/top_clay_plot.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(top_clay_plot)
dev.off()

top_silt_plot <- ggplot(data = dsc_seges, aes(x = HUMUS, y = silt, col = db)) +
  geom_point() +
  geom_smooth() +
  xlab("SOM (%)") +
  ylab("Silt (%)")

tiff(
  paste0(dir_fig, "/top_silt_plot.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(top_silt_plot)
dev.off()

top_fines_plot <- ggplot(data = dsc_seges, aes(x = HUMUS, y = fine_sum, col = db)) +
  geom_point() +
  geom_smooth() +
  xlab("SOM (%)") +
  ylab("Clay [%] + Silt [%]")

tiff(
  paste0(dir_fig, "/top_fines_plot.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(top_fines_plot)
dev.off()

top_clayfines_plot <- ggplot(data = dsc_seges, aes(x = HUMUS, y = clay_fines, col = db)) +
  geom_point() +
  geom_smooth() +
  xlab("SOM (%)") +
  ylab("Clay [%] / (Clay [%] + Silt [%])")

tiff(
  paste0(dir_fig, "top_clayfines_plot.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(top_clayfines_plot)
dev.off()

mean(dsc_selected$clay_fines)
mean(seges_selected$clay_fines)

difference_dsc_seges <- mean(seges_selected$clay_fines) / mean(dsc_selected$clay_fines)

difference_dsc_seges

# Nearest points comparison

near_dsc <- nearest(dsc_selected, seges_selected)
near_seges <- nearest(seges_selected, dsc_selected)

pairs_1 <- bind_spat_cols(
  dsc_selected,
  seges_selected[near_dsc$to_id]
  )
pairs_2 <- bind_spat_cols(
  dsc_selected[near_seges$to_id],
  seges_selected
)

allpairs <- bind_spat_rows(pairs_1, pairs_2) %>%
  distinct()

allpairs %>% names()

allpairs %<>%
  mutate(
    clay_norem = clay...8,
    silt_norem = silt...9,
    finesum_norem = fine_sum...17,
    clayfines_norem = clay_fines...18,
    SOM_norem = HUMUS...16,
    clay_rem = clay...26,
    silt_rem = silt...27,
    finesum_rem = fine_sum...36,
    clayfines_rem = clay_fines...37,
    SOM_mean = (HUMUS...16 + HUMUS...35) / 2,
    clay_dif = clay_rem - clay_norem,
    clay_rel = clay_rem / clay_norem,
    silt_dif = silt_rem - silt_norem,
    finesum_dif = finesum_rem - finesum_norem,
    clayfines_dif = clayfines_rem - clayfines_norem,
    SOM_cat = cut(SOM_norem, c(0, 4, 10))
  ) %>%
  filter(
    clay_norem < 20
  )

ggplot(allpairs, aes(x = clay_norem, y = clay_rem, col = SOM_cat)) +
  geom_point() +
  geom_smooth(method=lm) +
  coord_equal() +
  geom_abline()

ggplot(allpairs, aes(x = silt_norem, y = silt_rem, col = SOM_cat)) +
  geom_point() +
  geom_smooth(method=lm) +
  coord_equal() +
  geom_abline()

ggplot(allpairs, aes(x = finesum_norem, y = finesum_rem, col = SOM_cat)) +
  geom_point() +
  geom_smooth(method=lm) +
  coord_equal() +
  geom_abline()

ggplot(allpairs, aes(x = clayfines_norem, y = clayfines_rem)) +
  geom_point() +
  coord_equal() +
  geom_abline()

ggplot(allpairs, aes(x = SOM_norem, y = clay_dif)) +
  geom_point() +
  geom_smooth()

mean(allpairs$clay_dif)
sd(allpairs$clay_dif)

mean(allpairs$clayfines_dif)
sd(allpairs$clayfines_dif)

mean(allpairs$clayfines_dif) / mean(allpairs$clayfines_norem)

plot(allpairs$clayfines_norem, allpairs$clayfines_dif)

plot(
  x = allpairs$clay_norem / allpairs$SOM_norem,
  y = allpairs$clayfines_dif,
)

allpairs %>%
  select(
    c(clay_norem, silt_norem, finesum_norem, clayfines_norem, clay_rem, silt_rem, finesum_rem, clayfines_rem, SOM_mean, clay_dif, silt_dif, finesum_dif, clayfines_dif, SOM_norem, clay_rel)
  ) %>%
  values() %>%
  cor()

claymodel <- lm(
  "clay_rem ~ clay_norem + silt_norem + SOM_norem",
  allpairs
  )
claymodel
summary(claymodel)

claymodel2 <- lm(
  "clay_rel ~ clay_norem + silt_norem + SOM_norem",
  allpairs
)
claymodel2
summary(claymodel2)

siltmodel <- lm(
  "silt_rem ~ clay_norem + silt_norem + SOM_norem",
  allpairs
)
siltmodel
summary(siltmodel)

diffmodel <- lm(
  "clayfines_dif ~ clay_norem + silt_norem + SOM_mean",
  allpairs
)
diffmodel
summary(diffmodel)

# Heatmaps
allpairs %>%
  values() %>%
  ggplot(aes(x = clay_norem, y = silt_norem, z = clay_rem)) +
  stat_summary_2d(binwidth = 1)

allpairs %>%
  values() %>%
  ggplot(aes(x = clay_norem, y = silt_norem, z = clay_rel)) +
  stat_summary_2d(binwidth = 1)

allpairs %>%
  values() %>%
  ggplot(aes(x = SOM_mean, y = silt_norem, z = silt_rem)) +
  stat_summary_2d(binwidth = 1)

# Test 3: Profiles vs DSC (More than 5% SOM)

t3_dsc_analysis <- dsc_analysis %>%
  filter(HUMUS > 5)

profiles_top <- profiles_texture %>%
  filter(upper < 25) %>%
  mutate(
    HUMUS = SOC / 0.587,
    fine_sum = clay + silt,
    clay_fines = clay / (fine_sum)
  ) %>%
  filter(
    HUMUS < 10,
    is.finite(clay_fines)
  )

profiles_top_somrem <- profiles_top %>%
  filter(HUMUS > 5)

distances_profdsc <- distance(profiles_top_somrem, t3_dsc_analysis)

t3_mindist_prof <- apply(distances_profdsc, 1, min)
t3_mindist_dsc <- apply(distances_profdsc, 2, function(x) min(x, na.rm = TRUE))

t3_dsc_selected <- t3_dsc_analysis %>%
  filter(t3_mindist_dsc < 1000)
t3_prof_selected <- profiles_top_somrem %>%
  filter(t3_mindist_prof < 1000)

t3_dsc_summary <- t3_dsc_selected %>%
  values %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
  )
t3_prof_summary <- t3_prof_selected %>%
  values %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
  )
t3_topsoil_dif_rel <- t3_prof_summary / t3_dsc_summary
t3_topsoil_dif_abs  <- t3_prof_summary - t3_dsc_summary
round(t3_dsc_summary, 3)
round(t3_prof_summary, 3)
round(t3_topsoil_dif_rel, 3) - 1
round(t3_topsoil_dif_abs, 3)

# Test 4: SEGES vs profiles (less than 5% SOM)

t4_seges_analysis <- seges_analysis %>%
  filter(HUMUS < 5)

profiles_top_norem <- profiles_top %>%
  filter(HUMUS < 5)

distances_profseges <- distance(profiles_top_norem, t4_seges_analysis)

t4_mindist_prof <- apply(distances_profseges, 1, function(x) min(x, na.rm = TRUE))
t4_mindist_seges <- apply(distances_profseges, 2, function(x) min(x, na.rm = TRUE))

t4_seges_selected <- t4_seges_analysis %>%
  filter(t4_mindist_seges < 1000)
t4_prof_selected <- profiles_top_norem %>%
  filter(t4_mindist_prof < 1000)

t4_seges_summary <- t4_seges_selected %>%
  values %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
  )
t4_prof_summary <- t4_prof_selected %>%
  values %>%
  reframe(
    clay_q = quantile(clay, c(0.25, 0.5, 0.75)),
    silt_q = quantile(silt, c(0.25, 0.5, 0.75)),
    fines_q = quantile(fine_sum, c(0.25, 0.5, 0.75)),
    clayfines_q = quantile(clay_fines, c(0.25, 0.5, 0.75)),
    q = c(0.25, 0.5, 0.75)
  )
t4_topsoil_dif_rel <- t4_seges_summary / t4_prof_summary
t4_topsoil_dif_abs  <- t4_seges_summary - t4_prof_summary
round(t4_prof_summary, 3)
round(t4_seges_summary, 3)
round(t4_topsoil_dif_rel, 3) - 1
round(t4_topsoil_dif_abs, 3)

# Effect of clay increase

source("f_classify_soil_JB.R")

class_asis <- classify_soil_JB(
  clay = dsc_analysis$clay,
  silt = dsc_analysis$silt,
  sand_f = dsc_analysis$fine_sand,
  SOM = dsc_analysis$HUMUS,
  CaCO3 = dsc_analysis$CaCO3
) %>% table()
 
class_plus_0_5 <- classify_soil_JB(
  clay = dsc_analysis$clay + 0.5,
  silt = dsc_analysis$silt,
  sand_f = dsc_analysis$fine_sand,
  SOM = dsc_analysis$HUMUS,
  CaCO3 = dsc_analysis$CaCO3
) %>% table()

class_plus_1_0 <- classify_soil_JB(
  clay = dsc_analysis$clay + 1,
  silt = dsc_analysis$silt,
  sand_f = dsc_analysis$fine_sand,
  SOM = dsc_analysis$HUMUS,
  CaCO3 = dsc_analysis$CaCO3
) %>% table()

class_plus_1_5 <- classify_soil_JB(
  clay = dsc_analysis$clay + 1.5,
  silt = dsc_analysis$silt,
  sand_f = dsc_analysis$fine_sand,
  SOM = dsc_analysis$HUMUS,
  CaCO3 = dsc_analysis$CaCO3
) %>% table()

JB_change_barplot <- bind_rows(
  class_asis,
  class_plus_0_5,
  class_plus_1_0,
  class_plus_1_5
) %>%
  mutate(
    clay_increase = as.factor(c(0, 0.5, 1, 1.5))
  ) %>%
  pivot_longer(
    cols = -clay_increase,
    values_to = "N",
    names_to = "JB_class"
  ) %>%
  mutate(
    JB_class = factor(JB_class, levels = as.character(1:10), labels = as.character(1:10))
  ) %>%
  ggplot(aes(x = JB_class, y = N, fill = clay_increase)) +
  geom_bar(stat = "identity", position = position_dodge())
  

tiff(
  paste0(dir_fig, "JB_change_barplot.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
print(JB_change_barplot)
dev.off()

class_nochange <- classify_soil_JB(
  clay = dsc_analysis$clay,
  silt = dsc_analysis$silt,
  sand_f = dsc_analysis$fine_sand,
  SOM = dsc_analysis$HUMUS,
  CaCO3 = dsc_analysis$CaCO3
)

classchange <- sapply(
  seq(0, 2, by = 0.01),
  function(x) {
    class_new <- classify_soil_JB(
      clay = dsc_analysis$clay + x,
      silt = dsc_analysis$silt,
      sand_f = dsc_analysis$fine_sand,
      SOM = dsc_analysis$HUMUS,
      CaCO3 = dsc_analysis$CaCO3
    )
    out <- mean(class_new != class_nochange)
    return(out)
  }
)

tiff(
  paste0(dir_fig, "changeplot_line.tiff"),
  width = 16,
  height = 10,
  units = "cm",
  res = 300
)
plot(
  x = seq(0, 2, by = 0.01),
  y = classchange*100,
  xlab = "Clay increase (%)",
  ylab = "Changed JB class (%)",
  type = "l"
)
dev.off()


# END