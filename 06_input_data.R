# 02: Process and merge datasets

# To do :
# 1: Load datasets
# - Danish Soil Classification [DSC] (OK)
# - SEGES (OK)
# - SINKS (OK)
# - Profiles: Texture (OK)
# - Profiles: Water retention (later)
# - Profiles: Drainage (later)
# - Lucas database?
# - Forest samples
# 2: Standardize
# - Uniform column names (OK)
# - Scale texture fractions to a sum of 100 (OK)

# What fields should the final dataset contain?
# - General ID, "ID_new" (OK)
# - Database origin, "db" (OK)
# - Original ID "ID_old" (OK)
# - Date "date" (OK)
# - UTMX (OK)
# - UTMY (OK)
# - Upper boundary "upper" (OK)
# - Lower boundary "lower" (OK)
# - Clay (OK)
# - Silt (OK)
# - Fine sand (OK)
# - Coarse sand (OK)
# - SOC (OK)
# - CaCO3 (OK)
# - Indicate removal of SOM for texture analyses (OK)
# - BD (later)
# - pH (OK)
# - Nutrients: Pt, Kt, Mgt, CUT, TOTALN TOTALP ORGAP    K   NA   CA   MG (later)
# - Water retention (later)
# - Texture sub-fractions: Gsilt, GfinSD, GRVSILT S63 S125 S200 S500 (later)

mycolnames <- c(
  "ID_new",
  "db",
  "ID_old",
  "date",
  "UTMX",
  "UTMY",
  "upper",
  "lower",
  "clay",
  "silt",
  "fine_sand",
  "coarse_sand",
  "SOC",
  "CaCO3",
  "SOM_removed",
  "pH",
  "N",
  "BD",
  "CEC",
  "imputed"
)

# Start up

library(terra)
library(magrittr)
library(RODBC)
library(dplyr)
library(tidyr)
library(stringi)
library(readxl)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")


# 1.1: Load DSC observations

dsc <- dir_dat %>%
  paste0(
    .,
    "/observations/DanishSoilClassification/DLJ/DLJmDecimaler_DKBH.shp"
  ) %>%
  vect()

# 1.2: Load SEGES samples

SEGES <- dir_dat %>%
  paste0(
    .,
    "/observations/SEGES_samples/SEGES_samples_cleaned.csv"
  ) %>%
  read.table(
    sep = ";",
    header = TRUE
  )

# 1.3: Load SINKS data
# In the long run I should use all the data.
# However, I will use the topsoil until the dataset has been fixed.
# Use pivot_longer to get a row for each sample.

SINKS_db <- dir_dat %>%
  paste0(., "/observations/SINKS/SinksDatabase_V2.mdb")

con2 <- odbcConnectAccess2007(SINKS_db)

SINKS <- sqlFetch(con2, "bMHG_RODALTDETDUVIL_01JAN2011")

# 1.4: Profiles:

# 1.4.1: Profiles coordinates

profiles_shp <- dir_dat %>%
  paste0(
    .,
    "/observations/profiles/Profiles_coordinates_new/Profiles_coordinates_new.shp"
  ) %>%
  vect()

# 1.4.2: Profiles: Texture

profiles_db <- dir_dat %>%
  paste0(., "/observations/profiles/DDJD2023.accdb")

con3 <- odbcConnectAccess2007(profiles_db)

profiles_texture <- sqlFetch(con3, "ANALYSE")

# 1.4.3: Profiles: Horizons

profiles_horizons <- dir_dat %>%
  paste0(
    .,
    "/observations/profiles/DDJD_horizons/DDJD_HORISONT_FRA_TIL.xlsx"
  ) %>%
  read_excel()

# 1.5: Profiles: Water retention (later)

profiles_retention <- sqlFetch(con3, "VANDRETENTION")

# 1.6: Profiles: Drainage (later)
# 1.7: Lucas database?

# 1.8: Forest samples

forests_tax <- dir_dat %>%
  paste0(
    .,
    "/observations/forest_samples/Taksation_ETRS89.xlsx"
  ) %>%
  read_excel(.name_repair = "universal")

forests_dsc <- dir_dat %>%
  paste0(
    .,
    "/observations/forest_samples/Skovjord_ETRS89.xlsx"
  ) %>%
  read_excel(.name_repair = "universal")

# 2: Data standardization
# 2.1: DSC

dsc_coords <- geom(dsc) %>%
  as.data.frame() %>%
  select(c(x, y))

dsc_oldcoords <- dsc %>%
  values() %>%
  select(c(UTMX, UTMY))

dsc2 <- dsc %>%
  values() %>%
  mutate(
    across(
      where(is.numeric),
      function(x) replace(x, x == -1, NA)
    ) # Remove negative texture fraction measurements
  ) %>%
  bind_cols(dsc_coords) %>%
  select(-c(UTMX, UTMY)) %>%
  mutate(
    db = "Danish Soil Classification",
    UTMX = x,
    UTMY = y,
    ID_old = provenr,
    date = paste0("19", Dato),
    upper = DybFra,
    lower = DybTil,
    clay = Ler * 100 / (Ler + Silt + FinSD + GrovSD),
    silt = Silt * 100 / (Ler + Silt + FinSD + GrovSD),
    fine_sand = FinSD * 100 / (Ler + Silt + FinSD + GrovSD),
    coarse_sand = GrovSD * 100 / (Ler + Silt + FinSD + GrovSD),
    SOC = Humus * 0.587,
    SOM_removed = 0,
    imputed = FALSE
  ) %>%
  select(any_of(mycolnames))

# 2.2 SEGES

seges2 <- SEGES %>%
  mutate(
    db = "SEGES",
    date = as.character(date),
    ID_old = as.character(ID),
    UTMX = UTME,
    UTMY = UTMN,
    upper = 0,
    lower = 25,
    tsum = LerPct + SiltPct + FinsandPct + GrovsandPct,
    clay = case_when(
      is.na(tsum) ~ LerPct,
      !is.na(tsum) ~ LerPct * 100 / tsum
    ),
    silt = case_when(
      is.na(tsum) ~ SiltPct,
      !is.na(tsum) ~ SiltPct * 100 / tsum
    ),
    fine_sand = case_when(
      is.na(tsum) ~ FinsandPct,
      !is.na(tsum) ~ FinsandPct * 100 / tsum
    ),
    coarse_sand = case_when(
      is.na(tsum) ~ GrovsandPct,
      !is.na(tsum) ~ GrovsandPct * 100 / tsum
    ),
    SOM_removed = 1,
    pH = Rt - 0.5,
    N = TotalNPct,
    imputed = FALSE
  ) %>%
  select(any_of(mycolnames))

# 2.3 SINKS

depths_A <- c("A1", "A2", "A3", "A4")
depths_S <- c("S1", "S2", "S3", "S4")
depths_D <- c("D1", "D2", "D3", "D4")

new_names <- stri_replace_all_fixed(
  names(SINKS),
  c(depths_A, depths_S),
  c(depths_D, depths_D),
  vectorize_all = FALSE
)

new_names <- stri_replace_all_fixed(
  new_names,
  "ID_DJF",
  "IDDJF",
  vectorize_all = FALSE
)

new_names <- stri_replace_all_fixed(new_names, "__", "_", vectorize_all = FALSE)

SINKS_newnames <- SINKS
names(SINKS_newnames) <- new_names

DX <- paste0(depths_D, "_")

D_cols <- depths_D %>%
  paste(collapse = "|") %>%
  grep(
    names(SINKS_newnames),
    value = TRUE
  ) %>%
  unique() %>%
  setdiff(DX)

newcols <- D_cols %>%
  strsplit("_") %>%
  unlist() %>%
  unique() %>%
  setdiff(depths_D)

sinks2 <- SINKS_newnames %>%
  select(-any_of(DX)) %>%
  pivot_longer(
    all_of(D_cols),
    names_to = c("depth", ".value"),
    names_sep = "_"
  )

rownas <- apply(sinks2, 1, function(x) sum(is.na(x)))

lapply(sinks2, function(x) sum(!is.na(x))) %>% unlist()
lapply(sinks2, function(x) length(unique(x))) %>% unlist()

sinks2 %>%
  filter(!is.na(KunPro)) %>%
  select(depth, PH, TOC1, TOC2, TOC3, NTotal)
# Field "KunPro" is the closest thing to a unique ID for each depth sample

sinks3 <- sinks2 %>%
  filter(!is.na(KunPro))

sinks3 %>%
  group_by(depth) %>%
  summarise(
    mean_TOC1 = mean(TOC1, na.rm = TRUE),
    mean_TOC2 = mean(TOC2, na.rm = TRUE),
    mean_TOC3 = mean(TOC3, na.rm = TRUE)
  )

sinks3 %>% filter(!is.na(TOC2) & !is.na(TOC3))
# ^ None

sinks3 %>%
  filter(!is.na(TOC1) & !is.na(TOC3)) %>%
  select(depth, TOC1, TOC2, TOC3, NTotal) %>%
  print(n = 100)
# ^ All depth one
# C/N OK, but one C value may be N

sinks3 %>%
  filter(!is.na(TOC1) & !is.na(TOC2)) %>%
  select(depth, TOC1, TOC2, TOC3, NTotal) %>%
  print(n = 100)
# ^ Never depth one
# C/N OK, but one C value may be N

sinks3 %>%
  filter(TOC1 != TOC3) %>%
  select(depth, TOC1, TOC2, TOC3, NTotal) %>%
  print(n = 100)
# Only two values, TOC3 is N

sinks3 %>%
  filter(TOC1 != TOC2) %>%
  select(depth, TOC1, TOC2, TOC3, NTotal) %>%
  print(n = 100)
# Only one value, TOC2 is N

sinks3 %>%
  filter(
    TOC1 < NTotal | TOC2 < NTotal | TOC3 < NTotal
  ) %>%
  select(depth, TOC1, TOC2, TOC3, NTotal) %>%
  print(n = 100)

n_unique <- sinks3 %>%
  select(TOC1, TOC2, TOC3, NTotal) %>%
  apply(1, function(x) length(unique(x)))
# All points have three unique values

maxvalue <- sinks3 %>%
  select(TOC1, TOC2, TOC3, NTotal) %>%
  apply(1, function(x) max(x, na.rm = TRUE))

minvalue <- sinks3 %>%
  select(TOC1, TOC2, TOC3, NTotal) %>%
  apply(1, function(x) min(x, na.rm = TRUE))

plot(log(minvalue), log(maxvalue))

plot(maxvalue / minvalue)

# I assume that maxvalue is always SOC, and minvalue is always N
sinks3$SOC <- maxvalue
sinks3$N <- minvalue

# Impute values for missing depth intervals
sinks2_filled <- sinks2 

maxvalue2 <- pmax(
  sinks2_filled$TOC1,
  sinks2_filled$TOC2,
  sinks2_filled$TOC3,
  sinks2_filled$NTotal,
  na.rm = TRUE
)

minvalue2 <- pmin(
  sinks2_filled$TOC1,
  sinks2_filled$TOC2,
  sinks2_filled$TOC3,
  sinks2_filled$NTotal,
  na.rm = TRUE
)

sinks2_filled$SOC <- maxvalue2
sinks2_filled$N <- minvalue2

sinks2_filled %<>%
  mutate(
    PH = case_when(
      PH == 0 ~ NA,
      .default = PH
    )
  ) %>%
  group_by(IDBID) %>%
  fill(
    any_of(c("SOC", "N", "PH")),
    .direction = "down"
  ) %>%
  ungroup() %>%
  mutate(
    imputed = case_when(
      is.na(KunPro) ~ TRUE,
      .default = FALSE
    ),
    SOC = case_when(
      SOC > 12 & imputed == TRUE ~ 12,
      .default = SOC
    )
  )

# end of imputation

sinks3 <- sinks2_filled

sinks3$date <- sinks3$SD_Dato %>%
  as.character() %>%
  stri_replace_all_fixed("-", "") %>%
  substr(., 1, 8)

uppers <- c(0:3) * 34
lowers <- uppers + 30

sinks_upperlower <- data.frame(
  depth = depths_D,
  upper = uppers,
  lower = lowers
)

sinks4 <- sinks3 %>%
  mutate(
    db = "SINKS",
    ID_old = paste0(IDBID, "_", depth),
    UTMX = SD_X,
    UTMY = SD_Y,
    pH = PH
  ) %>%
  right_join(sinks_upperlower) %>%
  select(any_of(mycolnames))

# 2.4: Profiles texture

profiles_texture %>% filter(PRFRA > PRTIL)
profiles_texture %>% filter(TEKTURSUM < 0)

# Replace negative values (ok)
# Fix PRFRA PRTIL (ok)
# (Remove samples where upper and lower are still missing)

# Remove false zeroes
# - texture
# - BD (OK)
# - pH (OK)
# - SOC
# - CaCO3
# - CEC (OK)
# Standardize texture sums

tex_negs <- profiles_texture %>%
  select(-c(PRFRA, PRTIL)) %>%
  select_if(is.numeric) %>%
  as.matrix() %>%
  apply(., 1, function(x) {
    out <- sum(x < 0, na.rm = TRUE)
    return(out)
  })

profiles_texture[tex_negs > 0, ]
profiles_texture %>% filter((PRTIL - PRFRA) > 100)

check_negs <- profiles_texture %>%
  select(-c(PRFRA, PRTIL)) %>%
  select_if(is.numeric) %>%
  colnames()

profiles_tex2 <- profiles_texture %>%
  mutate(
    across(
      any_of(check_negs),
      function(x) replace(x, x < 0, NA)
    ) # Remove negative measurements
  ) %>%
  mutate(
    PRFRA = case_when(
      PRFRA > PRTIL & PROFILNR == 2517 ~ PRFRA * -1,
      .default = PRFRA
    ),
    PRTIL = case_when(
      PRFRA > PRTIL & PROFILNR == 2531 ~ 170,
      (PRTIL - PRFRA) > 100 & PROFILNR == 2539 ~ 30,
      PRFRA > PRTIL & PROFILNR %in% c(2528, 2529, 2530) ~ PRFRA,
      .default = PRTIL
    )
  )

profiles_tex2 %>%
  filter(PROFILNR %in% c(2517, 2528, 2529, 2530, 2531, 2539)) %>%
  select(c(PROFILNR, HORISONTNR, HORISONTPR, PRFRA, PRTIL))

plot(profiles_tex2$PRFRA, profiles_tex2$PRTIL)

# Analyse SOM removal
# SOM only removed for soils with more than 5% SOM

# profiles_tex3 <- profiles_tex2 %>%
#   mutate(claysilt = LER/(LER + SILT))
#
# plot(
#   profiles_tex3$HUMUS,
#   profiles_tex3$claysilt,
#   xlim = c(0, 10),
#   abline(v = 5, col = "blue")
#   )

# Remove false zeroes

plot(profiles_tex2$LER) # ok
plot(sqrt(profiles_tex2$HUMUS)) # ok
plot(sqrt(profiles_tex2$CACO3)) # ok
plot(profiles_tex2$PHCACL2) # ok
plot(profiles_tex2$CEC) # ok
plot(profiles_tex2$VOLVGT) # ok

profiles_tex2 %<>%
  mutate(
    PHCACL2 = case_when(
      PHCACL2 < 1 ~ NA,
      .default = PHCACL2
    ),
    CEC = case_when(
      CEC == 0 ~ NA,
      .default = CEC
    ),
    VOLVGT = case_when(
      VOLVGT == 0 ~ NA,
      .default = VOLVGT
    )
  )

profiles_tex2 %<>%
  rowwise() %>%
  mutate(
    tminsum = sum(LER, SILT, GRVSILT, S63, S125, S200, S500, na.rm = TRUE)
  ) %>%
  ungroup()

plot(profiles_tex2$tminsum)

profiles_tex2 %>%
  as.data.frame() %>%
  filter(TEKSTURART == "Kalk")

# TEKSTURART analysis:
# "Humus"         # More than 400 samples, all texture fractions are false zeroes
# "tekstur"       # 996 samples, false zeroes where all texture fractions are 0
# "ingen analyse" # One sample, texture fractions are false zeroes,
#                 # CaCO3 is a "true" 0 (i.e. not measured)
# "Ingen analyse" # 198 samples, all zero rows are false zeroes
# "Kalk"          # 27 samples, texture fractions and humus are false zeroes
# "Humus/kalk"    # More than 100 samples, all texture fractions are false zeroes
# NA              # More than 1000 samples, only false zeroes where all measurements are 0/NA
# "0,0"           # Only two, no analysis, all zeroes are false

# "humus"         # 11 samples, no false zeroes
# "Tekstur"       # More than 9000 samples, no false zeroes
# "kalk"          # 20 samples, no false zeroes, but "TOTAL KULSTOF" is unreliable
# "tekstur/kalk"  # Three samples, no false zeroes
# " "             # 11 samples, no false zeroes

# Summary:
# If texture, humus and CaCO3 are all zero, they are false zeroes
# If texture fractions are all zero, they are false zeroes
# If CaCO3 is >0, and texture fractions + humus are zero, humus is a false zero
# If one texture fraction is >0, humus and CaCO3 may be true zeroes
# If one texture fraction is >0, the other fractions may be true zeroes
# If humus or texture fractions are >0, CaCO3 may be true zeroes

profiles_tex2 %<>%
  rowwise() %>%
  mutate(
    TEKTURSUM_NY = sum(LER, SILT, GRVSILT, S63, S125, S200, S500, HUMUS, CACO3,
      na.rm = TRUE
    ),
    TEKSTUR_NZERO = sum(c(LER, SILT, GRVSILT, S63, S125, S200, S500) == 0,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

profiles_tex2 %>%
  as.data.frame() %>%
  filter(tminsum > 101) %>%
  arrange(-tminsum)
profiles_tex2 %>%
  as.data.frame() %>%
  filter(TEKTURSUM > 101) %>%
  arrange(-TEKTURSUM)
profiles_tex2 %>%
  as.data.frame() %>%
  filter(TEKTURSUM_NY < 99 & tminsum > 0) %>%
  arrange(-TEKTURSUM_NY)

profiles_tex2 %>%
  filter(TEKSTUR_NZERO == 2) %>%
  as.data.frame()
profiles_tex2 %>% filter(TEKSTUR_NZERO == 1 & TEKTURSUM_NY < 50)

# Mineral texture sums above 101 are generally due to repeated values for one
# of the sand fractions.
# For one sample (profile 975, horizon 4), S500 seems to be the sum of the other
# texture fractions.
# If two or more texture fractions are zero, they are false zeroes

# Check if any profiles consistently lack SOM measurements (inconclusive)

# profiles_zerosom <- profiles_tex2 %>%
#   group_by(PROFILNR) %>%
#   summarise(humsum = sum(HUMUS, na.rm = TRUE)) %>%
#   filter(humsum == 0)
#
# profiles_tex2 %>%
#   filter(
#   PROFILNR %in% profiles_zerosom$PROFILNR & TEKSTURART != "Ingen analyse")

# Decisions:
# 1: Remove repeated texture values where the sum is higher than 101% (remove
# last).
# 2: If two or more texture fractions are zero, set them to NA. (OK)
# 3: For profile 975, horizon 4, recalculate S500. (OK)
# 4: If all fractions, including humus and CaCO3 are zero/NA, set them all to
# NA. (OK)
# 5: Set humus zeroes to NA if CaCO3 > 0 and texture is zero/NA. (OK)

t_colnames <- c("LER", "SILT", "GRVSILT", "S63", "S125", "S200", "S500")

profiles_tex4 <- profiles_tex2 %>%
  rowwise() %>%
  mutate(
    across(
      c(SILT, GRVSILT, S63, S125, S200, S500),
      ~ case_when(
        TEKTURSUM_NY > 101 & c_across(LER:S500)[match(cur_column(), t_colnames) - 1] == . ~ NA,
        .default = .
      )
    )
  ) %>%
  mutate(
    across(
      c(LER, SILT, GRVSILT, S63, S125, S200, S500),
      ~ case_when(
        TEKSTUR_NZERO > 1 & . == 0 ~ NA,
        .default = .
      )
    )
  ) %>%
  mutate(
    S500 = case_when(
      PROFILNR == 975 & HORISONTNR == 4 ~ 100 - sum(
        LER, SILT, GRVSILT, S63, S125, S200, HUMUS, CACO3,
        na.rm = TRUE
      ),
      .default = S500
    )
  ) %>%
  ungroup() %>%
  mutate(
    across(
      c(LER, SILT, GRVSILT, S63, S125, S200, S500, HUMUS, CACO3),
      ~ case_when(
        TEKTURSUM_NY == 0 ~ NA,
        .default = .
      )
    )
  ) %>%
  mutate(
    HUMUS = case_when(
      HUMUS == 0 & CACO3 > 0 & tminsum == 0 ~ NA,
      HUMUS == 0 & CACO3 > 0 & is.na(tminsum) ~ NA,
      .default = HUMUS
    )
  )

profiles_tex4 %<>%
  rowwise() %>%
  mutate(
    TEKTURSUM_NY = sum(LER, SILT, GRVSILT, S63, S125, S200, S500, HUMUS, CACO3,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

plot(profiles_tex4$TEKTURSUM_NY)

profiles_tex4 %>%
  filter(TEKTURSUM_NY > 101) %>%
  as.data.frame()

# profile 3280 is ok (sum only 102%)
# profile 3180 horizon 6 sample 2 seems to be a duplicate of the topsoil sample
# with S500 as the only exception. The column "TOTAL KULSTOF" seems to contain
# the only valid measurement from this sample. The texture fractions for this
# sample should therefore be set to NA.

plot(profiles_tex4$HUMUS, profiles_tex4$`TOTAL KULSTOF`)
cor(profiles_tex4$HUMUS, profiles_tex4$`TOTAL KULSTOF`, use = "pairwise.complete.obs")
plot(log(profiles_tex4$HUMUS), log(profiles_tex4$`TOTAL KULSTOF`))
cor(log(profiles_tex4$HUMUS), log(profiles_tex4$`TOTAL KULSTOF`), use = "pairwise.complete.obs", method = "spearman")

is_naturalnumber <- function(x, tol = .Machine$double.eps^0.5) {
  out <- x > tol & abs(x - round(x)) < tol
  return(out)
}

profiles_tex4 %>%
  filter(!(JBNR %in% 1:12)) %>%
  filter(!is.na(JBNR)) %>%
  filter(is_naturalnumber(`TOTAL KULSTOF`)) %>%
  as.data.frame()

plot(profiles_tex4$PHH2O) # simply remove zeroes.
plot(profiles_tex4$PHKCL) # simply remove zeroes.
plot(profiles_tex4$POROES) # simply remove zeroes.
plot(profiles_tex4$TOTALP) # simply remove zeroes.
plot(profiles_tex4$TOTALN) # simply remove zeroes.
plot(profiles_tex4$ORGAP) # simply remove zeroes.
plot(profiles_tex4$`UORGANISK P`) # simply remove zeroes.
plot(profiles_tex4$`CITRATOPL P`) # simply remove zeroes.
plot(profiles_tex4$`TOTAL KULSTOF`) # simply remove zeroes.

library(ggplot2)

profiles_tex4 %>%
  filter(JBNR %in% 1:12) %>%
  filter(HUMUS != 0) %>%
  filter(`TOTAL KULSTOF` != 0) %>%
  ggplot(aes(x = HUMUS, y = `TOTAL KULSTOF`)) +
  geom_point()

profiles_tex4 %>%
  filter(JBNR %in% 1:12) %>%
  filter(HUMUS != 0) %>%
  filter(`TOTAL KULSTOF` != 0) %>%
  filter(`TOTAL KULSTOF` > HUMUS) %>%
  as.data.frame()

# The columns "UORGANISK P", "CITRATOPL P", "TOTAL KULSTOF", "JBNR" seem to
# share many false zeroes. These should be set to NA.
# The columns "TOTALN", "TOTALP", "ORGAP", "UORGANISK P" share many false
# zeroes. These should be set to NA.
# In the columns K, NA, CA, MG, BASER, SURION, CEC, BASEMAETN the zeroes are
# false if CEC is either NA or zero.
# If SURION and BASEMAETN are both zero, they are false zeroes, irrespective of
# CEC.
# JBNR should only ever be an integer between 1 and 12.
# `TOTAL KULSTOF` values that are integers seem to contain the JBNR in cases
# Where the latter is missing.

# Decisions:
# 6: Remove all invalid JB numbers (not in 1:12). (ok)
# 7: Where JBNR is missing and total C is an integer, use total C for JBNR and
# delete the value from the total C column. Remove any invalid JBNR values
# again. (ok)
# 8: For profile 3180 horizon 6 sample 2 set all texture fractions to NA. (ok)
# 9: Set all zeroes to NA for PHH2O, PHKCL, POROES, TOTALP, TOTALN, ORGAP,
# UORGANISK P, CITRATOPL P, TOTAL KULSTOF. (ok)
# 10: If CEC is NA or zero, set 0 to NA for the columns K, NA, CA, MG, BASER,
# SURION, CEC, BASEMAETN. (ok)
# 11: If SURION and BASEMAETN are both zero, set them to NA. (ok)
# 12: For profile 2532 horizon 3, HUMUS is 100 due to an error caused by
# missing texture fractions (not measured). According to the field report (MH
# Greve & P Sørensen, 1990: Lokalitetskortlægning på marginaljord - et pilot-
# projekt, 2. revideret udgave, bilag), the HUMUS content for this horizon is
# 2.66%. Correct this manually. (ok)

profiles_tex4 %<>%
  mutate(
    JBNR = case_when(
      !(JBNR %in% 1:12) ~ NA,
      .default = JBNR
    ),
    JBNR = case_when(
      is.na(JBNR) & is_naturalnumber(`TOTAL KULSTOF`) ~ `TOTAL KULSTOF`,
      .default = JBNR
    ),
    `TOTAL KULSTOF` = case_when(
      is_naturalnumber(`TOTAL KULSTOF`) & `TOTAL KULSTOF` == JBNR ~ NA,
      .default = `TOTAL KULSTOF`
    ),
    JBNR = case_when(
      !(JBNR %in% 1:12) ~ NA,
      .default = JBNR
    )
  ) %>%
  mutate(
    across(
      c(LER, SILT, GRVSILT, S63, S125, S200, S500, HUMUS, CACO3),
      ~ case_when(
        PROFILNR == 3180 & HORISONTNR == 6 & HORISONTPR == 2 ~ NA,
        .default = .
      )
    ),
    across(
      c(
        PHH2O, PHKCL, POROES, TOTALP, TOTALN, ORGAP, `UORGANISK P`,
        `CITRATOPL P`, `TOTAL KULSTOF`
      ),
      ~ case_when(
        . == 0 ~ NA,
        .default = .
      )
    ),
    across(
      c(K, `NA`, CA, MG, BASER, SURION, BASEMAETN),
      ~ case_when(
        is.na(CEC) & . == 0 ~ NA,
        .default = .
      )
    )
  ) %>%
  mutate(
    SURION = case_when(
      PROFILNR == 3181 & SURION == 0 ~ NA,
      .default = SURION
    ),
    BASEMAETN = case_when(
      PROFILNR == 3181 & BASEMAETN == 0 ~ NA,
      .default = BASEMAETN
    ),
    HUMUS = case_when(
      PROFILNR == 2532 & HORISONTNR == 3 ~ 2.66,
      .default = HUMUS
    )
  )

profiles_tex4 %<>%
  rowwise() %>%
  mutate(
    TEKTURSUM = sum(LER, SILT, GRVSILT, S63, S125, S200, S500, HUMUS, CACO3,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    TEKTURSUM = case_when(
      TEKTURSUM == 0 ~ NA,
      .default = TEKTURSUM,
    )
  )

# Export corrected table

profiles_analyse_corrected <- profiles_tex4 %>%
  select(colnames(profiles_texture))

write.table(
  profiles_analyse_corrected,
  paste0(dir_dat, "/observations/profiles/ANALYSE_corrected_20230606.csv"),
  sep = ";",
  na = "",
  row.names = FALSE
)

# Impute missing horizon boundaries
# Get boundaries from horizons table
# Get upper and lower from horizons table when they are missing for the sample
# (ok)

profiles_horizons_small <- profiles_horizons %>%
  mutate(
    HOR_ID = paste0(PROFILNR, "_", HORISONTNR)
  ) %>%
  select(-c(PROFILNR, HORISONTNR, HORISONT))

profiles_tex5 <- profiles_analyse_corrected %>%
  mutate(
    PRFRA = case_when(
      is.na(PRFRA) ~ PRTIL,
      .default = PRFRA
    ),
    PRTIL = case_when(
      is.na(PRTIL) ~ PRFRA,
      .default = PRTIL
    )
  ) %>%
  mutate(
    HOR_ID = paste0(PROFILNR, "_", HORISONTNR)
  ) %>%
  left_join(
    profiles_horizons_small,
    "HOR_ID"
  ) %>%
  mutate(
    PRFRA = case_when(
      is.na(PRFRA) ~ FRA,
      .default = PRFRA
    ),
    PRTIL = case_when(
      is.na(PRTIL) ~ TIL,
      .default = PRTIL
    )
  ) %>%
  select(-c(HOR_ID, FRA, TIL))

# Standardization
# Only standardize texture fractions if the combined sum, including humus and
# CaCO3 are more than 90%.

profiles_tex5 %>%
  filter(is.na(HUMUS) & !is.na(LER)) %>%
  as.data.frame()

profiles_tex5 %>%
  filter(!is.na(LER)) %>%
  filter(TEKTURSUM < 99) %>%
  filter(is.na(HUMUS)) %>%
  filter(is.na(CACO3)) %>%
  slice_sample(n = 20) %>%
  as.data.frame()

# If TEKTURSUM is more than 90, use sum of mineral fractions
# If not, and if HUMUS or CaCO3 are available, use 100 - (HUMUS + CaCO3)
# Otherwise, do not standardize

profiles_tex5 %<>%
  mutate(
    db = "Profile database",
    ID_old = paste0(PROFILNR, "_", HORISONTNR, "_", HORISONTPR),
    upper = PRFRA,
    lower = PRTIL,
    SOC = case_when(
      !is.na(HUMUS) ~ HUMUS * 0.587,
      is.na(HUMUS) & is.na(CACO3) ~ `TOTAL KULSTOF`,
      is.na(HUMUS) & CACO3 == 0 ~ `TOTAL KULSTOF`,
      .default = NA
    ),
    SOM_removed = case_when(
      HUMUS >= 5 ~ 1,
      HUMUS < 5 ~ 0,
      .default = NA
    ),
    CaCO3 = CACO3,
    BD = VOLVGT,
    pH = PHCACL2,
    N = TOTALN,
    imputed = FALSE
  ) %>%
  rowwise() %>%
  mutate(
    tminsum = sum(LER, SILT, GRVSILT, S63, S125, S200, S500, na.rm = TRUE),
    not_SOM_CaCO3 = case_when(
      !is.na(sum(HUMUS, CACO3)) ~ 100 - (HUMUS + CACO3),
      !is.na(HUMUS) & is.na(CACO3) ~ 100 - HUMUS,
      is.na(HUMUS) & !is.na(CACO3) ~ 100 - CACO3,
      .default = NA
    ),
    clay = case_when(
      TEKTURSUM >= 90 ~ LER * 100 / (tminsum),
      TEKTURSUM < 90 & !is.na(not_SOM_CaCO3) ~ LER * 100 / (not_SOM_CaCO3),
      .default = LER
    ),
    silt = case_when(
      TEKTURSUM >= 90 ~ SILT * 100 / (tminsum),
      TEKTURSUM < 90 & !is.na(not_SOM_CaCO3) ~ SILT * 100 / (not_SOM_CaCO3),
      .default = SILT
    ),
    fsand_raw = case_when(
      sum(!is.na(c_across(GRVSILT:S125))) > 0 ~ sum(GRVSILT, S63, S125, na.rm = TRUE),
      .default = NA
    ),
    fine_sand = case_when(
      TEKTURSUM >= 90 ~ fsand_raw * 100 / (tminsum),
      TEKTURSUM < 90 & !is.na(not_SOM_CaCO3) ~ fsand_raw * 100 / (not_SOM_CaCO3),
      .default = fsand_raw
    ),
    csand_raw = case_when(
      sum(!is.na(c_across(S200:S500))) > 0 ~ sum(S200, S500, na.rm = TRUE),
      .default = NA
    ),
    coarse_sand = case_when(
      TEKTURSUM >= 90 ~ csand_raw * 100 / (tminsum),
      TEKTURSUM < 90 & !is.na(not_SOM_CaCO3) ~ csand_raw * 100 / (not_SOM_CaCO3),
      .default = csand_raw
    )
  ) %>%
  ungroup() %>%
  select(any_of(c(mycolnames, "PROFILNR")))


# Get coordinates and date

profiles_xy_date <- profiles_shp %>%
  values() %>%
  mutate(
    date = paste0(
      AARSTAL,
      formatC(MND, 2, 2, flag = "0"),
      formatC(DAG, 2, 2, flag = "0")
    ),
    UTMX = x,
    UTMY = y
  ) %>%
  select(-c(AARSTAL, MND, DAG, x, y)) %>%
  arrange(PROFILNR)

profiles_tex5 %<>%
  left_join(profiles_xy_date) %>%
  as.data.frame()


# 2.8: Forest samples

forests_tax2 <- forests_tax %>%
  mutate(
    db = "Forest evaluation",
    ID_old = as.character(key),
    date = "1980",
    upper = 0,
    lower = 20,
  ) %>%
  filter(
    !is.na(Humus) & !is.na(UTMN) & !is.na(UTMØ)
  )

forests_dsc2 <- forests_dsc %>%
  mutate(
    across(
      where(is.numeric),
      function(x) replace(x, x == -1, NA)
    ) # Remove negative values
  ) %>%
  arrange(Lokalitet) %>%
  group_by(Lokalitet) %>%
  arrange(Dybde.fra) %>%
  mutate(
    db = "Forests DSC",
    ID_old = paste0(Lokalitet, "_", row_number()),
    date = "1991",
    upper = Dybde.fra,
    lower = Dybde.til,
    CaCO3 = case_when(
      is.na(CaCO3) ~ 0,
      .default = CaCO3
    ),
    Finsand = Grovsilt + ..63.mym + ..125.mym,
    Grovsand = ..200.mym + ..500.mym
  ) %>%
  ungroup() %>%
  arrange(Lokalitet)

forests_all <- bind_rows(forests_tax2, forests_dsc2) %>%
  mutate(
    UTMX = UTMØ,
    UTMY = UTMN,
    SOM_removed = 0,
    SOC = Humus * 0.586
  ) %>%
  rowwise() %>%
  mutate(
    tsum = sum(Ler, Silt, Finsand, Grovsand, na.rm = TRUE),
    tsum_all = sum(Ler, Silt, Finsand, Grovsand, Humus, CaCO3, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    tsum = case_when(
      tsum == 0 ~ NA,
      .default = tsum
    ),
    clay = Ler * 100 / tsum,
    silt = Silt * 100 / tsum,
    fine_sand = Finsand * 100 / tsum,
    coarse_sand = Grovsand * 100 / tsum,
    imputed = FALSE
  ) %>%
  select(any_of(mycolnames))

# Write files

dir_obs_proc <- dir_dat %>%
  paste0(., "/observations/processed/") %T>%
  dir.create()

write.table(
  dsc2,
  paste0(dir_obs_proc, "dsc.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  seges2,
  paste0(dir_obs_proc, "seges.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  sinks4,
  paste0(dir_obs_proc, "SINKS.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  profiles_tex5,
  paste0(dir_obs_proc, "profiles_texture.csv"),
  row.names = FALSE,
  sep = ";"
)

write.table(
  forests_all,
  paste0(dir_obs_proc, "forest_samples.csv"),
  row.names = FALSE,
  sep = ";"
)

# END
