# Process and merge datasets

# To do :
# 1: Load datasets
# - Danish Soil Classification [DSC] (OK)
# - SEGES (OK)
# - SINKS (OK)
# - Profiles: Texture (not right now)
# - Profiles: Water retention (later)
# - Profiles: Drainage (later)
# - Lucas database?
# 2: Standardize
# - Uniform column names
# - Scale texture fractions to a sum of 100.
# 3: Merge
# 4: Extract covariates
# 5: Get folds
# 6: Write to table

# What fields should the final dataset contain?
# - General ID, "ID_new"
# - Database origin, "db"
# - Original ID "ID_old"
# - Date "date"
# - UTMX
# - UTMY
# - Upper boundary "upper"
# - Lower boundary "lower"
# - Clay
# - Silt
# - Fine sand
# - Coarse sand
# - SOC
# - CaCO3
# - Indicate removal of SOM for texture analyses
# - BD (later)
# - pH (later)
# - Nutrients: Pt, Kt, Mgt, CUT (later)
# - Water retention (later)
# - Texture sub-fractions: Gsilt, GfinSD (later)

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
  "SOM_removed"
)

# Start up

library(terra)
library(magrittr)
library(RODBC)
library(dplyr)
library(tidyr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/digijord_data/")


# 1.1: Load DSC observations

dsc <- dir_dat %>%
  paste0(
    .,
    "/observations/DanishSoilClassification/DLJ/DLJmDecimaler_DKBH.shp"
  ) %>%
  vect

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

# 1.4: Profiles: Texture (not right now)
# 1.5: Profiles: Water retention (later)
# 1.6: Profiles: Drainage (later)
# 1.7 Lucas database?

# 2: Data standardization
# 2.1: DSC

# Remove negative texture fraction measurements

dsc2 <- dsc %>%
  values() %>%
  mutate(
    db = "Danish Soil Classification",
    ID_old = provenr,
    date = paste0("19", Dato),
    upper = DybFra,
    lower = DybTil,
    clay = Ler * 100 / (Ler + Silt + FinSD + GrovSD),
    silt = Silt * 100 / (Ler + Silt + FinSD + GrovSD),
    fine_sand = FinSD * 100 / (Ler + Silt + FinSD + GrovSD),
    coarse_sand = GrovSD * 100 / (Ler + Silt + FinSD + GrovSD),
    logSOM = log(Humus),
    logCaCO3 = log(CaCO3)
  )

# END