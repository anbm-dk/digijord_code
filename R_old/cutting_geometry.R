# Cutting geometry
# Anders Bjørn Møller

library(raster)
library(tidyr)
library(magrittr)
library(snow)

getwd()
# [1] "D:/anbm/digijord"

root <- getwd()

dem <- 'D:/anbm/digijord/cov_basic/DHM2015_terraen_10m.tif' %>% raster()

f <- function(x)
{
  out <- x*0 + 1
  return(out)
}

beginCluster(12)

clusterR(
  dem
  , calc
  , args = list(fun = f)
  , filename = paste0(root, '/Coast/dem_coast.tif')
  , datatype = 'INT2U'
  , overwrite = TRUE
  , format = 'GTiff'
)

endCluster()

# END
