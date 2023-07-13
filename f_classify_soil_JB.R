# Function to calculate JB code

classify_soil_JB <- function(clay, silt, sand_f, SOM, CaCO3, SOM_factor = 1) {
  out <- rep(0, length(clay))
  out[CaCO3 > 10] <- 12
  out[out == 0 & SOM * SOM_factor > 10] <- 11
  out[out == 0 & clay < 5 & silt < 20] <- 1
  out[out == 1 & sand_f > 50] <- 2
  out[out == 0 & clay < 5] <- 10
  out[out == 0 & clay < 10 & silt < 25] <- 3
  out[out == 3 & sand_f > 40] <- 4
  out[out == 0 & clay < 10] <- 10
  out[out == 0 & clay < 15 & silt < 30] <- 5
  out[out == 5 & sand_f > 40] <- 6
  out[out == 0 & clay < 15] <- 10
  out[out == 0 & clay < 25 & silt < 35] <- 7
  out[out == 0 & clay < 25] <- 10
  out[out == 0 & clay < 45 & silt < 45] <- 8
  out[out == 0 & clay < 45] <- 10
  out[out == 0 & silt < 50] <- 9
  out[out == 0 & clay < 50 & silt > 20] <- 10
  out[out == 0] <- NA
  return(out)
}

# END
