# Function to classify SOC contents


classify_SOC <- function(SOC) {
  out <- rep(0, length(SOC))
  out[SOC > 12] <- 60
  out[out == 0 & SOC > 6] <- 12
  out[out == 0 & SOC > 3] <- 6
  out[out == 0 & !is.na(SOC)] <- 3
  out[out == 0] <- NA
  return(out)
}

# set.seed(1)
# 
# testdata <- runif(100)*100
# 
# testdata[1] <- NA
# 
# classify_SOC(testdata)
# 
# plot(testdata, classify_SOC(testdata))
# cbind(testdata, classify_SOC(testdata))

# END