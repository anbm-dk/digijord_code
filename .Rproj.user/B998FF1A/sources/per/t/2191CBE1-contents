# Function to produce predictions where some, but not all covariates are NA

predict_passna <- function(mod, dat, ...) {
  library(caret)
  library(Cubist)
  
  rfun2 <- function(mod2, dat2, ...) {
    notallnas <- rowSums(is.na(dat2)) < (ncol(dat2) - 2)  # NB: Edit this line
    out2 <- rep(NA, nrow(dat2))
    if (sum(notallnas) > 0) {
      out2[notallnas] <- predict(
        object = mod2,
        newdata = dat2[notallnas, ],
        na.action = na.pass,
        ...
      )
    }
    return(out2)
  }
  
  out <- rfun2(mod, dat, ...)
  return(out)
}

# END