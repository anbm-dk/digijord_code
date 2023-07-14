# Code for weighted cubist


# NB: Weighted cubist

l <- getModelInfo("cubist")

cubist_weighted <- l$cubist

cubist_weighted$label <- "cubist_weighted"

cubist_weighted$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  if (!is.null(wts)) {
    out <- Cubist::cubist(x,
      y,
      committees = param$committees,
      weights = wts,
      ...
    )
  } else {
    out <- Cubist::cubist(x, y, committees = param$committees, ...)
  }
  if (last) out$tuneValue$neighbors <- param$neighbors
  out
}

cubist_weighted$predict <- function(modelFit, newdata, submodels = NULL) {
  out <- predict(modelFit,
    as.data.frame(newdata),
    neighbors = modelFit$tuneValue$neighbors
  )
  if (!is.null(submodels)) {
    tmp <- vector(mode = "list", length = nrow(submodels) + 1)
    tmp[[1]] <- out

    for (j in seq(along = submodels$neighbors)) {
      tmp[[j + 1]] <- predict(modelFit,
        as.data.frame(newdata),
        neighbors = submodels$neighbors[j]
      )
    }

    out <- tmp
  }
  out
}

cubist_weighted$tags <- c(l$cubist$tags, "Accepts Case Weights")


# END
