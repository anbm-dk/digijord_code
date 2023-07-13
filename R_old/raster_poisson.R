library(terra)

# Create a SpatRaster from scratch
x <- rast(nrows = 108, ncols = 21, xmin = 0, xmax = 10)

# Create a SpatRaster from a file
f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)

fun1 <- function(x) {
  out <- x * 0 + rpois(1, 1)
  return(out)
}

r_poisson <- app(r, fun = fun1, cores = 2)

plot(r_poisson)
