# Load and stack all rasters in a directory

loadandstack <- function(dir = NULL) {
  rlist <- list.files(dir, pattern = "tif$", full.names = TRUE)
  for (r in rlist)
  {
    name <- unlist(strsplit(r, "[.]"))[length(unlist(strsplit(r, "[.]"))) - 1]
    assign(name, raster(r))
  }
  output <- stack(rlist)
  return(output)
}
