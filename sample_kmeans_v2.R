# Function to sample a raster layer using kmeans

sample_kmeans <- function(input          = NULL   # Raster object
                          , clusters     = NULL   # Number of clusters
                          , ncells       = NULL   # Number of cells to extract, otherwise uses all cells
                          , use_xy       = FALSE  # Add xy coordinates as variables for clustering
                          , only_xy      = FALSE  # Use only xy coordinates
                          , weights      = NULL   # Raster layer with weights between 0 and 1
                          , scale        = TRUE   # Center and scale variables
                          , pca          = FALSE  # Use principal component analysis on variables
                          , tol_pca      = 0      # Tolerance for pca (remove PCs below threshold)
                          , n_pcs        = NULL   # Maximum number of principal components
                          , num_init     = 1      # See KMeans_rcpp
                          , max_iters    = 10     # See KMeans_rcpp
                          , initializer  = NULL   # See KMeans_rcpp
                          , CENTROIDS    = NULL   # See KMeans_rcpp
                          , tol_kmeans   = 1e-04  # See KMeans_rcpp
                          , tol_opt      = 0.3    # See KMeans_rcpp
                          , seed         = NULL   # See KMeans_rcpp
                          , MiniBatch    = FALSE  # Use MiniBatchKmeans (fast, less accurate)
                          , batch_size   = 10     # See MiniBatchKmeans
                          , init_frac    = 1      # See MiniBatchKmeans
                          , early_stop   = 10     # See MiniBatchKmeans
                          , filename_cl  = NULL   # File names for rasters with clusters (1) and distances (2)
                          , args_cl      = NULL   # List with arguments for writing raster
                          , filename_d   = NULL
                          , args_d       = NULL
                          # , round_d      = NULL   # Round distance raster to this number
                          , sp_pts       = FALSE  # Output locations as spatial points
                          , filename_pts = NULL   # Filename for output locations
                          , shp          = FALSE  # Write output locations as a shapefile
                          , args_pts     = NULL   # Arguments for writing output pointsx
                          , cores        = NULL   # Number of cpu cores to use
                          , m            = 2      # Number of blocks for clusterR
                          , verbose      = FALSE  # Print messages during processing
)
{
  backup_options <- options()
  options(error = traceback)  # Did this make it work? Find out how to reset options
  
  if (is.null(input) & is.null(weights)) { stop('No input data.') }
  
  require(magrittr)
  require(future.apply)
  require(fields)
  require(ClusterR)  # Not to be confused with the clusterR function
  
  if (verbose == TRUE) { message('Preparing input data.')  }
  
  if (is.null(input) & !is.null(weights))
  {
    message('No input variables. Using only coordinates.')
    only_xy <- TRUE
    input <- weights
  }
  
  # Extraction for rasters
  itisasraster <- is(input, "RasterLayer") |
    is(input, "SpatRaster") |
    is(input, "stars") |
    is(input, "RasterStack") |
    is(input, "RasterBrick")
  if (itisasraster)
  {
    if (!is.null(input) & !is.null(weights))
    {
      if (compareRaster(input, weights) == FALSE)
      {
        stop('Input and weights rasters do not match')
      }
      input <- stack(input, weights)
    }
    if (is.null(ncells)) {  # Use all cells if ncells is NULL
      df <- as.data.frame(input, xy = TRUE, na.rm = TRUE)
      ncells <- nrow(df)
    } else {
      if (ncells < clusters) {
        stop('ncells must be larger than the number of clusters')
      }
      df <- sampleRandom(input, ncells, xy = TRUE) %>% as.data.frame
      ncells <- nrow(df)
    }
    if (!is.null(weights))  # Weighted sampling
    {
      sampled <- sample(nrow(df)
                        , ncells
                        , prob = df[, ncol(df)]
                        , replace = TRUE)
      w <- df[sampled, ncol(df)]
      df <- df[sampled, -ncol(df)]
    }
  }
  
  # Extraction for spatial points
  if (class(input) == "SpatialPointsDataFrame")
  {
    if (!is.null(input) & !is.null(weights))
    {
      if (all.equal(input, weights) == FALSE)
      {
        stop('Input and weights coordinates do not match')
      }
      input@data <- cbind(input@data, weights@data)
    }
    df <- cbind(coordinates(input), input@data)
    if (is.null(ncells)) {  # Use all points if ncells is NULL
      ncells <- nrow(df)
    } else {
      if (ncells < clusters) {
        stop('ncells must be larger than the number of clusters')
      }
    }
    if (is.null(weights))
    {
      if (ncells > nrow(df))
      {
        message('ncells is larger than the number of input points. Using all input points instead.')
      } else {
        # standard sampling
        sampled <- sample(nrow(df)
                          , ncells
                          , replace = FALSE)
        df <- df[sampled, ]
      }
    } else {
      # weighted sampling
      sampled <- sample(nrow(df)
                        , ncells
                        , prob = df[, ncol(df)]
                        , replace = TRUE)
      w <- df[sampled, ncol(df)]
      df <- df[sampled, -ncol(df)]
    }
    
  }
  
  xy <- df[, 1:2]  # Extract coordinates
  
  if (only_xy == TRUE)
  {
    use_xy <- TRUE
    df <- xy
  }
  
  if (use_xy == FALSE) { df <- df[, -c(1:2)] }
  
  if (scale == TRUE)  # Scaling
  {
    if (verbose == TRUE) { message('Scaling input variables.') }
    if (!is.data.frame(df)) { df %<>% as.data.frame }
    means <- apply(df, 2, mean)
    sds <- apply(df, 2, sd)
    if(use_xy == TRUE)
    {
      sds[1:2] <- max(sds[1:2])
    }
    sds[sds == 0] <- 1
    if (verbose == TRUE)
    {
      scaling <- rbind(means, sds)
      rownames(scaling) <- c('Mean', 'SD')
      print(scaling)
    }
    df %<>%
      sweep(MARGIN = 2, STATS = means, check.margin = FALSE) %>%
      sweep(MARGIN = 2, FUN = "/", STATS = sds, check.margin = FALSE)
  }
  if (pca == TRUE)  # Principal components analysis
  {
    if (verbose == TRUE) { message('Conducting principal components analysis.') }
    
    if (is.null(n_pcs)) { n_pcs <- ncol(df) }
    pcs <- prcomp(df
                  , scale. = FALSE
                  , tol = tol_pca
                  , retx = TRUE
                  , rank. = n_pcs
    )
    df <- pcs$x
  }
  
  if (!is.data.frame(df)) { df %<>% as.data.frame }  # Make sure it's a dataframe
  
  out <- list()
  
  # Run kmeans
  if (verbose == TRUE) { message('Running k-means.') }
  
  if (is.null(seed)) { seed <- sample(10000, 1) }
  if (is.null(initializer)) { initializer  = "optimal_init" }
  
  if (MiniBatch == FALSE)
  {
    myclusters <- KMeans_rcpp(df
                              , clusters = clusters
                              , num_init = num_init
                              , max_iters = max_iters
                              , initializer = initializer
                              , CENTROIDS = CENTROIDS
                              , tol = tol_kmeans
                              , tol_optimal_init = tol_opt
                              , seed = seed
                              , verbose = verbose
    )
  } else {
    myclusters <- MiniBatchKmeans(df
                                  , clusters = clusters
                                  , num_init = num_init
                                  , max_iters = max_iters
                                  , initializer = initializer
                                  , CENTROIDS = CENTROIDS
                                  , tol = tol_kmeans
                                  , tol_optimal_init = tol_opt
                                  , seed = seed
                                  , batch_size = batch_size
                                  , init_fraction = init_frac
                                  , early_stop_iter = early_stop
                                  , verbose = verbose
    )
  }
  
  # Remove empty centroids
  if (myclusters$centroids %>% complete.cases %>% sum < clusters)
  {
    mycentroids <- myclusters$centroids %>%
      as.data.frame %>%
      drop_na
    missing <- clusters - nrow(mycentroids)
    
    message(paste0('A number of clusters ('
                   , missing
                   , ') contained no data. The function only produced '
                   , clusters - missing
                   , ' points.'))
  } else {
    mycentroids <- myclusters$centroids %>% as.data.frame
  }
  
  # Functions to map clusters
  if (pca == FALSE & scale == FALSE)
  {
    map_clusters_fun <- function(x)
    {
      if (x %>% sum %>% is.na) {
        out <- c(NA, NA)
      } else {
        dist <- x %>%
          matrix(1) %>%
          data.frame %>%
          rdist(mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
      
    }
  }
  if (pca == FALSE & scale == TRUE)
  {
    map_clusters_fun <- function(x)
    {
      if (x %>% sum %>% is.na) {
        out <- c(NA, NA)
      } else {
        dist <- x %>%
          '-'(means) %>%
          '/'(sds) %>%
          matrix(1) %>%
          data.frame %>%
          rdist(mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  if (pca == TRUE & scale == FALSE)
  {
    map_clusters_fun <- function(x)
    {
      if (x %>% sum %>% is.na) {
        out <- c(NA, NA)
      } else {
        x %<>%
          matrix(1) %>%
          data.frame
        colnames(x) <- pcs$rotation %>% rownames
        dist <- x %>%
          predict(pcs, .) %>%
          rdist(mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  if (pca == TRUE & scale == TRUE)
  {
    map_clusters_fun <- function(x)
    {
      if (x %>% sum %>% is.na) {
        out <- c(NA, NA)
      } else {
        x %<>%
          '-'(means) %>%
          '/'(sds) %>%
          matrix(1) %>%
          data.frame
        colnames(x) <- pcs$rotation %>% rownames
        dist <- x %>%
          predict(pcs, .) %>%
          rdist(mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  
  # Mapping procedure for rasters
  if (itisasraster)
  {
    # Create coordinate raster
    if (verbose == TRUE) { message('Producing coordinate raster.') }
    xy_r <- input %>%
      brick(values = FALSE, nl = 2) %>%
      setValues(values = coordinates(.))
    if (use_xy == TRUE) {
      if (only_xy == TRUE) {
        input <- xy_r
      } else {
        input <- raster::stack(xy_r, input)
      }
    }
    
    # Drop weights raster from the input
    if (!is.null(weights) & only_xy == FALSE)
    {
      input <- dropLayer(input, nlayers(input))
    }
    
    out$points <- NA
    
    if (verbose == TRUE) { message('Mapping clusters.') }
    if (is.null(cores))
    {
      out$clusters <- calc(input
                           , fun = map_clusters_fun
                           , forceapply = TRUE
      )
    } else {
      
      require(snow)
      
      closeAllConnections()  # Remove any previous connections in case they are still open.
      
      print(input) # debug
      print(map_clusters_fun) # debug
      
      beginCluster(cores)
      
      cl <- getCluster()
      
      export_this <- c('mycentroids')
      
      print(mycentroids) # debug
      
      if (exists('pcs')) { export_this <- c(export_this, 'pcs') 
      print(pcs) # debug
      }
      if (exists('means')) { export_this <- c(export_this, 'means', 'sds')
      print(means) # debug
      print(sds) # debug
      }
      
      print(export_this) # debug
      
      snow::clusterExport(cl
                          , export_this
                          , envir = environment())
      
      print(cl)
      
      out$clusters <- clusterR(input
                               , calc
                               , args = list(fun = map_clusters_fun
                                             , forceapply = TRUE)
                               , m = m
                               , verbose = TRUE
      )
      
      endCluster()
    }
    
    # Calculate weighted distances
    out$distances <- out$clusters[[2]]
    if (!is.null(weights))
    {
      if (verbose == TRUE) { message('Calculating weighted distances.') }
      s <- raster::stack(out$distances, weights)
      calc_wdist <- function(x)
      {
        if (x %>% sum %>% is.na)
        {
          out <- NA
        } else {
          if (x[2] == 0)
          {
            out <- NA
          } else {
            out <- x[1]/x[2]
          }
        }
        return(out)
      }
      
      if (is.null(cores))
      {
        wdist <- calc(s, fun = calc_wdist, forceapply = TRUE) 
      } else {
        closeAllConnections()
        
        beginCluster(cores)
        
        wdist <- clusterR(s
                          , overlay
                          , args = list(fun = calc_wdist
                                        , forceapply = TRUE))
        
        endCluster()
      }
      out$distances <- wdist
    }
    
    # Write clusters and distances to files if requested
    if (!is.null(filename_d))
    {
      if (verbose == TRUE) { message('Writing distance map to file.') }
      do.call(writeRaster
              , args = c(list(x = out$distances
                              , filename = filename_d)
                         , args_d)
      )
      out$distances <- raster(filename_d)
    }
    if (is.null(filename_cl))
    {
      out$clusters <- out$clusters[[1]]
    } else {
      if (verbose == TRUE) { message('Writing cluster map to file.') }
      do.call(writeRaster
              , args = c(list(x = out$clusters[[1]]
                              , filename = filename_cl)
                         , args_d)
      )
      out$clusters <- raster(filename_cl)
    }
    
    # Find the cluster centers
    if (verbose == TRUE) { message('Identifying cluster centers.') }
    zs <- zonal(out$distances, out$clusters, min)
    s <- stack(out$clusters, out$distances)
    
    findpoint <- function(x)
    {
      if (x %>% sum %>% is.na)
      {
        out <- NA
      } else {
        ismin <- zs[as.integer(x[1]), 2] == x[2]
        if (!ismin) { out <- NA } else { out <- as.integer(x[1]) }
      }
      return(out)
    }
    
    if (is.null(cores))
    {
      pts <- calc(s, fun = findpoint, forceapply = TRUE)
    } else {
      closeAllConnections()  # Remove any previous connections in case they are still open.
      
      beginCluster(cores)
      
      cl <- getCluster()
      
      snow::clusterExport(cl
                          , 'zs'
                          , envir = environment())
      
      pts <- clusterR(s
                      , calc
                      , args = list(fun = findpoint
                                    , forceapply = TRUE)
                      , m = m
                      , verbose = TRUE
      )
      
      endCluster()
    }
    
    out$points <- raster::rasterToPoints(pts, spatial = FALSE) %>%
      as.data.frame %>%
      dplyr::arrange(layer) %>%
      dplyr::rename(ID = layer)
    
    out$points <- out$points[!duplicated(out$points$ID), ]
  }
  
  # Mapping procedure for spatial points
  if (class(input) == "SpatialPointsDataFrame")
  {
    if (use_xy == TRUE) {
      if (only_xy == TRUE) {
        input@data <- coordinates(input)
      } else {
        input@data <- cbind(coordinates(input), input@data)
      }
    }
    
    # Drop weights from the input
    if (!is.null(weights) & only_xy == FALSE)
    {
      input@data <- input@data[, -ncol(input@data)]
    }
    
    out$points <- NA
    
    if (verbose == TRUE) { message('Mapping clusters.') }

    out$clusters <- apply(input@data
                          , 1
                          , FUN = map_clusters_fun
                          ) %>% t

    # Calculate weighted distances
    out$distances <- out$clusters[, 2] %>% unname
    out$clusters <- out$clusters[, 1] %>% unname
    if (!is.null(weights))
    {
      if (verbose == TRUE) { message('Calculating weighted distances.') }
      out$distances <- out$distances/weights@data[, 1]
    }
    
    # Find the cluster centers
    if (verbose == TRUE) { message('Identifying cluster centers.') }
    zs1 <- out$clusters %>%
      unique %>%
      sort
    
    zs2 <- sapply(zs1, function(x)
    {
      cl_min <- min(out$distances[out$clusters == x], na.rm = TRUE)
      return(cl_min)
    }
    )
    
    zs <- cbind(zs1, zs2)
    s <- cbind(out$clusters, out$distances)
    
    findpoint <- function(x)
    {
      if (x %>% sum %>% is.na)
      {
        out <- NA
      } else {
        ismin <- zs[as.integer(x[1]), 2] == x[2]
        if (!ismin) { out <- NA } else { out <- as.integer(x[1]) }
      }
      return(out)
    }
    
    pts <- apply(s, 1, FUN = findpoint)
    
    out$points <- cbind(coordinates(input), pts, c(1:nrow(input))) %>%
      as.data.frame %>%
      drop_na
    
    colnames(out$points) <- c('x', 'y', 'ID', 'Index')
    
    out$points <- out$points[!duplicated(out$points$ID), ] %>%
      arrange(ID)
  } 
    
  print(out$points)  # Debug
  
  # Write points to file if requested
  if (!is.null(filename_pts))
  {
    if (verbose == TRUE) { message('Writing points to file.') }
    require(tools)
    if (file_ext(filename_pts) == 'shp') { shp <- TRUE }
    if (shp == FALSE)
    {
      do.call(write.table
              , c(list(x = out$points
                       , file = filename_pts
              )
              , args_pts))
    } else {
      points_sp <- out$points
      coordinates(points_sp) <- ~ x + y
      proj4string(points_sp) <- crs(input)
      do.call(shapefile
              , c(list(x = points_sp
                       , filename = filename_pts
              )
              , args_pts))
    }
  }
  if (sp_pts == TRUE)
  {
    coordinates(out$points) <- ~ x + y
    proj4string(out$points) <- crs(input)
  }
  
  options(backup_options)
  
  return(out)
}

# END