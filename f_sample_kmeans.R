# Function to sample a raster layer using kmeans

sample_kmeans <- function(
    input = NULL # Raster object
    , clusters = 3 # Number of clusters
    , ncells = NULL # Number of cells to extract, otherwise uses all cells
    , use_xy = FALSE # Add xy coordinates as variables for clustering
    , only_xy = FALSE # Use only xy coordinates
    , weights = NULL # Raster layer or numeric vector with weights between 0 and 1
    , scale = TRUE # Center and scale variables
    , pca = FALSE # Use principal component analysis on variables
    , tol_pca = 0 # Tolerance for pca (remove PCs below threshold)
    , n_pcs = NULL # Maximum number of principal components
    , num_init = 1 # See KMeans_rcpp
    , max_iters = 10 # See KMeans_rcpp
    , initializer = NULL # See KMeans_rcpp
    , CENTROIDS = NULL # See KMeans_rcpp
    , tol_kmeans = 1e-04 # See KMeans_rcpp
    , tol_opt = 0.3 # See KMeans_rcpp
    , seed = NULL # See KMeans_rcpp
    , MiniBatch = FALSE # Use MiniBatchKmeans (fast, less accurate)
    , batch_size = 10 # See MiniBatchKmeans
    , init_frac = 1 # See MiniBatchKmeans
    , early_stop = 10 # See MiniBatchKmeans
    , filename_cl = NULL # File names for rasters with clusters (1) and distances (2)
    , args_cl = NULL # List with arguments for writing raster
    , filename_d = NULL,
    args_d = NULL
    # , round_d      = NULL   # Round distance raster to this number
    , sp_pts = FALSE # Output locations as spatial points
    , filename_pts = NULL # Filename for output locations
    , shp = FALSE # Write output locations as a shapefile
    , args_pts = NULL # Arguments for writing output pointsx
    , cores = NULL # Number of cpu cores to use
    # , m = 2 # Number of blocks for clusterR
    , verbose = FALSE # Print messages during processing
    ) {
  backup_options <- options()
  options(error = traceback) # Did this make it work? Find out how to reset options

  if (is.null(input) & is.null(weights)) {
    stop("No input data.")
  }

  require(magrittr)
  require(future.apply)
  require(fields)
  require(ClusterR) # Not to be confused with the clusterR function

  if (verbose == TRUE) {
    message("Preparing input data.")
  }

  if (is.null(input) & !is.null(weights)) {
    message("No input variables. Using only coordinates.")
    only_xy <- TRUE
    input <- weights
  }

  # Extraction for rasters
  itisasraster <- is(input, "SpatRaster")
  
  itispoints <- FALSE
  if (is(input, "SpatVector")) {
    if (terra::geomtype(input) == "points") {
      itispoints <- TRUE
    }
  }
  
  if (only_xy) {
    n_vars <- 2
  } else {
    if (itisasraster) {
      n_vars <- terra::nlyr(input)
    }
    if (itispoints) {
      n_vars <- ncol(input)
    }
    if (use_xy) {
      n_vars <- n_vars + 2
    }
  }
  
  if (itisasraster) {
    if (!is.null(input) & !is.null(weights)) {
      if (terra::compareGeom(input, weights) == FALSE) {
        stop("Input and weights rasters do not match")
      }
      input <- c(input, weights)
    }
    if (is.null(ncells)) { # Use all cells if ncells is NULL
      df <- as.data.frame(input, xy = TRUE, na.rm = TRUE)
      ncells <- nrow(df)
    } else {
      if (ncells < clusters) {
        stop("ncells must be larger than the number of clusters")
      }
      df <- terra::spatSample(input, ncells, xy = TRUE, as.df = TRUE)
      ncells <- nrow(df)
    }
    if (!is.null(weights)) # Weighted sampling
      {
        sampled <- sample(nrow(df),
          ncells,
          prob = df[, ncol(df)],
          replace = TRUE
        )
        w <- df[sampled, ncol(df)]
        df <- df[sampled, -ncol(df)]
      }
  }

  # Extraction for spatial points
  if (itispoints) {
    if (!is.null(input) & !is.null(weights)) {
      if (length(input) != length(weights)) {
        stop("Input and weights coordinates do not match")
      }
      input$weights <- weights
    }
    df <- cbind(terra::crds(input), terra::values(input))
    if (is.null(ncells)) { # Use all points if ncells is NULL
      ncells <- nrow(df)
    } else {
      if (ncells < clusters) {
        stop("ncells must be larger than the number of clusters")
      }
    }
    if (is.null(weights)) {
      if (ncells > nrow(df)) {
        message("ncells is larger than the number of input points. Using all input points instead.")
      } else {
        # standard sampling
        sampled <- sample(nrow(df),
          ncells,
          replace = FALSE
        )
        df <- df[sampled, ]
      }
    } else {
      # weighted sampling
      sampled <- sample(nrow(df),
        ncells,
        prob = df[, ncol(df)],
        replace = TRUE
      )
      w <- df[sampled, ncol(df)]
      df <- df[sampled, -ncol(df)]
    }
  }

  xy <- df[, 1:2] # Extract coordinates

  if (only_xy == TRUE) {
    use_xy <- TRUE
    df <- xy
  }

  if (use_xy == FALSE) {
    df <- df[, -c(1:2)]
  }

  if (scale == TRUE) # Scaling
    {
      if (verbose == TRUE) {
        message("Scaling input variables.")
      }
      if (!is.data.frame(df)) {
        df %<>% as.data.frame
      }
      means <- apply(df, 2, mean)
      sds <- apply(df, 2, sd)
      if (use_xy == TRUE) {
        sds[1:2] <- max(sds[1:2])
      }
      sds[sds == 0] <- 1
      if (verbose == TRUE) {
        scaling <- rbind(means, sds)
        rownames(scaling) <- c("Mean", "SD")
        print(scaling)
      }
      df %<>%
        sweep(MARGIN = 2, STATS = means, check.margin = FALSE) %>%
        sweep(MARGIN = 2, FUN = "/", STATS = sds, check.margin = FALSE)
    }
  if (pca == TRUE) # Principal components analysis
    {
      if (verbose == TRUE) {
        message("Conducting principal components analysis.")
      }

      if (is.null(n_pcs)) {
        n_pcs <- ncol(df)
      }
      pcs <- prcomp(df,
        scale. = FALSE,
        tol = tol_pca,
        retx = TRUE,
        rank. = n_pcs
      )
      df <- pcs$x
    }

  if (!is.data.frame(df)) {
    df %<>% as.data.frame
  } # Make sure it's a dataframe

  out <- list()

  # Run kmeans
  if (verbose == TRUE) {
    message("Running k-means.")
  }

  if (is.null(seed)) {
    seed <- sample(10000, 1)
  }
  if (is.null(initializer)) {
    initializer <- "optimal_init"
  }

  if (MiniBatch == FALSE) {
    myclusters <- ClusterR::KMeans_rcpp(
      df,
      clusters = clusters,
      num_init = num_init,
      max_iters = max_iters,
      initializer = initializer,
      CENTROIDS = CENTROIDS,
      tol = tol_kmeans,
      tol_optimal_init = tol_opt,
      seed = seed,
      verbose = verbose
    )
  } else {
    myclusters <- ClusterR::MiniBatchKmeans(
      df,
      clusters = clusters,
      num_init = num_init,
      max_iters = max_iters,
      initializer = initializer,
      CENTROIDS = CENTROIDS,
      tol = tol_kmeans,
      tol_optimal_init = tol_opt,
      seed = seed,
      batch_size = batch_size,
      init_fraction = init_frac,
      early_stop_iter = early_stop,
      verbose = verbose
    )
  }

  # Remove empty centroids
  if (myclusters$centroids %>% complete.cases() %>% sum() < clusters) {
    mycentroids <- myclusters$centroids %>%
      as.data.frame() %>%
      drop_na()
    missing <- clusters - nrow(mycentroids)

    message(paste0(
      "A number of clusters (",
      missing,
      ") contained no data. The function only produced ",
      clusters - missing,
      " points."
    ))
  } else {
    mycentroids <- myclusters$centroids %>% as.data.frame()
  }

  # Functions to map clusters
  if (pca == FALSE & scale == FALSE) {
    map_clusters_fun <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- c(NA, NA)
      } else {
        dist <- x %>%
          matrix(1) %>%
          data.frame() %>%
          fields::rdist(., mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  if (pca == FALSE & scale == TRUE) {
    map_clusters_fun <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- c(NA, NA)
      } else {
        dist <- x %>%
          "-"(means) %>%
          "/"(sds) %>%
          matrix(1) %>%
          data.frame() %>%
          fields::rdist(., mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  if (pca == TRUE & scale == FALSE) {
    map_clusters_fun <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- c(NA, NA)
      } else {
        x %<>%
          matrix(1) %>%
          data.frame()
        colnames(x) <- pcs$rotation %>% rownames()
        dist <- x %>%
          predict(pcs, .) %>%
          fields::rdist(., mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }
  if (pca == TRUE & scale == TRUE) {
    map_clusters_fun <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- c(NA, NA)
      } else {
        x %<>%
          "-"(means) %>%
          "/"(sds) %>%
          matrix(1) %>%
          data.frame()
        colnames(x) <- pcs$rotation %>% rownames()
        dist <- x %>%
          predict(pcs, .) %>%
          fields::rdist(., mycentroids)
        out <- c(which.min(dist), min(dist, na.rm = TRUE))
      }
      return(out)
    }
  }

  # Mapping procedure for rasters
  if (itisasraster) {
    # Create coordinate raster
    if (verbose == TRUE) {
      message("Producing coordinate raster.")
    }
    xy_r <- c(
      terra::init(input, "x"),
      terra::init(input, "y")
    )
    if (use_xy == TRUE) {
      if (only_xy == TRUE) {
        input <- xy_r
      } else {
        input <- c(xy_r, input)
      }
    }

    # Drop weights raster from the input
    if (terra::nlyr(input) > n_vars) {
      input <- terra::subset(input, c(1:n_vars))
    }

    out$points <- NA

    if (verbose == TRUE) {
      message("Mapping clusters.")
    }
    if (is.null(cores)) {
      out$clusters <- terra::app(
        input,
        fun = map_clusters_fun
      )
    } else {
      showConnections()
      
      cl <- parallel::makeCluster(cores)
      
      parallel::clusterEvalQ(
        cl,
        {
          library(fields)
          library(magrittr)
        }
      )
      
      export_this <- c("mycentroids")

      if (exists("pcs")) {
        export_this <- c(export_this, "pcs")
      }
      if (exists("means")) {
        export_this <- c(export_this, "means", "sds")
      }

      parallel::clusterExport(
        cl,
        c(
          export_this
        ),
        envir = environment()
      )
      
      out$clusters <- terra::app(
        input,
        fun = map_clusters_fun,
        cores = cl
      )
      
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      rm(cl)
    }
    
    # Calculate weighted distances
    out$distances <- out$clusters[[2]]
    out$clusters <- out$clusters[[1]]
    if (!is.null(weights)) {
      if (verbose == TRUE) {
        message("Calculating weighted distances.")
      }
      s <- c(out$distances, weights)
      calc_wdist <- function(x) {
        if (x %>% sum() %>% is.na()) {
          out <- NA
        } else {
          if (x[2] == 0) {
            out <- NA
          } else {
            out <- x[1] / x[2]
          }
        }
        return(out)
      }

      if (is.null(cores)) {
        wdist <- terra::app(s, fun = calc_wdist)
      } else {
        cl <- parallel::makeCluster(cores)
        
        parallel::clusterEvalQ(
          cl,
          {
            library(magrittr)
          }
        )
        
        wdist <- terra::app(s, fun = calc_wdist, cores = cl)

        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
        rm(cl)
      }
      out$distances <- wdist
    }
    
    names(out$distances) <- "distance"
    names(out$clusters) <- "cluster"

    # Write clusters and distances to files if requested
    if (!is.null(filename_d)) {
      if (verbose == TRUE) {
        message("Writing distance map to file.")
      }
      do.call(
        terra::writeRaster,
        args = c(
          list(
            x = out$distances,
            filename = filename_d
          ),
          args_d
        )
      )
      out$distances <- terra::rast(filename_d)
    }
    if (is.null(filename_cl)) {
      out$clusters <- out$clusters[[1]]
    } else {
      if (verbose == TRUE) {
        message("Writing cluster map to file.")
      }
      do.call(
        terra::writeRaster,
        args = c(
          list(
            x = out$clusters[[1]],
            filename = filename_cl
          ),
          args_d
        )
      )
      out$clusters <- terra::rast(filename_cl)
    }

    # Find the cluster centers
    if (verbose == TRUE) {
      message("Identifying cluster centers.")
    }
    zs <- terra::zonal(out$distances, out$clusters, "min", na.rm = TRUE)
    s <- c(out$clusters, out$distances)

    findpoint <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- NA
      } else {
        ismin <- zs[as.integer(x[1]), 2] == x[2]
        if (!ismin) {
          out <- NA
        } else {
          out <- as.integer(x[1])
        }
      }
      return(out)
    }

    if (is.null(cores)) {
      pts <- terra::app(s, fun = findpoint)
    } else {
      showConnections()
      
      cl <- parallel::makeCluster(cores)
      
      parallel::clusterEvalQ(
        cl,
        {
          library(magrittr)
        }
      )

      parallel::clusterExport(cl,
        "zs",
        envir = environment()
      )

      pts <- terra::app(s, fun = findpoint, cores = cl)

      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      rm(cl)
    }
    
    names(pts) <- "ID"

    out$points <- as.data.frame(pts, xy = TRUE, na.rm = TRUE) %>%
      dplyr::arrange(ID)

    out$points <- out$points[!duplicated(out$points$ID), ]
  }

  # Mapping procedure for spatial points
  if (itispoints) {
    if (use_xy == TRUE) {
      if (only_xy == TRUE) {
        terra::values(input) <- terra::crds(input)
      } else {
        terra::values(input) <- cbind(terra::crds(input), terra::values(input))
      }
    }

    # Drop weights from the input
    if (ncol(terra::values(input)) > n_vars) {
      terra::values(input) <- terra::values(input)[, c(1:n_vars)]
    }

    out$points <- NA

    if (verbose == TRUE) {
      message("Mapping clusters.")
    }

    out$clusters <- apply(
      terra::values(input),
      1,
      FUN = map_clusters_fun
    ) %>% t()

    # Calculate weighted distances
    out$distances <- out$clusters[, 2] %>% unname()
    out$clusters <- out$clusters[, 1] %>% unname()
    if (!is.null(weights)) {
      if (verbose == TRUE) {
        message("Calculating weighted distances.")
      }
      out$distances <- out$distances / weights
    }

    # Find the cluster centers
    if (verbose == TRUE) {
      message("Identifying cluster centers.")
    }
    zs1 <- out$clusters %>%
      unique() %>%
      sort()

    zs2 <- sapply(zs1, function(x) {
      cl_min <- min(out$distances[out$clusters == x], na.rm = TRUE)
      return(cl_min)
    })

    zs <- cbind(zs1, zs2)
    s <- cbind(out$clusters, out$distances)

    findpoint <- function(x) {
      if (x %>% sum() %>% is.na()) {
        out <- NA
      } else {
        ismin <- zs[as.integer(x[1]), 2] == x[2]
        if (!ismin) {
          out <- NA
        } else {
          out <- as.integer(x[1])
        }
      }
      return(out)
    }

    pts <- apply(s, 1, FUN = findpoint)

    out$points <- cbind(terra::crds(input), pts, c(1:length(pts))) %>%
      as.data.frame() %>%
      stats::complete.cases()

    colnames(out$points) <- c("x", "y", "ID", "Index")

    out$points <- out$points[!duplicated(out$points$ID), ] %>%
      arrange(ID)
  }

  # Write points to file if requested
  if (!is.null(filename_pts)) {
    require(tools)
    if (tools::file_ext(filename_pts) == "shp") {
      shp <- TRUE
    }
  }
  
  if (shp == TRUE | sp_pts == TRUE) {
    points_sp <- vect(out$points, geom = c("x", "y"), crs(input))
  }
  
  if (!is.null(filename_pts)) {
    if (verbose == TRUE) {
      message("Writing points to file.")
    }
    
    if (shp == FALSE) {
      do.call(
        utils::write.table,
        c(
          list(
            x = out$points,
            file = filename_pts
          ),
          args_pts
        )
      )
    } else {
      
      do.call(
       terra::writeVector,
        c(
          list(
            x = points_sp,
            filename = filename_pts
          ),
          args_pts
        )
      )
    }
  }
  if (sp_pts == TRUE) {
    out$points <- points_sp
  }

  options(backup_options)

  return(out)
}

# END
