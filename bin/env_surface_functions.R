
## GloSED auxillary script
## Source functions used in the `environmental_surface.R` script
##
## This implementation is based on the Dissimilarity Index (DI) approach from the CAST package (Meyer & Pebesma, 2021).
## This implementation is optimized for chunk-based processing of large raster datasets
## (it pre-computes training data scaling/weighting to avoid repeated calculations)


## Function to prepare training data for environmental representativeness computation
prepare_train_data <- function(train_data, var_importance) {
  # train_data = data.table containing predictor variables for the set of observed (training) samples
  # var_importance = named numeric vector of variable importance values
  #   (NB! importances are used as-is and are NOT normalised internally)

  setDT(train_data)
  train_dt <- copy(train_data)  # to avoid modifying the original data

  ## Check for sufficient training samples
  if (nrow(train_dt) < 2) {
    stop("At least two training points are required to compute inter-sample distances")
  }

  ## Verify that var_importance is a named numeric vector
  if (!is.numeric(var_importance) || is.null(names(var_importance))) {
    stop("var_importance must be a named numeric vector")
  }

  ## Keep only the variables present in both train data and the importance vector
  vars <- intersect(names(var_importance), names(train_dt))
  if (length(vars) == 0) {
    stop("No matching variables found between train_data and var_importance")
  }

  ## Reorder the importance vector to match the variables order
  weights <- var_importance[ vars ]

  ## Extract matrix of predictor values
  train_mat <- train_dt[, ..vars]

  ## Exclude rows with NAs
  cc <- complete.cases(train_mat)
  if(any(!cc)){
    cat("WARNING: NAs found in training data - these rows will be excluded\n")
    cat("Number of rows with NAs:", sum(!cc), "\n")
    train_mat <- train_mat[ cc ]
  }

  ## Standardise training data (centre and scale)
  train_mean <- colMeans(train_mat, na.rm = TRUE)
  train_sd   <- apply(train_mat, 2, stats::sd, na.rm = TRUE)

  ## Avoid division by zero if a predictor has zero variance
  train_sd[ train_sd == 0 ] <- 1

  ## Scale training data
  train_scaled <- sweep(train_mat, 2, train_mean, FUN = "-")
  train_scaled <- sweep(train_scaled, 2, train_sd,   FUN = "/")

  ## Multiply predictors by their variable importance weights
  train_weighted <- sweep(train_scaled, 2, weights, FUN = "*")

  ## Pre-compute squared norms for training data
  train_norm2 <- rowSums(train_weighted^2)

  ## Compute mean inter-sample distance for normalization
  train_dist_vec     <- dist(train_weighted)
  trainDist_avrgmean <- mean(train_dist_vec)

  ## Return pre-processed training data object
  result <- list(
    train_weighted     = train_weighted,     # matrix of scaled and weighted training data
    train_norm2        = train_norm2,        # pre-computed squared norms for each training sample
    train_mean         = train_mean,         # column means used for scaling
    train_sd           = train_sd,           # column standard deviations used for scaling
    weights            = weights,            # variable importance weights in correct order
    vars               = vars,               # variable names in order used
    trainDist_avrgmean = trainDist_avrgmean, # mean inter-sample distance for normalization
    n_train            = nrow(train_mat)     # number of training samples
  )

  class(result) <- c("train_env", "list")
  return(result)
}



## Compute dissimilarity index (DI) for new data
compute_DI <- function(parquet_file, train_prep, just_DI = FALSE, to_file = TRUE) {
  # parquet_file = Name of a parquet file with numeric predictor variables for the locations where the representativeness should be evaluated
  # train_prep = pre-processed training data object created by `prepare_train_data()`
  # just_DI = logical. If FALSE, return the input data with an added DI column
  #                    If TRUE, return only DI values with coordinates
  # to_file = logical. If TRUE, write the result to parquet file

  ## NB! Rows with missing predictor values will receive `NA` 

  ## Create subdir for file results
  if(to_file == TRUE) {
    OUTDIR <- "tmp_DI_tiles"
    dir.create(path = OUTDIR, showWarnings = FALSE)
  }

  ## Verify train_prep is the correct object type
  if (!inherits(train_prep, "train_env")) {
    stop("train_prep must be an object created by prepare_train_data()")
  }

  ## Read parquet file with the new data
  new_data <- arrow::read_parquet(parquet_file)
  setDT(new_data)

  if(nrow(new_data) == 0) {
    stop("No data found in the new data file\n")
  }

  ## Check that required variables are present in new_data
  missing_vars <- setdiff(train_prep$vars, names(new_data))
  if (length(missing_vars) > 0) {
    stop("Variables missing in new_data: ", paste(missing_vars, collapse = ", "))
  }

  ## Extract matrix of predictor values in the correct order
  clz <- train_prep$vars

  ## Wrapper function to process new data
  process_new_data <- function(nd, just_DI = just_DI) {

    new_mat <- nd[, ..clz]

    ## Identify rows with complete cases (no NAs) to avoid propagating NAs
    cc <- complete.cases(new_mat)
    if(any(!cc)){
        cat("WARNING: NAs found in new data - these rows will be excluded\n")
        cat("Number of rows with NAs:", sum(!cc), "\n")
        new_mat <- new_mat[ cc ]
    }

    n_ok <- sum(cc)
    
    if(n_ok < 1){
        cat("WARNING: No complete cases found in new data\n")
        
        dummy_result <- nd[1, ]
        dummy_result[, DI := NA]
        dummy_result <- dummy_result[ -1, ]  # to create empty data.table with all required columns

        if(just_DI == TRUE) {
        clz <- c("x", "y", "DI")
        dummy_result <- dummy_result[ , ..clz ]
        }

        return(dummy_result)
    }

    ### TO DO - optimize sweep operations (use pure data.table to avoid copying data)
    ## Scale new data using pre-computed training parameters
    new_scaled <- sweep(new_mat, 2, train_prep$train_mean, FUN = "-")
    new_scaled <- sweep(new_scaled, 2, train_prep$train_sd, FUN = "/")

    ## Apply variable importance weights
    new_scaled <- sweep(new_scaled, 2, train_prep$weights, FUN = "*")

    ## Use Fast Nearest Neighbor search algo
    ## to find the minimum distance to any training point for each new point
    min_dist <- FNN::knnx.dist(
        data = as.matrix(train_prep$train_weighted),  # reference (training) data
        query = as.matrix(new_scaled),                # query (new) data  
        k = 1,                                        # find only the closest point
        algorithm = "brute"                           # use brute force algorithm
    )[, 1]

    ## Calculate the dissimilarity index by normalising the minimum distances
    DI <- min_dist / train_prep$trainDist_avrgmean

    ## Add DI to the data (only to complete cases)
    nd[ cc, DI := DI]

    ## Just keep coorinates if all predictors are not required
    if(just_DI == TRUE) {
      clz <- c("x", "y", "DI")
      nd <- nd[ , ..clz ]
    }

    return(nd)
  }                # end of `process_new_data`

  ## Split large datasets into chunks
  MAXN <- 50000
  num_chunks <- ceiling(nrow(new_data) / MAXN)
  chnk <- metagMisc::chunk_table(x = new_data, n = num_chunks, to_list = TRUE)
  
  ## Process chunks without writing to file, merge chunks into a single data.table
  if(to_file == FALSE) {

    res <- plyr::llply(.data = chnk, .fun = process_new_data, just_DI = just_DI)
    res <- rbindlist(res)
    return(res)

  ## Process chunks and write each chunk to a separate parquet file
  } else {

    if(length(chnk) == 1) { names(chnk) <- "1" }

    plyr::a_ply(
      .data = names(chnk),
      .margins = 1,
      .fun = function(x){
        rez <- process_new_data(nd = chnk[[ x ]], just_DI = just_DI)
        
        fname <- sub(pattern = ".parquet", replacement = "", x = basename(parquet_file))
        fname <- paste0("DI__", fname, "__", x, ".parquet")
        
        arrow::write_parquet(x = rez,
          sink = file.path(OUTDIR, fname),
          compression = "zstd",
          compression_level = 6)

        rm(rez)
        return(NULL)
      }, .progress = "none")

    return(NULL)
  }

}

