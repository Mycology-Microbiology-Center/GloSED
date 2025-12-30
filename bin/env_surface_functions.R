
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


