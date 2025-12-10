
## GloSED
## Script to compute environmental representativeness surface (dissimilarity index)

## Workflow:
# - Computes a dissimilarity index (DI) for each global grid cell relative to sampled training locations
# - Predictors are scaled (Z-score), weighted by variable importance, and Euclidean distance to the nearest
#   training point is normalised by the mean inter-sample distance
# - Low DI values indicate well-represented environments
#   high values indicate extrapolation beyond the sampled environmental space

## Inputs:
# - Training predictors and variable importance (`training_data.qs`, `variable_importance.qs`)
# - Global predictor tiles in Parquet format (`parquet/*.parquet`) with common extent, resolution and CRS
#   (e.g., Equal Earth), prepared beforehand by cropping, reprojecting and tiling rasters
# - Reference raster with target extent/CRS for output (e.g., `reprojected/Bio02.tif`)

## Outputs:
# - Per-tile DI parquet files written to `tmp_DI_tiles/`
# - Continuous DI raster `GloSED__DI.tif`
# - Binned novelty raster `GloSED__DI_binned.tif` (1 = inside training space, higher values = more novel)

## !NB. Approximate requirements for the run:
# - ~37 GB disk space (rasters, parquet tiles, intermediate files)
# - ~30–40 GB RAM (can be reduced via smaller tiles / fewer CPU threads)


## ----------------------------- Estimate dissimilarity index per tile, export to parquet

library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)
# library(plyr)
# library(FNN)
# library(qs)


## `prepare_train_data` and `compute_DI` functions
source("bin/env_surface_functions.R")

## Load training data and variable importance
EXTR <- qs::qread(file = "training_data.qs")
VI   <- qs::qread(file = "variable_importance.qs")

## Normalize importance weights
VI[ , Importance := XGImportance / max(XGImportance) ]

ggplot(data = VI, aes(y = Predictor, x = Importance)) +
  geom_point(size = 2.5) +
  geom_segment(aes(yend = Predictor, xend = 0))

## Convert importances to a named vector
var_importance <- VI$Importance
names(var_importance) <- VI$Predictor

## Prepare training data
train_prep <- prepare_train_data(
  train_data = EXTR,
  var_importance = var_importance)

rm(EXTR); gc()

qs::qsave(train_prep, file = "train_prep.qs",
  preset = "custom", algorithm = "zstd", compress_level = 16L, nthreads = 1L)


## Predictor tiles
prq <- list.files(
  path = "parquet/",
  pattern = "*.parquet", full.names = TRUE)


library(doFuture)
library(progressr)
handlers(global = TRUE)
registerDoFuture()
# plan(multisession, workers = 5)      # for RStudio
plan(multicore, workers = 5)           # will crash RStudio !!
options(future.globals.maxSize = 5e9)
options(doFuture.rng.onMisuse = "ignore")

## Process chunks (save results to files, as too large to fit in mem)
## By default, tiles are processed in chunks of 50,000 rows
plyr::a_ply(
  .data      = prq,
  .margins   = 1,
  .fun       = function(x){
    compute_DI(
      parquet_file = x,
      train_prep   = train_prep,
      just_DI = TRUE,
      to_file = TRUE)
  },
  .parallel = TRUE, .progress = "progressr")


## ----------------------------- Create DI raster


library(data.table)
library(arrow)
library(dplyr)
library(terra)

## Load results
RES <- arrow::open_dataset("tmp_DI_tiles/") %>%
  filter(! is.na(DI) ) %>% 
  collect() %>% 
  setDT()

## Reference raster
REF <- rast("reprojected/Bio02.tif")

## Re-create raster (this step can be RAM-intensive)
RST <- rast(x = RES, type = "xyz", crs = crs(REF), extent = ext(REF))
RST

writeRaster(
  x = RST,
  filename = "GloSED__DI.tif",
  overwrite = TRUE,
  gdal = c("COMPRESS=ZSTD", "PREDICTOR=2", "ZSTD_LEVEL=8"))



## ----------------------------- Bin DI index

library(data.table)
library(terra)

## Load the "training data"
train_prep <- qs::qread(file = "train_prep.qs")

## Load DI raster
RST <- rast("GloSED__DI.tif")

## Observed range of DI values
minmax(RST, compute = TRUE)

terra::global(x = RST, fun = quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

## AOA-like threshold
## For each neighboar of a training point, estimate DI (distance / divided by trainDist_avrgmean),
## thresold = Q3 + 1.5*IQR
DI_threshold <- function(train_prep) {
  stopifnot(inherits(train_prep, "train_env"))

  W <- as.matrix(train_prep$train_weighted)
  # k=2: neighbour 1 = self (0 dist), neighbour 2 = nearest other point
  nn <- FNN::get.knnx(data = W, query = W, k = 2)
  DI_train <- nn$nn.dist[, 2] / train_prep$trainDist_avrgmean
  summary(DI_train)
 
  threshold <- boxplot.stats(DI_train)$stats[5]  # the extreme of the upper whisker

  return(threshold)
}

THR <- DI_threshold(train_prep)


## Multi-band "novelty"
## 1 = inside, 2...5 = increasingly "novel"
breaks <- c(-Inf, THR, 2*THR, 4*THR, 6*THR, Inf)

RR <- classify(
  x = RST,
  rcl = cbind(breaks[-length(breaks)], breaks[-1], 1:(length(breaks)-1) ))

## Convert to factor
RR <- terra::as.factor(RR)

plot(RR)

writeRaster(
  x = RR,
  filename = "GloSED__DI_binned.tif",
  overwrite = TRUE,
  gdal = c("COMPRESS=ZSTD", "PREDICTOR=2", "ZSTD_LEVEL=8"))

