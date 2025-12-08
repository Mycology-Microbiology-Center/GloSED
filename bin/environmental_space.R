
## GloSED
## Script to assess environmental coverage of soil sampling locations

## Workflow description:
# - The script assesses how well the sampled environmental space covers global environmental conditions
#   by creating binned 2D visualizations of environmental variable pairs
# - For each pair of environmental variables (e.g., climate - Bio01/Bio12, and soil - pH/N),
#   the function bins the global environmental space into a regular grid and determines
#   which bins contain at least one sampled location
# - In the output, the script produces heatmaps showing the global area share of each environmental bin
#   overlaid with sampled points to visualize coverage gaps

## Input data:
# - Environmental raster files:
#   - `Bio01.tif` (Annual Mean Temperature)
#   - `Bio12.tif` (Annual Precipitation)
#   - `pH.tif`    (soil pH)
#   - `N.tif`     (soil nitrogen content)
# - Sample metadata: `GloSED__SampleMetadata.xlsx` with geographic coordinates (`Latitude`, `Longitude`)

## Output data:
# - Plots in PDF format: `EnvCoverage__Bio01-vs-Bio12.pdf` and `EnvCoverage__pH-vs-N.pdf`
#   showing binned environmental coverage with sample overlay
# - Coverage statistics: percentage of global environmental bins represented by samples

library(terra)
library(data.table)
library(ggplot2)

theme_set(theme_classic(base_size = 12))

## Load rasters
rstpath <- "/path/to/rasters/"
Bio01 <- rast( file.path(rstpath, "Bio01.tif") )   # Annual Mean Temperature
Bio12 <- rast( file.path(rstpath, "Bio12.tif") )   # Annual Precipitation
N     <- rast( file.path(rstpath, "N.tif") )       # Soil nitrogen content
pH    <- rast( file.path(rstpath, "pH.tif") )      # Soil pH

## Create a raster stack
RST <- c(Bio01, Bio12, pH, N)

## Load sample metadata
meta <- readxl::read_excel(
  path = "GloSED__SampleMetadata.xlsx",
  na = c("", " ", "NA", "#N/A"),
  guess_max = 6000)

setDT(meta)
metad <- meta[ , .(SampleID, Latitude, Longitude)]
metad <- unique(metad, by = c("Longitude", "Latitude"))

## Extract predictors
xy <- vect(metad, geom=c("Longitude", "Latitude"), crs="+proj=longlat +datum=WGS84")
PREDICTORS <- extract(RST, xy)
npreds <- ncol(PREDICTORS) - 1

## Remove rows with all NAs
PREDICTORS <- cbind(metad, PREDICTORS)
PREDICTORS <- PREDICTORS[ rowSums(is.na(PREDICTORS)) < npreds ]

## Function to compute coverage in binned 2D environmental space
coverage_binned_2d <- function(
    r1, r2, obs_dt,
    v1_col  = names(r1),
    v2_col  = names(r2),
    nbins   = 50,
    binning = c("equalwidth", "quantile"),
    trim_q  = c(0.01, 0.99),
    overlay_points = TRUE,
    label_mode = c("midpoint", "range"),
    digits = 2,
    tick_step = NULL  # e.g., 2 to label every 2nd bin
    ) {
  
  stopifnot(inherits(r1, "SpatRaster"), inherits(r2, "SpatRaster"))
  binning <- match.arg(binning)

  ## Align & mask rasters
  r <- c(r1, r2)
  names(r) <- c("v1", "v2")
  r <- mask(r, !is.na(r[[1]]) & !is.na(r[[2]]))

  ## Extract values + cell areas
  vals <- as.data.frame(r, xy = FALSE, na.rm = TRUE)
  setDT(vals)

  ## Trim extremes by weighted quantiles for stability
  q1 <- quantile(vals$v1, probs = c(trim_q[1], trim_q[2]))
  q2 <- quantile(vals$v2, probs = c(trim_q[1], trim_q[2]))
  vals <- vals[v1 >= q1[1] & v1 <= q1[2] & v2 >= q2[1] & v2 <= q2[2]]

  ## Build bin breaks
  if(binning == "equalwidth") {
    b1 <- seq(min(vals$v1), max(vals$v1), length.out = nbins + 1)
    b2 <- seq(min(vals$v2), max(vals$v2), length.out = nbins + 1)
  } else {
    b1 <- quantile(vals$v1, probs = seq(0, 1, length.out = nbins + 1))
    b2 <- quantile(vals$v2, probs = seq(0, 1, length.out = nbins + 1))
    ## ensure strictly increasing (can happen with ties)
    b1 <- unique(b1); b2 <- unique(b2)
    if (length(b1) < 3L || length(b2) < 3L)
      stop("Not enough unique breaks after quantile binning. Try equalwidth or fewer bins.")
  }

  ## Bin global cells
  vals[, b1i := cut(v1, breaks = b1, include.lowest = TRUE, right = TRUE)]
  vals[, b2i := cut(v2, breaks = b2, include.lowest = TRUE, right = TRUE)]
  vals <- vals[!is.na(b1i) & !is.na(b2i)]
 
  ## Global availability (without area weighting)
  g_bins <- vals[, .(Availability = .N), by = .(b1i, b2i) ]
  g_bins[, area_share := Availability / sum(Availability) ]

  ## Prepare observed table (expects values already present)
  stopifnot(all(c(v1_col, v2_col) %in% names(obs_dt)))
  obs <- as.data.table(obs_dt)
  obs <- obs[!is.na(get(v1_col)) & !is.na(get(v2_col))]
  ## Keep only obs inside the trimmed ranges
  obs <- obs[get(v1_col) >= q1[1] & get(v1_col) <= q1[2] &
             get(v2_col) >= q2[1] & get(v2_col) <= q2[2]]

  ## Map observations to same bins
  obs[, b1i := cut(get(v1_col), breaks = b1, include.lowest = TRUE, right = TRUE)]
  obs[, b2i := cut(get(v2_col), breaks = b2, include.lowest = TRUE, right = TRUE)]
  obs <- obs[!is.na(b1i) & !is.na(b2i)]

  ## Determine coverage: bins with at least one observation
  covered <- unique(obs[, .(b1i, b2i)])
  g_bins[, covered := as.integer(paste(b1i, b2i) %in% paste(covered$b1i, covered$b2i))]

  ## Coverage stats (bin-level and area-weighted)
  n_bins_total   <- nrow(g_bins)
  n_bins_covered <- g_bins[, sum(covered)]
  bin_coverage   <- 100 * n_bins_covered / n_bins_total

  ## Bin indices
  g_bins[, `:=`(
    x = as.integer(factor(b1i, levels = levels(vals$b1i))),
    y = as.integer(factor(b2i, levels = levels(vals$b2i)))
  )]

  ## Build axis labels from breaks
  b1_mid <- head(b1, -1) + diff(b1)/2
  b2_mid <- head(b2, -1) + diff(b2)/2

  if(label_mode == "midpoint") {
    x_labels <- format(round(b1_mid, digits), trim = TRUE)
    y_labels <- format(round(b2_mid, digits), trim = TRUE)
  } else { # "range"
    fmt <- function(lo, hi) paste0("[", format(round(lo, digits), trim = TRUE),
                                   "-", format(round(hi, digits), trim = TRUE), "]")
    x_labels <- fmt(head(b1, -1), tail(b1, -1))
    y_labels <- fmt(head(b2, -1), tail(b2, -1))
  }

  ## Thin ticks if required
  x_breaks <- seq_along(x_labels)
  y_breaks <- seq_along(y_labels)
  if(!is.null(tick_step) && tick_step > 1) {
    x_keep <- seq(1, length(x_breaks), by = tick_step)
    y_keep <- seq(1, length(y_breaks), by = tick_step)
    x_labels <- x_labels[x_keep]; x_breaks <- x_breaks[x_keep]
    y_labels <- y_labels[y_keep]; y_breaks <- y_breaks[y_keep]
  }

  ## Plot
  p <- ggplot(g_bins, aes(x, y)) +
    geom_tile(aes(fill = area_share)) +
    scale_fill_viridis_c(
        name = "Global area share", trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::label_percent(accuracy = 0.00001)) +
    scale_x_continuous(
      name = paste0(names(r1), " (binned)"),
      breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
    scale_y_continuous(
      name = paste0(names(r2), " (binned)"),
      breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
    coord_equal() +
    theme(axis.text.x=element_text(angle=90))

  ## Overlay observed points in the same binned coordinate system
  if(overlay_points && nrow(obs)) {
    obs_plot <- unique(
        obs[, 
          .(x = as.integer(factor(b1i, levels = levels(vals$b1i)))),
           by = .(b1i, b2i)])
    obs_plot[, y := as.integer(factor(b2i, levels = levels(vals$b2i))) ]
    p <- p + 
      geom_point(
        data = obs_plot, aes(x, y),
        inherit.aes = FALSE,
        size = 0.6, color = "black", alpha = 0.4)
  }

  list(
    plot = p,
    stats = list(
      n_bins_total     = n_bins_total,
      n_bins_covered   = n_bins_covered,
      bin_coverage_pct = bin_coverage),
    breaks = list(v1 = b1, v2 = b2)
  )
}

## Bio01 vs Bio12
C1 <- coverage_binned_2d(
    r1 = Bio01,
    r2 = Bio12,
    obs_dt = PREDICTORS,
    nbins = 50, binning = "equalwidth",
    label_mode = "range")

## pH vs N
C2 <- coverage_binned_2d(
    r1 = pH,
    r2 = N,
    obs_dt = PREDICTORS,
    nbins = 50, binning = "equalwidth",
    label_mode = "range")

## Export plots
dir.create("EnvCoverage", showWarnings = FALSE)

ggsave(
    filename = file.path("EnvCoverage", "EnvCoverage__Bio01-vs-Bio12.pdf"),
    plot = C1$plot,
    width = 15, height = 16)

ggsave(
    filename = file.path("EnvCoverage", "EnvCoverage__pH-vs-N.pdf"),
    plot = C2$plot,
    width = 15, height = 16)
