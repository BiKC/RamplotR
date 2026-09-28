# Pure functions used by Shiny and the scientific regression tests.
ram_reference_group <- function(resn, next_resn = rep(NA_character_, length(resn)),
                                bonded_to_next = rep(FALSE, length(resn))) {
  stopifnot(length(resn) == length(next_resn),
            length(resn) == length(bonded_to_next))
  group <- rep("General", length(resn))
  group[!is.na(resn) & resn == "PRO"] <- "PRO"
  group[!is.na(resn) & resn == "GLY"] <- "GLY"
  prepro <- !is.na(resn) & !resn %in% c("PRO", "GLY") &
    !is.na(next_resn) & next_resn == "PRO" &
    !is.na(bonded_to_next) & bonded_to_next
  group[prepro] <- "preProline"
  group
}

# Classification is independent of the selected display background. The
# explicitly labelled legacy mode keeps old figures reproducible.
ram_classify_torsions <- function(torsions, reference_dir, selected_reference,
                                  mode = "residue", threshold_fn) {
  stopifnot(mode %in% c("residue", "legacy"))
  result <- torsions
  n <- nrow(result)
  result$reference_group <- ram_reference_group(
    result$resn,
    if ("next_resn" %in% names(result)) result$next_resn else rep(NA_character_, n),
    if ("bonded_to_next" %in% names(result)) result$bonded_to_next else rep(FALSE, n)
  )
  result$reference_used <- if (mode == "residue") result$reference_group else
    rep("Selected plotting background (legacy)", n)
  result$region <- rep(NA_character_, n)
  result$density <- rep(NA_real_, n)
  if (n == 0L) return(result)
  for (group in unique(result$reference_used)) {
    rows <- which(result$reference_used == group)
    ref <- if (mode == "legacy") selected_reference else
      ram_read_reference(file.path(reference_dir, group))
    limits <- ram_density_thresholds(ref)
    px <- match(round(result$phi[rows]), ref$x)
    py <- match(round(result$psi[rows]), ref$y)
    valid <- which(!is.na(px) & !is.na(py))
    if (!length(valid)) next
    positions <- rows[valid]
    z <- ref$z[cbind(py[valid], px[valid])]
    region <- rep("Not allowed", length(z))
    region[z > limits[3]] <- "Generously allowed"
    region[z > limits[2]] <- "Allowed"
    region[z > limits[1]] <- "Favoured"
    result$region[positions] <- region
    result$density[positions] <- ram_density_ranks(ref, z)
  }
  result
}

# Find contour cutoffs by cumulative probability mass, without recursive search.
# Ties are kept together; the result is the closest attainable coverage.
ram_density_thresholds <- function(reference, percentages = c(85, 98, 99.95)) {
  stopifnot(all(is.finite(percentages) & percentages >= 0 & percentages <= 100))
  z <- as.numeric(reference$z)
  if (!length(z) || any(!is.finite(z)) || any(z < 0) || sum(z) <= 0) {
    stop("The reference density grid must contain finite, non-negative values with positive total mass.")
  }
  levels <- sort(unique(z), decreasing = TRUE)
  counts <- tabulate(match(z, levels), nbins = length(levels))
  mass <- cumsum(levels * counts) / sum(z) * 100
  # Candidate cutoffs lie exactly on density levels. Because classification
  # uses z > cutoff, coverage at level k equals mass from higher levels.
  attainable <- c(0, head(mass, -1L))
  vapply(percentages, function(target) {
    index <- which.min(abs(attainable - target))
    levels[index]
  }, numeric(1))
}

# Reference grids are immutable while the application is running.
ram_read_reference <- local({
  cache <- new.env(parent = emptyenv())
  function(path) {
    key <- normalizePath(path, mustWork = TRUE)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, readRDS(key), envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
})

ram_density_ranks <- function(reference, values) {
  z <- as.numeric(reference$z)
  if (!length(z) || any(!is.finite(z)) || any(z < 0) || sum(z) <= 0) {
    stop("Invalid reference density grid")
  }
  levels <- sort(unique(z))
  weights <- levels * tabulate(match(z, levels), nbins = length(levels))
  larger_mass <- rev(cumsum(rev(weights))) - weights
  # findInterval maps arbitrary observed density to the greatest lower grid
  # level. If no grid density is lower, every positive grid value is higher.
  idx <- findInterval(values, levels)
  result <- numeric(length(values))
  result[idx == 0] <- 100
  has_level <- idx > 0
  if (any(has_level)) {
    result[has_level] <- larger_mass[idx[has_level]] / sum(z) * 100
    strictly_between <- has_level & values != levels[pmax(1L, idx)]
    if (any(strictly_between)) {
      result[strictly_between] <- (larger_mass[idx[strictly_between]] +
        weights[idx[strictly_between]]) / sum(z) * 100
    }
  }
  result
}
