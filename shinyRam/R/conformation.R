# Coarse, comparison-oriented backbone conformational states.
#
# These labels are deliberately separate from RamplotR density regions,
# Rama8000 validation and DSSP/secondary-structure assignment. They provide a
# compact way to ask whether independent structures/models occupy the same
# broad area of phi/psi space.
#
# Each finite phi/psi pair is assigned to the nearest canonical prototype on
# the periodic Ramachandran torus, provided it lies within max_distance
# degrees in wrapped Euclidean phi/psi distance. More remote positions are
# labelled "Other"; missing angles remain NA.
ram_backbone_basin_centers <- data.frame(
  basin = c("Alpha-R", "Beta", "PPII", "Alpha-L"),
  phi = c(-63, -135, -75, 60),
  psi = c(-43, 135, 145, 40),
  stringsAsFactors = FALSE
)

ram_backbone_basin <- function(phi, psi, max_distance = 70) {
  phi <- as.numeric(phi)
  psi <- as.numeric(psi)
  if (length(phi) != length(psi))
    stop("Phi and psi must have the same length.", call. = FALSE)
  max_distance <- suppressWarnings(as.numeric(max_distance))
  if (length(max_distance) != 1L || !is.finite(max_distance) ||
      max_distance <= 0 || max_distance > 180)
    stop("Backbone-basin max_distance must be between 0 and 180 degrees.",
         call. = FALSE)

  out <- rep(NA_character_, length(phi))
  valid <- which(is.finite(phi) & is.finite(psi))
  if (!length(valid)) return(out)

  wrapped <- function(value, center)
    ((value - center + 180) %% 360) - 180

  distances <- vapply(seq_len(nrow(ram_backbone_basin_centers)), function(i) {
    dphi <- wrapped(phi[valid], ram_backbone_basin_centers$phi[[i]])
    dpsi <- wrapped(psi[valid], ram_backbone_basin_centers$psi[[i]])
    sqrt(dphi^2 + dpsi^2)
  }, numeric(length(valid)))

  if (length(valid) == 1L)
    distances <- matrix(distances, nrow = 1L)

  nearest <- max.col(-distances, ties.method = "first")
  nearest_distance <- distances[cbind(seq_along(valid), nearest)]
  out[valid] <- ram_backbone_basin_centers$basin[nearest]
  out[valid[nearest_distance > max_distance]] <- "Other"
  out
}
