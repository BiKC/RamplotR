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
      readRDS(file.path(reference_dir, group))
    limits <- vapply(c(85, 98, 99.95),
                     function(pct) threshold_fn(ref, pct), numeric(1))
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
    total <- sum(ref$z)
    if (is.finite(total) && total > 0) {
      result$density[positions] <- vapply(z, function(value) {
        100 * sum(ref$z[ref$z > value]) / total
      }, numeric(1))
    }
  }
  result
}
