# MolProbity/Phenix-style Rama8000 Ramachandran evaluation.
#
# The vendored score grids under static/rama8000 originate from the current
# cctbx/Phenix implementation. See static/rama8000/SOURCE.md for provenance
# and licensing. This file intentionally mirrors cctbx residue-class priority,
# periodic interpolation and score thresholds.

ram_rama8000_classes <- c(
  "general", "glycine", "cis-proline", "trans-proline",
  "pre-proline", "isoleucine or valine"
)

ram_rama8000_file <- c(
  "general" = "general.csv",
  "glycine" = "glycine.csv",
  "cis-proline" = "cis-proline.csv",
  "trans-proline" = "trans-proline.csv",
  "pre-proline" = "pre-proline.csv",
  "isoleucine or valine" = "ile-val.csv"
)

ram_rama8000_table <- local({
  cache <- new.env(parent = emptyenv())
  function(reference_dir, residue_class) {
    if (length(residue_class) != 1L || is.na(residue_class) ||
        !residue_class %in% ram_rama8000_classes) {
      stop("Unknown Rama8000 residue class.")
    }
    path <- file.path(reference_dir, unname(ram_rama8000_file[[residue_class]]))
    key <- normalizePath(path, mustWork = TRUE)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      table <- as.matrix(utils::read.csv(
        key, header = FALSE, check.names = FALSE
      ))
      storage.mode(table) <- "double"
      if (!identical(dim(table), c(180L, 180L)) || any(!is.finite(table))) {
        stop("Invalid Rama8000 reference table: ", basename(path))
      }
      assign(key, table, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
})

ram_rama8000_group <- function(resn, next_resn = rep(NA_character_, length(resn)),
                               bonded_to_next = rep(FALSE, length(resn)),
                               omega_prev = rep(NA_real_, length(resn))) {
  stopifnot(length(resn) == length(next_resn),
            length(resn) == length(bonded_to_next),
            length(resn) == length(omega_prev))
  resn <- toupper(as.character(resn))
  next_resn <- toupper(as.character(next_resn))
  out <- rep("general", length(resn))

  pro <- !is.na(resn) & resn == "PRO"
  # cctbx is_cislike_peptide(): -90 < omega < 90 is assigned cis-Pro.
  # Missing omega falls through to trans-Pro, matching current ramalyze.
  cis <- pro & is.finite(omega_prev) & omega_prev > -90 & omega_prev < 90
  out[pro] <- "trans-proline"
  out[cis] <- "cis-proline"

  gly <- !pro & !is.na(resn) & resn == "GLY"
  out[gly] <- "glycine"

  # cctbx gives pre-Pro priority over the Ile/Val class.
  prepro <- !pro & !gly & !is.na(next_resn) & next_resn == "PRO" &
    !is.na(bonded_to_next) & bonded_to_next
  out[prepro] <- "pre-proline"

  ileval <- !pro & !gly & !prepro & !is.na(resn) & resn %in% c("ILE", "VAL")
  out[ileval] <- "isoleucine or valine"
  out
}

ram_rama8000_bins <- function(value) {
  value <- as.numeric(value)
  while (value > 180) value <- value - 360
  while (value < -180) value <- value + 360
  lower <- floor(value)
  if ((as.integer(lower) %% 2L) == 0L) lower <- lower - 1
  higher <- ceiling(value)
  if ((as.integer(higher) %% 2L) == 0L) higher <- higher + 1
  if (lower == higher) higher <- higher + 2
  bin <- function(x) {
    b <- as.integer((x + 179) / 2)
    if (b > 179L) b <- b - 180L
    if (b < 0L) b <- b + 180L
    b + 1L
  }
  list(lower_value = lower, higher_value = higher,
       lower_bin = bin(lower), higher_bin = bin(higher),
       value = value)
}

ram_rama8000_score_one <- function(table, phi, psi) {
  if (!is.finite(phi) || !is.finite(psi)) return(NA_real_)
  px <- ram_rama8000_bins(phi)
  py <- ram_rama8000_bins(psi)
  x1 <- px$lower_value; x2 <- px$higher_value
  y1 <- py$lower_value; y2 <- py$higher_value
  x <- px$value; y <- py$value
  q11 <- table[px$lower_bin, py$lower_bin]
  q22 <- table[px$higher_bin, py$higher_bin]
  q12 <- table[px$lower_bin, py$higher_bin]
  q21 <- table[px$higher_bin, py$lower_bin]
  # Same bilinear interpolation used by cctbx rama_eval.h.
  ((q11 * (x2 - x) * (y2 - y)) +
   (q21 * (x - x1) * (y2 - y)) +
   (q12 * (x2 - x) * (y - y1)) +
   (q22 * (x - x1) * (y - y1))) / ((x2 - x1) * (y2 - y1))
}

ram_rama8000_region <- function(score, residue_class) {
  out <- rep(NA_character_, length(score))
  finite <- is.finite(score)
  out[finite & score >= 0.02] <- "Favored"
  allowed_cutoff <- rep(0.001, length(score))
  allowed_cutoff[residue_class == "general"] <- 0.0005
  allowed_cutoff[residue_class == "cis-proline"] <- 0.002
  out[finite & score < 0.02 & score >= allowed_cutoff] <- "Allowed"
  out[finite & score < allowed_cutoff] <- "Outlier"
  out
}

ram_rama8000_classify <- function(torsions, reference_dir) {
  result <- torsions
  n <- nrow(result)
  if (!"omega_prev" %in% names(result)) result$omega_prev <- rep(NA_real_, n)
  result$rama8000_group <- ram_rama8000_group(
    result$resn,
    if ("next_resn" %in% names(result)) result$next_resn else rep(NA_character_, n),
    if ("bonded_to_next" %in% names(result)) result$bonded_to_next else rep(FALSE, n),
    result$omega_prev
  )
  result$rama8000_score <- rep(NA_real_, n)
  if (n) {
    for (group in unique(result$rama8000_group)) {
      rows <- which(result$rama8000_group == group &
                    is.finite(result$phi) & is.finite(result$psi))
      if (!length(rows)) next
      table <- ram_rama8000_table(reference_dir, group)
      result$rama8000_score[rows] <- vapply(rows, function(i) {
        ram_rama8000_score_one(table, result$phi[[i]], result$psi[[i]])
      }, numeric(1))
    }
  }
  result$rama8000_region <- ram_rama8000_region(
    result$rama8000_score, result$rama8000_group
  )
  result
}
