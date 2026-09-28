# Run from the repository root: Rscript tests/profile.R
source(file.path("shinyRam", "R", "ramachandran.R"))
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
dir <- tempfile("ramplotr-profile-")
dir.create(dir)
on.exit(unlink(dir, recursive = TRUE), add = TRUE)

# Check exact grid-value percentile lookups on sparse, tied, zero-valued and
# larger synthetic density grids against the uncached calculations.
set.seed(142L)
datasets <- list(
  sparse = c(0, 0, 1, 3),
  tied = c(1, 1, 2, 6),
  dense = sample(c(0, 1e-6, 0.2, 1, 5), 16384L, replace = TRUE)
)
for (name in names(datasets)) {
  z <- datasets[[name]]
  side <- as.integer(sqrt(length(z)))
  reference <- list(x = seq_len(side), y = seq_len(side),
                    z = matrix(z, nrow = side))
  path <- file.path(dir, name)
  saveRDS(reference, path)
  profile <- ram_reference_profile(path)
  assert(identical(profile, ram_reference_profile(path)),
         paste("Reference profile cache changed for", name))
  assert(isTRUE(all.equal(profile$thresholds,
                          ram_density_thresholds(reference), tolerance = 1e-12)),
         paste("Cached density cutoffs differ for", name))
  values <- unique(z)
  observed <- profile$percentiles[match(values, profile$levels)]
  expected <- ram_density_ranks(reference, values)
  assert(isTRUE(all.equal(as.numeric(observed), as.numeric(expected),
                          tolerance = 1e-10)),
         paste("Cached percentile lookups differ for", name))
}

# The new cached path must reproduce the uncached, per-group classifications.
for (name in c("General", "GLY", "PRO", "preProline")) {
  file.copy(file.path(dir, "dense"), file.path(dir, name))
}
torsions <- data.frame(
  resn = c("ALA", "GLY", "PRO", "SER", "ALA"),
  next_resn = c(NA, NA, NA, "PRO", NA),
  bonded_to_next = c(FALSE, FALSE, FALSE, TRUE, FALSE),
  phi = c(2, 3, 4, 5, NA_real_),
  psi = c(2, 3, 4, 5, NA_real_)
)
selected_reference <- readRDS(file.path(dir, "General"))
res <- ram_classify_torsions(torsions, dir, selected_reference, "residue",
                             threshold_fn = ram_density_thresholds)
for (i in which(is.finite(torsions$phi) & is.finite(torsions$psi))) {
  ref <- readRDS(file.path(dir, res$reference_group[[i]]))
  value <- ref$z[torsions$psi[[i]], torsions$phi[[i]]]
  limits <- ram_density_thresholds(ref)
  expected_region <- if (value > limits[[1L]]) "Favoured" else
    if (value > limits[[2L]]) "Allowed" else
      if (value > limits[[3L]]) "Generously allowed" else "Not allowed"
  assert(identical(res$region[[i]], expected_region),
         "Cached classification changed a residue's region")
  assert(isTRUE(all.equal(res$density[[i]], ram_density_ranks(ref, value))),
         "Cached classification changed a residue's density percentile")
}
assert(is.na(res$region[[5L]]) && is.na(res$density[[5L]]),
       "Undefined torsion angles must remain unclassified")
unlink(dir, recursive = TRUE)
message("Cached reference profile tests passed")
