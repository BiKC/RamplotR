# Run from repository root: Rscript tests/classification.R
source(file.path("shinyRam", "R", "ramachandran.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

dir <- tempfile("ramplotr-reference-")
dir.create(dir)
reference <- function(z) list(x = c(0, 1), y = c(0, 1),
                               z = matrix(z, nrow = 2))
general <- reference(c(1, 1, 1, 97))
glycine <- reference(c(97, 1, 1, 1))
saveRDS(general, file.path(dir, "General"))
saveRDS(glycine, file.path(dir, "GLY"))
saveRDS(general, file.path(dir, "PRO"))
saveRDS(general, file.path(dir, "preProline"))

torsions <- data.frame(
  resn = c("ALA", "GLY", "ALA", "PRO"),
  chain = c("A", "A", "A", "B"),
  resi = c(1L, 2L, 3L, 1L),
  phi = c(1, 1, 1, NA_real_),
  psi = c(1, 1, 1, NA_real_),
  next_resn = c(NA_character_, NA_character_, "PRO", NA_character_),
  bonded_to_next = c(FALSE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

res <- ram_classify_torsions(torsions, dir, general, mode = "residue",
                             threshold_fn = function(reference, pct) 0)
assert(identical(res$reference_group,
                 c("General", "GLY", "preProline", "PRO")),
       "Residue grouping failed")
assert(identical(res$region, c("Favoured", "Not allowed", "Favoured", NA_character_)),
       "Residue-specific reference selection failed")
assert(is.na(res$density[4L]), "Missing torsions must remain unclassified")
assert(identical(res$reference_used[2L], "GLY"),
       "Scientific reference must be retained with the results")
legacy <- ram_classify_torsions(torsions, dir, general, mode = "legacy",
                               threshold_fn = function(reference, pct) 0)
assert(identical(legacy$region[2L], "Favoured"),
       "Legacy mode should keep the selected display reference")
assert(isTRUE(all.equal(unname(res$density[1L]), 0)),
       "High-density residue percentile should be zero")
assert(isTRUE(all.equal(unname(res$density[2L]), 97)),
       "Low-density glycine percentile should match the reference mass")

unlink(dir, recursive = TRUE)
message("Residue-aware classification tests passed")
