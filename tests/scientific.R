# Run from the repository root: Rscript tests/scientific.R
source(file.path("shinyRam", "R", "ramachandran.R"))
source(file.path("shinyRam", "R", "conformation.R"))
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
groups <- ram_reference_group(
  c("ALA", "GLY", "PRO", "SER", "THR"),
  c("PRO", "PRO", "ALA", "PRO", "PRO"),
  c(TRUE, TRUE, TRUE, FALSE, TRUE)
)
assert(identical(groups, c("preProline", "GLY", "PRO", "General", "preProline")),
       "Wrong residue-specific grouping")
empty <- data.frame(resn = character(), phi = numeric(), psi = numeric())
result <- ram_classify_torsions(empty, ".", NULL, "residue",
                               function(matrix, pct) 0)
assert(nrow(result) == 0L, "Empty structures should not fail")
assert(all(c("region", "reference_group", "reference_used", "density") %in%
           names(result)), "Classification output schema changed")
message("Residue-grouping regression tests passed")

# Ties, invalid references, empty input and deterministic coverages.
ref <- list(z = matrix(c(1, 1, 2, 6), nrow = 2), x = c(-1, 0), y = c(-1, 0))
targets <- ram_density_thresholds(ref, c(0, 60, 80, 100))
assert(identical(unname(targets), c(6, 2, 1, 1)),
       "Cumulative contour thresholds must preserve tied density levels")
assert(identical(unname(ram_density_thresholds(ref, c(60, 80))),
                 unname(ram_density_thresholds(ref, c(60, 80)))),
       "Contour calculations must be deterministic")
ranks <- ram_density_ranks(ref, c(6, 2, 1))
assert(isTRUE(all.equal(unname(ranks), c(0, 60, 80))),
       "Density percentile lookup must match cumulative mass")
invalid <- try(ram_density_thresholds(list(z = matrix(0, 2, 2))), silent = TRUE)
assert(inherits(invalid, "try-error"), "Zero-density references must be rejected")
message("Contour-threshold regression tests passed")

# Coarse backbone states are comparison labels, not validation categories.
states <- ram_backbone_basin(
  c(-63,-135,-75,60,0,NA),
  c(-43,135,145,40,0,0)
)
assert(identical(states,
  c("Alpha-R","Beta","PPII","Alpha-L","Other",NA_character_)),
  "Backbone-state prototypes or conservative fallback changed")
assert(identical(ram_backbone_basin(c(179,-179),c(179,-179)),
                 c("Beta","Beta")),
       "Backbone-state distance must remain periodic at the angle seam")
message("Backbone-state regression tests passed")
