# Run from the repository root: Rscript tests/scientific.R
source(file.path("shinyRam", "R", "ramachandran.R"))
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
