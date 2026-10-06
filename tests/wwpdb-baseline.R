# Exact-source reproducibility gate for the independent experimental cohort.
# Usage: Rscript tests/wwpdb-baseline.R ACCESSION OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("Usage: Rscript tests/wwpdb-baseline.R ACCESSION OUTPUT_DIR")

id <- toupper(args[[1L]])
folder <- args[[2L]]
geometry_expect <- read.csv("validation/baseline-2026-09-28.csv",
                            stringsAsFactors = FALSE)
standard_expect <- read.csv("validation/rama8000-baseline-2026-10-06.csv",
                            stringsAsFactors = FALSE)
hashes <- read.csv("validation/baseline-source-hashes-2026-09-28.csv",
                   stringsAsFactors = FALSE)
geometry_expect <- geometry_expect[geometry_expect$accession == id, , drop=FALSE]
standard_expect <- standard_expect[standard_expect$accession == id, , drop=FALSE]
hashes <- hashes[hashes$accession == id, , drop=FALSE]
if (nrow(geometry_expect) != 1L || nrow(standard_expect) != 1L ||
    nrow(hashes) != 1L)
  stop("No unique pinned baseline for the requested accession.")

actual <- read.csv(file.path(folder, "summary.csv"), stringsAsFactors = FALSE)
standard <- read.csv(file.path(folder, "rama8000-summary.csv"),
                     stringsAsFactors = FALSE)
sources <- read.csv(file.path(folder, "provenance.csv"),
                    stringsAsFactors = FALSE)
details <- read.csv(file.path(folder, "rama8000_comparison.csv"),
                    stringsAsFactors = FALSE)
if (nrow(actual) != 1L || nrow(standard) != 1L || nrow(sources) != 1L)
  stop("No unique measured summary and source provenance.")

if (sources$coordinate_sha256 != hashes$coordinate_sha256 ||
    sources$validation_sha256 != hashes$validation_sha256) {
  stop(paste("Independent source files for", id, "changed since the",
             "pinned baseline. Inspect/revalidate their exact contents and",
             "update the versioned source manifest in a reviewed commit."))
}

# Geometry correctness is independent of any classification scheme.
geometry_fields <- c(
  "ram_residues", "ram_finite_phi_psi", "wwpdb_identifier_matches",
  "matched_phi_psi", "matched_angle_coverage",
  "max_phi_difference", "max_psi_difference"
)
for (field in geometry_fields) {
  if (!isTRUE(all.equal(as.numeric(actual[[field]]),
                        as.numeric(geometry_expect[[field]]),
                        tolerance = 1e-6))) {
    stop(sprintf("%s no longer matches the pinned %s geometry baseline.",
                 field, id))
  }
}

# The current Rama8000 implementation is expected to reproduce official
# wwPDB categories exactly for the pinned corpus. This is a direct comparison,
# not a remapping of RamplotR's native four density regions.
standard_fields <- c(
  "finite_phi_psi", "comparable_rama8000", "comparison_coverage",
  "matching_rama8000", "differing_rama8000", "rama8000_agreement"
)
for (field in standard_fields) {
  if (!isTRUE(all.equal(as.numeric(standard[[field]]),
                        as.numeric(standard_expect[[field]]),
                        tolerance = 1e-12))) {
    stop(sprintf("%s no longer matches the pinned %s Rama8000 baseline.",
                 field, id))
  }
}
if (standard$rama8000_agreement[[1L]] != 1 ||
    standard$differing_rama8000[[1L]] != 0L)
  stop("Pinned Rama8000 categories no longer reproduce wwPDB exactly.")

official_outliers <- sum(details$wwpdb_region == "outlier", na.rm=TRUE)
standard_outliers <- sum(details$rama8000_region == "outlier", na.rm=TRUE)
if (official_outliers != geometry_expect$wwpdb_outliers ||
    standard_outliers != geometry_expect$wwpdb_outliers)
  stop("Official and Rama8000 outlier counts differ from the pinned snapshot.")

if (id == "2DQ4") {
  pinned <- read.csv("validation/2dq4-wwpdb-outlier-examples.csv",
                     stringsAsFactors = FALSE)
  external <- details[details$wwpdb_region == "outlier", , drop=FALSE]
  key <- function(chain, num, ins, name) {
    paste(chain, num, ifelse(is.na(ins), "", ins), name, sep=":")
  }
  selection <- match(
    key(pinned$chain, pinned$residue_number,
        pinned$insertion_code, pinned$residue_name),
    key(external$chain, external$resi,
        external$insertion_code, external$resn)
  )
  if (anyNA(selection) || length(selection) != nrow(external))
    stop("Pinned independent 2DQ4 outlier identities are no longer identical.")
  selected <- external[selection, , drop=FALSE]
  if (!identical(as.character(selected$ram_region_four),
                 as.character(pinned$ramplotr_region)) ||
      !identical(as.character(selected$rama8000_group),
                 as.character(pinned$rama8000_group)) ||
      !identical(as.character(selected$rama8000_region),
                 as.character(pinned$rama8000_region)) ||
      !identical(as.character(selected$wwpdb_region),
                 as.character(pinned$wwpdb_region))) {
    stop("2DQ4 direct standard/native classifications changed.")
  }
}

message("Pinned independent snapshot reproduced for ", id,
        " (source hashes, angles and exact Rama8000/wwPDB categories).")
