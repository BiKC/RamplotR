# Exact-source reproducibility gate for the independent experimental cohort.
# Usage: Rscript tests/wwpdb-baseline.R ACCESSION OUTPUT_DIR
# This is run after benchmarks/compare_wwpdb.R has saved its provenance,
# per-residue comparison and scientific summary.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("Usage: Rscript tests/wwpdb-baseline.R ACCESSION OUTPUT_DIR")
id <- toupper(args[[1L]])
folder <- args[[2L]]
expect <- read.csv("validation/baseline-2026-09-28.csv",
                   stringsAsFactors = FALSE)
hashes <- read.csv("validation/baseline-source-hashes-2026-09-28.csv",
                   stringsAsFactors = FALSE)
expect <- expect[expect$accession == id, , drop = FALSE]
hashes <- hashes[hashes$accession == id, , drop = FALSE]
if (nrow(expect) != 1L || nrow(hashes) != 1L)
  stop("No unique pinned baseline for the requested accession.")
actual <- read.csv(file.path(folder, "summary.csv"), stringsAsFactors = FALSE)
sources <- read.csv(file.path(folder, "provenance.csv"),
                    stringsAsFactors = FALSE)
details <- read.csv(file.path(folder, "residue_comparison.csv"),
                    stringsAsFactors = FALSE)
if (nrow(actual) != 1L || nrow(sources) != 1L)
  stop("No unique measured summary and source provenance.")
if (sources$coordinate_sha256 != hashes$coordinate_sha256 ||
    sources$validation_sha256 != hashes$validation_sha256) {
  stop(paste("Independent source files for", id, "changed since the",
             "pinned baseline. Inspect/revalidate their exact contents and",
             "update the versioned source manifest in a reviewed commit."))
}
numeric_fields <- c(
  "ram_residues", "ram_finite_phi_psi", "wwpdb_identifier_matches",
  "matched_phi_psi", "matched_angle_coverage",
  "max_phi_difference", "max_psi_difference", "comparable_labels",
  "matching_labels", "differing_labels", "label_agreement"
)
for (field in numeric_fields) {
  if (!isTRUE(all.equal(as.numeric(actual[[field]]),
                        as.numeric(expect[[field]]), tolerance = 1e-6))) {
    stop(sprintf("%s no longer matches the pinned %s baseline.", field, id))
  }
}
official_outliers <- sum(details$wwpdb_region == "outlier", na.rm = TRUE)
our_outliers <- sum(details$ram_region_three == "outlier", na.rm = TRUE)
if (official_outliers != expect$wwpdb_outliers ||
    our_outliers != expect$ramplotr_outliers_3way) {
  stop("Independent outlier counts changed since the pinned reference snapshot.")
}
if (id == "2DQ4") {
  pinned <- read.csv("validation/2dq4-wwpdb-outlier-examples.csv",
                     stringsAsFactors = FALSE)
  external <- details[details$wwpdb_region == "outlier", , drop = FALSE]
  key <- function(chain, num, ins, name) {
    paste(chain, num, ifelse(is.na(ins), "", ins), name, sep = ":")
  }
  selection <- match(key(pinned$chain, pinned$residue_number,
                         pinned$insertion_code, pinned$residue_name),
                     key(external$chain, external$resi,
                         external$insertion_code, external$resn))
  if (anyNA(selection) || length(selection) != nrow(external))
    stop("Pinned independent 2DQ4 outlier identities are no longer identical.")
  selected <- external[selection, , drop = FALSE]
  if (!identical(as.character(selected$ram_region_four),
                 as.character(pinned$ramplotr_region)) ||
      !identical(as.character(selected$wwpdb_region),
                 as.character(pinned$wwpdb_region))) {
    stop("2DQ4 residue classifications changed relative to the pinned fixture.")
  }
}
message("Pinned independent snapshot reproduced for ", id,
        " (input SHA256 hashes, angles, labels and outlier counts).")
