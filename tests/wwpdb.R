# Pure tests for independent wwPDB report ingestion; no Bio3D or web needed.
source("benchmarks/wwpdb_helpers.R")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
fixture <- "validation/fixtures/synthetic-validation.xml"
data <- ram_read_wwpdb_report(fixture)
assert(nrow(data) == 7L, "XML reader dropped residue records.")
assert(all(c("favored", "allowed", "outlier") %in% data$wwpdb_region),
       "wwPDB category normalization failed.")
assert(identical(ram_wwpdb_report_url("1CRN"),
  "https://files.rcsb.org/pub/pdb/validation_reports/cr/1crn/1crn_validation.xml.gz"),
  "Public wwPDB validation report path is incorrect.")
stopifnot(isTRUE(all.equal(
  ram_circular_angle_difference(c(179.9, -179.9), c(-179.8, 179.8)),
  c(0.3, 0.3), tolerance = 1e-10
)))
ours <- data.frame(
  chain = c("A", "A", "B", "B"), resi = c(1L, 2L, 10L, 11L),
  insertion_code = c("", "", "A", ""),
  resn = c("ALA", "GLY", "PRO", "SER"),
  phi = c(179.9, -70, -70, -50), psi = c(100, 80, 145, -40),
  region = c("Favoured", "Allowed", "Generously allowed", "Not allowed"),
  reference_group = c("General", "GLY", "PRO", "General"),
  stringsAsFactors = FALSE
)
comp <- ram_compare_wwpdb(ours, data)
assert(nrow(comp) == 4L && all(comp$identifier_matched),
       "Stable chain/residue/insertion joins failed.")
assert(isTRUE(all.equal(comp$phi_abs_degrees[1L], 0.3)),
       "Circular angle near +/-180 degrees was computed incorrectly.")
assert(all(comp$labels_agree), "Crosswalk should agree with this fixture.")
assert(comp$wwpdb_phi[2L] == -70,
       "Unlabelled primary conformer must precede alternate B.")
assert(comp$wwpdb_region[3L] == "allowed",
       "Residue insertion codes must remain distinguishable.")
summary <- ram_wwpdb_summary(comp, accession = "SYN1")
assert(summary$matched_phi_psi == 4L && summary$comparable_labels == 4L &&
         summary$label_agreement == 1,
       "Coverage or class summary incorrect.")
cross <- ram_wwpdb_contingency(comp)
assert(sum(cross$residues) == 4L &&
         sum(cross$ramplotr == cross$wwpdb & cross$residues > 0L) == 4L,
       "Contingency matrix silently changed categories.")
# Differential tests: real disagreement is counted but does not fail the
# scientific-angle gate. Gross angle differences or missing wwPDB data do.
different <- data
different$wwpdb_region[different$resi == 10L] <- "outlier"
changed <- ram_compare_wwpdb(ours, different)
assert(sum(!changed$labels_agree) == 1L,
       "Scientific classification disagreement must remain visible.")
assert(ram_wwpdb_summary(changed, "SYN1")$differing_labels == 1L,
       "Disagreements should be quantified, not used as an equality assertion.")
incorrect <- data
incorrect$phi[incorrect$resi == 1L & incorrect$model == "1"] <- -160
bad_angle <- try(ram_wwpdb_summary(
  ram_compare_wwpdb(ours, incorrect), "SYN1"), silent = TRUE)
assert(inherits(bad_angle, "try-error"),
       "Large independent angle differences must fail validation.")
missing <- data
missing$phi[missing$model == "1"] <- NA_real_
bad_coverage <- try(ram_wwpdb_summary(
  ram_compare_wwpdb(ours, missing), "SYN1"), silent = TRUE)
assert(inherits(bad_coverage, "try-error"),
       "Missing reference angles must fail the coverage gate.")
message("wwPDB XML parsing, identity, circular-angle and class tests passed.")
