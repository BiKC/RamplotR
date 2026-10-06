# Usage from repository root:
# Rscript benchmarks/compare_wwpdb.R 1CRN original benchmarks/output/wwpdb/1CRN \
#   benchmarks/output/wwpdb/input/1crn_validation.xml.gz \
#   benchmarks/output/wwpdb/input/1CRN.cif
#
# The final two input paths are optional. Supplying both pins the exact
# coordinates and official report and avoids silently changing source data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 5L) {
  stop("Usage: Rscript benchmarks/compare_wwpdb.R PDB_ID REFERENCE_SET OUTPUT_DIR [VALIDATION_XML_GZ] [PDB_CIF]")
}
accession <- toupper(args[[1L]])
reference_set <- args[[2L]]
output <- args[[3L]]
if (!grepl("^[A-Z0-9]{4}$", accession))
  stop("Expected a four-character PDB accession.")
known <- c("original", "alphafold", "alphafold_filtered",
           "astral2.08", "custom_high_resolution")
if (!reference_set %in% known) stop("Unrecognized reference dataset.")
source("shinyRam/R/io.R")
source("shinyRam/R/backbone.R")
source("shinyRam/R/ramachandran.R")
source("shinyRam/R/rama8000.R")
source("benchmarks/wwpdb_helpers.R")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

xml_path <- if (length(args) >= 4L) args[[4L]] else
  file.path(output, paste0(tolower(accession), "_validation.xml.gz"))
cif_path <- if (length(args) >= 5L) args[[5L]] else
  file.path(output, paste0(accession, ".cif"))
validation_url <- ram_wwpdb_report_url(accession)
coordinate_url <- sprintf("https://files.rcsb.org/download/%s.cif", accession)
if (!file.exists(xml_path))
  utils::download.file(validation_url, xml_path, mode = "wb", quiet = FALSE)
if (!file.exists(cif_path))
  utils::download.file(coordinate_url, cif_path, mode = "wb", quiet = FALSE)
if (file.info(xml_path)$size < 100L || file.info(cif_path)$size < 100L)
  stop("Downloaded validation report or coordinate file is unexpectedly small.")

# Explicitly select model 1 for multi-model entries. Later models require
# model-matched wwPDB reports and a new per-model comparison.
pdb <- ram_model_at(
  ram_load_structure(path = cif_path, original_name = paste0(accession, ".cif")),
  1L
)
torsions <- ram_extract_torsions(pdb)
if (!nrow(torsions)) stop("No protein backbone coordinates were extracted.")
refdir <- file.path("shinyRam", "static", reference_set)
classified <- ram_classify_torsions(
  torsions, refdir, NULL, "residue", threshold_fn = ram_density_thresholds
)
classified <- ram_rama8000_classify(
  classified, file.path("shinyRam", "static", "rama8000")
)
external <- ram_read_wwpdb_report(xml_path)
details <- ram_compare_wwpdb(classified, external, model = 1L)
standard <- ram_compare_rama8000_wwpdb(classified, external, model = 1L)
detail_path <- file.path(output, "residue_comparison.csv")
write.csv(details, detail_path, row.names = FALSE, na = "")
cross <- ram_wwpdb_contingency(details)
write.csv(cross, file.path(output, "group_contingency.csv"), row.names = FALSE)
write.csv(standard, file.path(output, "rama8000_comparison.csv"),
          row.names = FALSE, na = "")
write.csv(ram_rama8000_wwpdb_contingency(standard),
          file.path(output, "rama8000_contingency.csv"), row.names = FALSE)

# Persist source hashes and the full comparison even if a data/angle quality
# gate fails; a disagreement in region LABELS is never a CI-failure criterion.
source_info <- data.frame(
  accession = accession, reference_set = reference_set,
  coordinates = coordinate_url, validation_xml = validation_url,
  coordinate_sha256 = digest::digest(file = cif_path, algo = "sha256"),
  validation_sha256 = digest::digest(file = xml_path, algo = "sha256"),
  extraction_model = 1L,
  regions_ramplotr = "Favoured|Allowed|Generously allowed|Not allowed",
  regions_wwpdb = "Favored|Allowed|OUTLIER",
  ramplotr_group_count = 4L, wwpdb_group_count = 6L,
  analysis_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  repository_commit = tryCatch(
    trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[1L]),
    error = function(e) "unavailable"),
  stringsAsFactors = FALSE
)
write.csv(source_info, file.path(output, "provenance.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()),
             paste("xml2", as.character(utils::packageVersion("xml2"))),
             paste("digest", as.character(utils::packageVersion("digest")))),
           file.path(output, "r-session.txt"))
summary <- ram_wwpdb_summary(details, accession,
                             min_coverage = 0, max_angle_deviation = 180)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
standard_summary <- ram_rama8000_wwpdb_summary(
  standard, accession, min_coverage = 0, min_agreement = 0)
write.csv(standard_summary, file.path(output, "rama8000-summary.csv"),
          row.names = FALSE)
print(summary, row.names = FALSE)
print(standard_summary, row.names = FALSE)
print(cross[cross$residues > 0L, , drop = FALSE], row.names = FALSE)
print(ram_rama8000_wwpdb_contingency(standard), row.names = FALSE)
# Independent geometry agreement must be established before making
# interpretive claims from the classification comparison.
gate <- tryCatch(
  ram_wwpdb_summary(details, accession, min_coverage = 0.90,
                    max_angle_deviation = 1.5),
  error = function(e) e
)
if (inherits(gate, "error")) {
  writeLines(conditionMessage(gate), file.path(output, "GATE_FAILURE.txt"))
  stop(conditionMessage(gate), call. = FALSE)
}
standard_gate <- tryCatch(
  ram_rama8000_wwpdb_summary(
    standard, accession, min_coverage = 0.90, min_agreement = 0.99),
  error = function(e) e
)
if (inherits(standard_gate, "error")) {
  writeLines(conditionMessage(standard_gate),
             file.path(output, "RAMA8000_GATE_FAILURE.txt"))
  stop(conditionMessage(standard_gate), call. = FALSE)
}
message("wwPDB angle gate and direct Rama8000 category gate passed.")
