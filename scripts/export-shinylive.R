# Usage: Rscript scripts/export-shinylive.R bikc.be https://bikc.be/RamplotR/reference-data
# Run from repository root. Source scientific files are never modified.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L ||
    !grepl("^https?://", args[[2L]])) {
  stop("Usage: Rscript scripts/export-shinylive.R OUTPUT_DIR ABSOLUTE_REFERENCE_URL")
}
if (!requireNamespace("shinylive", quietly = TRUE)) {
  stop("Install the shinylive R package")
}
target <- normalizePath(args[[1L]], mustWork = FALSE)
base_url <- sub("/+$", "", args[[2L]])
source_app <- normalizePath("shinyRam", mustWork = TRUE)
reference_dirs <- c("original", "alphafold", "alphafold_filtered",
                    "astral2.08", "custom_high_resolution")
stage <- tempfile("ramplotr-thin-")
dir.create(stage)
on.exit(unlink(stage, recursive = TRUE), add = TRUE)
thin <- file.path(stage, "shinyRam")
dir.create(thin)
dir.create(file.path(thin, "static"))
# Copy all app logic and web assets, but not the large reference datasets.
for (entry in c("app.R", "R", "www")) {
  if (!file.copy(file.path(source_app, entry), thin, recursive = TRUE)) {
    stop("Failed to copy ", entry)
  }
}
# Preserve only the small manifest inside the app; the original RDS bytes
# remain unchanged and are published alongside the static website.
public <- file.path(target, "RamplotR", "reference-data")
dir.create(public, recursive = TRUE, showWarnings = FALSE)
for (dataset in reference_dirs) {
  source_dir <- file.path(source_app, "static", dataset)
  names <- list.files(source_dir, all.files = FALSE)
  names <- names[file.info(file.path(source_dir, names))$isdir %in% FALSE]
  if (!length(names)) stop("Missing reference dataset: ", dataset)
  from <- file.path(source_dir, names)
  manifest <- data.frame(file = names, md5 = unname(tools::md5sum(from)))
  index_dir <- file.path(thin, "static", dataset)
  dir.create(index_dir)
  utils::write.table(manifest, file.path(index_dir, "reference-index.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  public_dir <- file.path(public, dataset)
  dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, public_dir, overwrite = TRUE)
  if (!all(ok)) stop("Failed to publish reference files: ", dataset)
  published <- file.path(public_dir, names)
  if (!identical(unname(tools::md5sum(published)), manifest$md5)) {
    stop("Published references differ: ", dataset)
  }
}
# Set the public base URL inside the staged app only.
app <- file.path(thin, "app.R")
lines <- readLines(app, warn = FALSE)
lines <- append(lines, paste0("options(ramplotr.reference_base_url = ",
                             deparse(base_url), ")"), after = 0L)
writeLines(lines, app)
# Export into a fresh directory so old packed RDS files cannot survive.
# Keep published reference-data aside while removing an earlier app export.
tmp_public <- tempfile("ramplotr-reference-data-")
dir.create(tmp_public)
if (!file.copy(public, tmp_public, recursive = TRUE)) {
  stop("Cannot preserve reference-data during rebuild")
}
unlink(file.path(target, "RamplotR"), recursive = TRUE)
shinylive::export(thin, target, subdir = "RamplotR")
dir.create(file.path(target, "RamplotR"), recursive = TRUE,
           showWarnings = FALSE)
restored <- file.path(tmp_public, "reference-data")
if (!file.copy(restored, file.path(target, "RamplotR"), recursive = TRUE)) {
  stop("Failed to restore reference-data after export")
}
message("Thin Shinylive export complete: ", file.path(target, "RamplotR"))
message("Publish the full output directory so reference-data/ is accessible at ",
        base_url)
