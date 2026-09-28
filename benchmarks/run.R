# Reproducible benchmark:
#   Rscript benchmarks/run.R 1CRN 3 original benchmarks/results.csv
#   Rscript benchmarks/run.R local_complex.cif 3 original benchmarks/results.csv
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[[1L]] else "1CRN"
repeats <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3L
reference_set <- if (length(args) >= 3L) args[[3L]] else "original"
output <- if (length(args) >= 4L) args[[4L]] else
  file.path("benchmarks", "results.csv")
if (!is.finite(repeats) || repeats < 1L || repeats > 100L) {
  stop("Repetitions must be between 1 and 100.")
}
if (!reference_set %in% c("original", "alphafold", "alphafold_filtered",
                          "astral2.08", "custom_high_resolution")) {
  stop("Unknown bundled reference dataset.")
}
source(file.path("shinyRam", "R", "io.R"))
source(file.path("shinyRam", "R", "backbone.R"))
source(file.path("shinyRam", "R", "ramachandran.R"))
refdir <- file.path("shinyRam", "static", reference_set)

started <- proc.time()[["elapsed"]]
structure <- if (file.exists(input)) {
  ram_load_structure(path = input, original_name = basename(input))
} else {
  ram_load_structure(pdb_id = toupper(input))
}
parse_seconds <- proc.time()[["elapsed"]] - started
general <- ram_read_reference(file.path(refdir, "General"))
# Warm up the immutable reference cache separately from timed calculations.
invisible(lapply(c("General", "GLY", "PRO", "preProline"),
                 function(name) ram_read_reference(file.path(refdir, name))))
record <- vector("list", repeats)
for (iteration in seq_len(repeats)) {
  gc()
  tic <- proc.time()[["elapsed"]]
  torsions <- ram_extract_torsions(structure)
  torsion_seconds <- proc.time()[["elapsed"]] - tic

  tic <- proc.time()[["elapsed"]]
  classified <- ram_classify_torsions(
    torsions, refdir, general, mode = "residue",
    threshold_fn = ram_density_thresholds
  )
  classification_seconds <- proc.time()[["elapsed"]] - tic
  record[[iteration]] <- data.frame(
    input = input, reference_set = reference_set, iteration = iteration,
    atoms = nrow(structure$atom), residues = nrow(torsions),
    classified = sum(!is.na(classified$region)),
    parser_seconds = parse_seconds,
    torsions_seconds = torsion_seconds,
    classification_seconds = classification_seconds,
    total_processing_seconds = torsion_seconds + classification_seconds,
    pdb_object_bytes = as.numeric(object.size(structure))
  )
}
results <- do.call(rbind, record)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
writeLines(c(
  paste("input:", input), paste("reference_set:", reference_set),
  paste("timestamp_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  capture.output(sessionInfo())
), paste0(output, ".session.txt"))
print(results, row.names = FALSE)
