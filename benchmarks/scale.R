# Measure real-structure and replicated-complex scaling in fresh R processes.
#
# From the repository root:
#   Rscript benchmarks/scale.R benchmarks/output/source-6vxx.rds 3 3 original benchmarks/output/scaling/scaling-3.csv
#
# The multiplier duplicates the same experimental structure with unique chain
# names. It measures computational scaling, not independent biological data.
# Invoke each multiplier in a separate /usr/bin/time process to collect RSS.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript benchmarks/scale.R structure.rds multiplier repeats reference_set output.csv")
}
input <- args[[1L]]
multiplier <- suppressWarnings(as.integer(args[[2L]]))
repeats <- suppressWarnings(as.integer(args[[3L]]))
reference_set <- args[[4L]]
output <- args[[5L]]
if (length(multiplier) != 1L || is.na(multiplier) ||
    multiplier < 1L || multiplier > 20L) {
  stop("Multiplier must be an integer between 1 and 20.")
}
if (length(repeats) != 1L || is.na(repeats) ||
    repeats < 1L || repeats > 20L) {
  stop("Repeats must be an integer between 1 and 20.")
}
allowed_references <- c(
  "original", "alphafold", "alphafold_filtered",
  "astral2.08", "custom_high_resolution"
)
if (!reference_set %in% allowed_references) stop("Unknown reference dataset.")
if (!file.exists(input) || !grepl("\\.rds$", input, ignore.case = TRUE)) {
  stop("Provide an existing .rds file created from a single loaded structure.")
}
source(file.path("shinyRam", "R", "backbone.R"))
source(file.path("shinyRam", "R", "ramachandran.R"))

loaded <- readRDS(input)
if (!is.list(loaded) || !is.data.frame(loaded$atom) || !nrow(loaded$atom)) {
  stop("The input does not contain a nonempty Bio3D atom table.")
}
original_atoms <- loaded$atom
if (!"chain" %in% names(original_atoms)) {
  stop("The atom table has no chain column.")
}
original_chains <- as.character(original_atoms$chain)
original_chains[is.na(original_chains)] <- ""
# Prefix chain labels so repeated chains cannot acquire peptide bonds across
# copies even if residue numbering or coordinates are identical.
copies <- lapply(seq_len(multiplier), function(k) {
  atoms <- original_atoms
  atoms$chain <- paste0("copy", sprintf("%03d", k), "_", original_chains)
  atoms
})
loaded$atom <- do.call(rbind, copies)
rm(copies, original_atoms)
gc()

refdir <- file.path("shinyRam", "static", reference_set)
general <- ram_read_reference(file.path(refdir, "General"))
results <- vector("list", repeats)
last_torsions <- NULL
last_classified <- NULL
for (iteration in seq_len(repeats)) {
  gc()
  start <- proc.time()[["elapsed"]]
  torsions <- ram_extract_torsions(loaded)
  torsion_seconds <- proc.time()[["elapsed"]] - start

  start <- proc.time()[["elapsed"]]
  classified <- ram_classify_torsions(
    torsions, refdir, general,
    mode = "residue", threshold_fn = ram_density_thresholds
  )
  classification_seconds <- proc.time()[["elapsed"]] - start

  results[[iteration]] <- data.frame(
    structure_rds = basename(input),
    reference_set = reference_set,
    multiplier = multiplier,
    iteration = iteration,
    cache_state = if (iteration == 1L) "cold_profile" else "warm_profile",
    atoms = nrow(loaded$atom),
    residues = nrow(torsions),
    classified = sum(!is.na(classified$region)),
    torsions_seconds = torsion_seconds,
    classification_seconds = classification_seconds,
    total_processing_seconds = torsion_seconds + classification_seconds,
    pdb_object_bytes = as.numeric(object.size(loaded)),
    git_sha = Sys.getenv("GITHUB_SHA", unset = "local"),
    stringsAsFactors = FALSE
  )
  last_torsions <- torsions
  last_classified <- classified
}

# Scientific consistency check is outside timed regions. Every replicated
# copy must yield the same phi, psi, region and density as the first copy.
if (nrow(last_torsions) %% multiplier != 0L) {
  stop("Unequal residue counts after structural replication.")
}
per_copy <- nrow(last_torsions) %/% multiplier
if (per_copy == 0L) stop("No residues available for comparison.")
first <- seq_len(per_copy)
if (multiplier > 1L) {
  for (k in 2L:multiplier) {
    rows <- ((k - 1L) * per_copy + 1L):(k * per_copy)
    if (!identical(last_torsions$resn[first], last_torsions$resn[rows]) ||
        !isTRUE(all.equal(last_torsions$phi[first], last_torsions$phi[rows])) ||
        !isTRUE(all.equal(last_torsions$psi[first], last_torsions$psi[rows])) ||
        !identical(last_classified$region[first], last_classified$region[rows]) ||
        !isTRUE(all.equal(last_classified$density[first],
                          last_classified$density[rows]))) {
      stop(sprintf("Replicated copy %d differs from the first copy.", k))
    }
  }
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), output, row.names = FALSE)
writeLines(c(
  paste("structure_rds:", normalizePath(input, mustWork = TRUE)),
  paste("multiplier:", multiplier),
  paste("reference_set:", reference_set),
  paste("timestamp_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Peak RSS must be measured from an external process monitor.",
  capture.output(sessionInfo())
), paste0(output, ".session.txt"))
print(do.call(rbind, results), row.names = FALSE)
message("All replicated chains produced identical torsions and classifications.")
