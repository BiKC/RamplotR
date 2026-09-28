# Independent angle comparison: Rscript benchmarks/validate_vs_bio3d.R 1CRN
# Requires bio3d and internet access to retrieve the requested PDB entry.
args <- commandArgs(trailingOnly = TRUE)
accession <- if (length(args) >= 1L) toupper(trimws(args[[1L]])) else "1CRN"
output <- if (length(args) >= 2L) args[[2L]] else
  file.path("benchmarks", "validation-results.csv")
source(file.path("shinyRam", "R", "io.R"))
source(file.path("shinyRam", "R", "backbone.R"))

pdb <- ram_load_structure(pdb_id = accession)
ours <- ram_extract_torsions(pdb)
reference <- bio3d::torsion.pdb(pdb)$tbl
if (!nrow(ours) || !nrow(reference)) stop("No comparable residues")

# Bio3D's table uses row names of the form residue.chain.name[.insertion].
# Compare only rows with no insertion code so that atom numbering cannot
# ambiguously match two distinct residues in this first cross-check.
plain <- !is.na(ours$insertion_code) & ours$insertion_code == ""
ours <- ours[plain, , drop = FALSE]
key <- paste(ours$resi, ours$chain, ours$resn, sep = ".")
reference_index <- match(key, rownames(reference))
angle_difference <- function(a, b) abs(((a - b + 180) %% 360) - 180)

results <- lapply(c("phi", "psi"), function(angle) {
  i <- which(!is.na(reference_index) & is.finite(ours[[angle]]))
  if (length(i)) {
    i <- i[is.finite(reference[reference_index[i], angle])]
  }
  if (!length(i)) stop(sprintf("No comparable %s angles in %s", angle, accession))
  differences <- angle_difference(ours[[angle]][i],
                                  as.numeric(reference[reference_index[i], angle]))
  data.frame(
    accession = accession, angle = angle, matched = length(i),
    max_abs_degrees = max(differences),
    median_abs_degrees = median(differences),
    mismatches_over_0.5_degrees = sum(differences > 0.5),
    stringsAsFactors = FALSE
  )
})
report <- do.call(rbind, results)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(report, output, row.names = FALSE)
writeLines(c(
  paste("accession:", accession),
  paste("timestamp_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  capture.output(sessionInfo())
), paste0(output, ".session.txt"))
print(report, row.names = FALSE)
if (any(report$mismatches_over_0.5_degrees != 0)) {
  stop("Independent angle comparison exceeded the 0.5-degree tolerance")
}
message("Bio3D torsion comparison passed")
