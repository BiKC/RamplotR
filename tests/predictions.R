# Pure-R predicted-structure confidence regression tests.
source(file.path("shinyRam", "R", "predictions.R"))
assert <- function(condition, message) if (!isTRUE(condition))
  stop(message, call. = FALSE)

torsions <- data.frame(
  chain = c("A", "A", "A"), resi = 1:3,
  insertion_code = c("", "", ""), resn = c("ALA", "GLY", "PRO"),
  region = c("Favoured", "Not allowed", "Favoured"),
  stringsAsFactors = FALSE
)
atom <- data.frame(
  chain = c("A", "A", "A", "A", "A", "A"),
  resno = c(1, 1, 2, 2, 3, 3), insert = "",
  elety = rep(c("N", "CA"), 3L), b = c(91, 93, 78, 77, 22, 21),
  stringsAsFactors = FALSE
)
pdb <- list(atom = atom)
experimental <- try(ram_prediction_from_atoms(pdb, torsions, "experimental"),
                    silent = TRUE)
assert(inherits(experimental, "try-error"),
       "Experimental crystallographic B factors must not be used as confidence")
for (source in c("alphafold_db", "alphafold2", "alphafold3", "esmfold")) {
  prediction <- ram_prediction_from_atoms(pdb, torsions, source)
  assert(identical(unname(prediction$plddt), c(93, 77, 21)),
         paste("Incorrect pLDDT for", source))
}
assert(identical(ram_plddt_category(c(NA, 32, 60, 85, 95)),
                 c("Unavailable", "Very low", "Low", "Confident", "Very high")),
       "Confidence categories or inclusive boundaries changed")

json <- list(predicted_aligned_error = list(
  list(0, 1, 2), list(2, 0, 4), list(3, 5, 0)
))
af2 <- ram_prediction_json(json, torsions, source = "alphafold2")
assert(identical(dim(af2$pae), c(3L, 3L)) &&
         identical(unname(af2$pae[2L, 1L]), 2),
       "AlphaFold DB PAE must retain alignment direction")
assert(identical(af2$pae_rows, 1:3), "AFDB PAE index mapping changed")

af3 <- list(
  pae = list(list(0, 4, 10, 11), list(3, 0, 8, 9),
             list(10, 8, 0, 2), list(11, 9, 2, 0)),
  token_chain_ids = as.list(c("A", "A", "X", "A")),
  token_res_ids = as.list(c(1, 2, 1, 3)),
  atom_chain_ids = as.list(rep("A", nrow(atom))),
  atom_plddts = as.list(c(91, 95, 78, 79, 82, 83)),
  ptm = 0.67
)
mapped <- ram_prediction_json(af3, torsions, atoms = atom,
                              source = "alphafold3")
assert(identical(dim(mapped$pae), c(3L, 3L)) &&
         identical(mapped$pae_rows, 1:3) &&
         isTRUE(all.equal(mapped$plddt, c(95, 79, 83))) &&
         mapped$ptm == 0.67,
       "AlphaFold 3 protein-token or per-atom mapping changed")
assert(mapped$pae[3, 1] == 11 && mapped$pae[1, 3] == 11,
       "AF3 PAE must retain token orientation and original numeric data")

altered <- af3
altered$atom_chain_ids[[1L]] <- "B"
mismatch <- ram_prediction_json(altered, torsions, atoms = atom,
                                source = "alphafold3")
assert(all(is.na(mismatch$plddt)) && length(mismatch$notes) > 0L,
       "Mismatched AF3 atom chains must never silently map confidence")

insertion <- torsions
insertion$insertion_code[2] <- "A"
inserted <- ram_prediction_json(af3, insertion, atoms = atom,
                                source = "alphafold3")
assert(is.null(inserted$pae), "Insertion-code ambiguity must disable AF3 PAE")
multimer <- torsions
multimer$chain[3] <- "B"
ambiguous <- ram_prediction_json(json, multimer, source = "alphafold2")
assert(is.null(ambiguous$pae) && length(ambiguous$notes) > 0L,
       "Unmapped multimer PAE must not be silently treated as sequential")

invalid <- try(ram_square_pae(list(c(0, 1), c(2))), silent = TRUE)
assert(inherits(invalid, "try-error"), "Ragged PAE must be rejected")
invalid <- try(ram_square_pae(list(c(0, -1), c(2, 0))), silent = TRUE)
assert(inherits(invalid, "try-error"), "Negative PAE must be rejected")
invalid <- try(ram_square_pae(rep(list(0), 1601)), silent = TRUE)
assert(inherits(invalid, "try-error"), "Unbounded PAE dimensions must be rejected")

prediction <- list(residues = ram_prediction_from_atoms(pdb, torsions, "esmfold"))
joined <- ram_apply_prediction(torsions, prediction)
assert(identical(unname(joined$plddt), c(93, 77, 21)) &&
         !anyDuplicated(ram_prediction_key(joined$chain, joined$resi,
                                           joined$insertion_code)),
       "Confidence join must preserve residue identifiers and row order")
assert(identical(ram_confidence_review(joined),
                 c("High confidence · geometry in range", "Other",
                   "Lower-confidence prediction")),
       "Review status must remain independent from scientific class")
missing_angle <- joined
missing_angle$region[[1L]] <- NA_character_
assert(identical(ram_confidence_review(missing_angle)[[1L]],
                 "Confidence available · backbone geometry unassessed"),
       "Missing phi/psi must never be described as acceptable geometry")

# Exercise the same API used by the Shiny upload workflow, not only its
# constituent parsers. The base ESMFold path requires no optional JSON package.
uploaded <- ram_prepare_prediction(pdb, torsions, "esmfold")
assert(identical(unname(uploaded$residues$plddt), c(93, 77, 21)) &&
         is.null(uploaded$pae) && identical(uploaded$source, "esmfold"),
       "Declared ESMFold structures must expose pLDDT without fabricating PAE")


mock_api <- function(url) {
  assert(grepl("https://alphafold.ebi.ac.uk/api/prediction/P12345",
               url, fixed = TRUE), "AFDB URL changed unexpectedly")
  list(list(modelEntityId = "AF-P12345-F1",
            pdbUrl = "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb",
            paeDocUrl = "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-predicted_aligned_error_v4.json"))
}
entry <- ram_afdb_entry("p12345", fetch = mock_api)
assert(identical(entry$accession, "P12345") &&
         identical(entry$structure_format, "pdb") &&
         !is.null(entry$pae_url),
       "AlphaFold DB model link parsing changed")
bad <- try(ram_afdb_entry("p12345", fetch = function(url)
  list(list(pdbUrl = "https://evil.example/pdb"))), silent = TRUE)
assert(inherits(bad, "try-error"), "Never fetch arbitrary model URLs")
# Even when visualising a large AF multimer at reduced resolution, show
# the residues immediately before and after each protein-chain boundary.
large_torsions <- data.frame(chain = c(rep("A", 250), rep("B", 260)),
  resi = c(seq_len(250), seq_len(260)), insertion_code = "",
  stringsAsFactors = FALSE)
large_prediction <- list(
  pae = matrix(2.5, 510, 510), pae_rows = seq_len(510)
)
overview <- ram_pae_plot_data(large_prediction, large_torsions,
                              max_display = 400L)
indices <- vapply(overview$residues, function(x) {
  if (x$chain == "A") x$resi else 250L + x$resi
}, integer(1))
assert(length(overview$labels) == 400L &&
         all(c(1L, 250L, 251L, 510L) %in% indices) &&
         overview$downsampled,
       "Downsampled PAE must retain all chain boundaries and termini")

message("AlphaFold and ESMFold prediction-format tests passed")
