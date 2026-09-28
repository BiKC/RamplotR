# Normalize residue keys before comparing independent torsion results.
# Do not remove insertion codes: that would merge distinct PDB residues.
ram_match_reference_keys <- function(query_keys, reference_keys) {
  normalized <- trimws(reference_keys)
  if (anyDuplicated(normalized)) {
    stop("Duplicate Bio3D residue identifiers after whitespace normalization.")
  }
  match(trimws(query_keys), normalized)
}
