# Run from repository root: Rscript tests/peptide.R
source(file.path("shinyRam", "R", "backbone.R"))
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
# Synthetic chains include an insertion code and a second disconnected chain.
make_atoms <- function(chain, number, insertion, resn, shift = 0) {
  data.frame(
    chain = chain, resno = number, insert = insertion,
    resid = resn, elety = c("N", "CA", "C"), alt = "",
    x = c(shift, shift + 0.6, shift + 1.3),
    y = c(0, 0.9, 0.7), z = c(0, 0.2, 1.0),
    stringsAsFactors = FALSE
  )
}
a <- rbind(
  make_atoms("A", 1, "", "ALA", 0),
  make_atoms("A", 1, "A", "PRO", 2.5),
  make_atoms("A", 2, "", "SER", 20),
  make_atoms("B", 1, "", "PRO", 4.9)
)
p <- ram_extract_torsions(list(atom = a))
assert(nrow(p) == 4L, "Insertion codes or duplicate chain numbering lost")
assert(identical(p$insertion_code, c("", "A", "", "")),
       "Insertion codes must be preserved")
assert(identical(p$bonded_to_next, c(TRUE, FALSE, FALSE, FALSE)),
       "Only geometric peptide neighbours in the same chain count")
assert(identical(p$next_resn[1L], "PRO"), "Pre-proline neighbour not detected")
assert(is.na(p$next_resn[3L]), "A broken peptide chain cannot have a neighbour")
assert(is.finite(p$psi[1L]) && is.finite(p$phi[2L]),
       "Connected complete residues should have phi/psi torsions")
assert(is.na(p$psi[2L]) && is.na(p$phi[3L]),
       "Disconnected residue torsions must be missing")
assert(is.na(ram_dihedral(c(0, 0, 0), c(1, 0, 0),
                         c(2, 0, 0), c(3, 0, 0))),
       "Collinear coordinates should not produce fabricated torsion angles")
message("Peptide continuity and insertion-code tests passed")

# Many short disconnected chains exercise the preallocated output path.
copies <- lapply(seq_len(250L), function(k) {
  make_atoms(sprintf("chain%03d", k), 1L, "", "ALA", shift = k * 100)
})
many <- ram_extract_torsions(list(atom = do.call(rbind, copies)))
assert(nrow(many) == 250L, "Every distinct chain must be retained")
assert(all(!many$bonded_to_next), "Output allocation must not cross chains")
assert(all(is.na(many$phi) & is.na(many$psi)),
       "Disconnected single-residue chains have undefined torsions")
message("Preallocated output regression tests passed")
