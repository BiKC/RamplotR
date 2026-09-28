# Run from repository root: Rscript tests/input.R
source(file.path("shinyRam", "R", "io.R"))
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
assert(identical(ram_detect_format("protein.PDB"), "pdb"), "PDB format")
assert(identical(ram_detect_format("model.ent"), "pdb"), "ENT format")
assert(identical(ram_detect_format("protein.CIF"), "cif"), "CIF format")
assert(identical(ram_detect_format("protein.mmcif"), "cif"), "mmCIF format")
assert(inherits(try(ram_detect_format("protein.xyz"), silent = TRUE), "try-error"),
       "Invalid formats must fail")

uploaded <- tempfile()
file.create(uploaded)
mock_pdb <- function(path, ...) list(atom = data.frame(chain = "A"))
mock_cif <- function(path, ...) list(atom = data.frame(chain = "B"))
assert(identical(ram_load_structure(uploaded, "protein.pdb",
  read_pdb = mock_pdb, read_cif = mock_cif)$atom$chain, "A"), "PDB upload")
assert(identical(ram_load_structure(uploaded, "protein.mmcif",
  read_pdb = mock_pdb, read_cif = mock_cif)$atom$chain, "B"), "CIF upload")
assert(inherits(try(ram_load_structure(uploaded, "bad.xyz",
  read_pdb = mock_pdb, read_cif = mock_cif), silent = TRUE), "try-error"),
  "Unknown extension must fail before parsing")
assert(inherits(try(ram_load_structure(pdb_id = "../bad",
  read_pdb = mock_pdb, read_cif = mock_cif), silent = TRUE), "try-error"),
  "Invalid identifier must fail")
assert(identical(ram_load_structure(pdb_id = "1BBB",
  read_pdb = mock_pdb, read_cif = mock_cif)$atom$chain, "B"), "CIF first")
cif_error <- function(...) stop("CIF unavailable")
assert(identical(ram_load_structure(pdb_id = "1BBB",
  read_pdb = mock_pdb, read_cif = cif_error)$atom$chain, "A"),
  "Fallback to PDB")
unlink(uploaded)
message("Structure input tests passed")
