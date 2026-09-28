# Parse structure input independently of the Shiny interface.
ram_detect_format <- function(filename) {
  if (!is.character(filename) || length(filename) != 1L ||
      is.na(filename) || !nzchar(filename)) {
    stop("Choose a .pdb, .ent, .cif or .mmcif structure file.")
  }
  extension <- tolower(tools::file_ext(filename))
  if (extension %in% c("pdb", "ent")) return("pdb")
  if (extension %in% c("cif", "mmcif", "mcif")) return("cif")
  stop("Unsupported structure format. Upload a PDB or mmCIF file.")
}

ram_load_structure <- function(path = NULL, original_name = NULL,
                               pdb_id = NULL, read_pdb = bio3d::read.pdb,
                               read_cif = bio3d::read.cif) {
  if (!is.null(path)) {
    if (!file.exists(path)) stop("The uploaded structure is unavailable.")
    format <- ram_detect_format(original_name)
    if (format == "pdb") {
      structure <- read_pdb(path, rm.insert = FALSE, rm.alt = FALSE,
                            ATOM.only = TRUE, verbose = FALSE)
    } else {
      structure <- read_cif(path, rm.insert = FALSE, rm.alt = FALSE,
                            verbose = FALSE)
    }
  } else {
    if (!is.character(pdb_id) || length(pdb_id) != 1L ||
        is.na(pdb_id) || !grepl("^[[:alnum:]]{4}$", pdb_id)) {
      stop("Enter a valid four-character PDB accession.")
    }
    structure <- tryCatch(
      read_cif(pdb_id, rm.insert = FALSE, rm.alt = FALSE, verbose = FALSE),
      error = function(cif_error) {
        tryCatch(
          read_pdb(pdb_id, rm.insert = FALSE, rm.alt = FALSE,
                   ATOM.only = TRUE, verbose = FALSE),
          error = function(pdb_error) {
            stop(sprintf("Could not retrieve %s: %s",
                         pdb_id, conditionMessage(pdb_error)), call. = FALSE)
          }
        )
      }
    )
  }
  if (!is.list(structure) || !is.data.frame(structure$atom) ||
      nrow(structure$atom) == 0L) {
    stop("This structure contains no atom coordinates.")
  }
  structure
}
