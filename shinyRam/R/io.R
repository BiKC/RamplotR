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
      structure <- read_pdb(path, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                            ATOM.only = TRUE, verbose = FALSE)
    } else {
      structure <- read_cif(path, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                            verbose = FALSE)
    }
  } else {
    if (!is.character(pdb_id) || length(pdb_id) != 1L ||
        is.na(pdb_id) || !grepl("^[[:alnum:]]{4}$", pdb_id)) {
      stop("Enter a valid four-character PDB accession.")
    }
    structure <- tryCatch(
      read_cif(pdb_id, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE, verbose = FALSE),
      error = function(cif_error) {
        tryCatch(
          read_pdb(pdb_id, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
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


# Bio3D's multi=TRUE keeps consistent atom records in $atom, with one
# coordinate row per model in $xyz. Never claim multi-model support when atom
# counts differ or when Bio3D cannot return a complete coordinate matrix.
ram_model_count <- function(pdb) {
  xyz <- pdb$xyz
  if (is.matrix(xyz) && ncol(xyz) == 3L * nrow(pdb$atom))
    max(1L, nrow(xyz))
  else 1L
}

ram_model_at <- function(pdb, index = 1L) {
  index <- suppressWarnings(as.integer(index))
  count <- ram_model_count(pdb)
  if (length(index) != 1L || is.na(index) || index < 1L || index > count)
    stop("The requested structural model is unavailable.")
  if (count == 1L) return(pdb)
  coords <- matrix(as.numeric(pdb$xyz[index, ]), ncol = 3L, byrow = TRUE)
  if (nrow(coords) != nrow(pdb$atom))
    stop("Model atom count differs from the reference model.")
  selected <- pdb
  selected$atom$x <- coords[, 1L]
  selected$atom$y <- coords[, 2L]
  selected$atom$z <- coords[, 3L]
  # Prevent downstream consumers from accidentally treating the selected
  # model as a new multi-model ensemble.
  selected$xyz <- matrix(as.vector(t(coords)), nrow = 1L)
  selected
}
