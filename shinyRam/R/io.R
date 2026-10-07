# Parse structure input independently of the Shiny interface.
ram_nonwater_hetero_atoms <- function(atoms) {
  if (!is.data.frame(atoms) || !nrow(atoms) || !"type" %in% names(atoms))
    return(if (is.data.frame(atoms)) atoms[0,,drop=FALSE] else data.frame())
  type <- toupper(trimws(as.character(atoms$type)))
  resid <- if ("resid" %in% names(atoms))
    toupper(trimws(as.character(atoms$resid))) else rep("",nrow(atoms))
  waters <- c("HOH","WAT","DOD","H2O","SOL","TIP","TIP3","TIP3P")
  keep <- type %in% c("HETATM","HET") & !resid %in% waters
  if (all(c("x","y","z") %in% names(atoms)))
    keep <- keep & is.finite(atoms$x) & is.finite(atoms$y) & is.finite(atoms$z)
  atoms[keep,,drop=FALSE]
}

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
      # Keep the analysis atom table protein-only, but retain non-water HETATM
      # records separately for residue-level spatial context. This avoids
      # changing backbone/model bookkeeping while making ligands/cofactors/ions
      # available to the inspector.
      full_structure <- tryCatch(
        read_pdb(path, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                 ATOM.only = FALSE, verbose = FALSE),
        error = function(e) NULL
      )
      structure$hetero_atom <- if (is.null(full_structure))
        structure$atom[0,,drop=FALSE]
        else ram_nonwater_hetero_atoms(full_structure$atom)
      structure$hetero_context_model <- 1L
    } else {
      structure <- read_cif(path, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                            verbose = FALSE)
      structure$hetero_atom <- ram_nonwater_hetero_atoms(structure$atom)
      structure$hetero_context_model <- 1L
    }
  } else {
    if (!is.character(pdb_id) || length(pdb_id) != 1L ||
        is.na(pdb_id) || !grepl("^[[:alnum:]]{4}$", pdb_id)) {
      stop("Enter a valid four-character PDB accession.")
    }
    structure <- tryCatch({
      value <- read_cif(pdb_id, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                        verbose = FALSE)
      value$hetero_atom <- ram_nonwater_hetero_atoms(value$atom)
      value$hetero_context_model <- 1L
      value
    }, error = function(cif_error) {
      tryCatch({
        value <- read_pdb(pdb_id, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                          ATOM.only = TRUE, verbose = FALSE)
        full_value <- tryCatch(
          read_pdb(pdb_id, multi = TRUE, rm.insert = FALSE, rm.alt = FALSE,
                   ATOM.only = FALSE, verbose = FALSE),
          error = function(e) NULL
        )
        value$hetero_atom <- if (is.null(full_value))
          value$atom[0,,drop=FALSE]
          else ram_nonwater_hetero_atoms(full_value$atom)
        value$hetero_context_model <- 1L
        value
      }, error = function(pdb_error) {
        stop(sprintf("Could not retrieve %s: %s",
                     pdb_id, conditionMessage(pdb_error)), call. = FALSE)
      })
    })
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
