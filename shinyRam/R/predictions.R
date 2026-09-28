# Optional prediction confidence; never interpret experimental B factors as pLDDT.
# Supports AlphaFold DB/AF2, AlphaFold 3 and ESMFold independently of Shiny.

ram_prediction_key <- function(chain, resi, insertion_code = "") {
  paste(as.character(chain), as.character(resi),
        ifelse(is.na(insertion_code), "", as.character(insertion_code)), sep = "\r")
}

ram_plddt_category <- function(value) {
  result <- rep("Unavailable", length(value))
  result[is.finite(value) & value < 50] <- "Very low"
  result[is.finite(value) & value >= 50 & value < 70] <- "Low"
  result[is.finite(value) & value >= 70 & value < 90] <- "Confident"
  result[is.finite(value) & value >= 90 & value <= 100] <- "Very high"
  result
}

ram_prediction_from_atoms <- function(pdb, torsions, source) {
  permitted <- c("alphafold_db", "alphafold2", "alphafold3", "esmfold",
                 "other_prediction")
  if (length(source) != 1L || !source %in% permitted)
    stop("Declare the prediction source before interpreting B factors as pLDDT.")
  atoms <- pdb$atom
  field <- intersect(c("b", "B", "b_iso_or_equiv"), names(atoms))
  scores <- rep(NA_real_, nrow(torsions))
  if (length(field)) {
    b <- suppressWarnings(as.numeric(atoms[[field[[1L]]]]))
    if (length(b) == nrow(atoms) && all(is.na(b) | is.finite(b) &
                                       b >= 0 & b <= 100)) {
      if (!"insert" %in% names(atoms)) atoms$insert <- ""
      atoms$insert[is.na(atoms$insert)] <- ""
      ids <- ram_prediction_key(atoms$chain, atoms$resno, atoms$insert)
      row_ids <- ram_prediction_key(torsions$chain, torsions$resi,
                                    torsions$insertion_code)
      # Carbon-alpha scores represent protein residues for AF2 and ESMFold.
      # For AF3 the full JSON can supply separate per-atom confidence.
      candidates <- which(atoms$elety == "CA" & is.finite(b) &
                           (if ("alt" %in% names(atoms))
                              is.na(atoms$alt) | atoms$alt %in% c("", "A")
                            else TRUE))
      best <- candidates[!duplicated(ids[candidates])]
      scores <- b[best][match(row_ids, ids[best])]
    }
  }
  data.frame(chain = as.character(torsions$chain),
             resi = as.integer(torsions$resi),
             insertion_code = as.character(torsions$insertion_code),
             plddt = as.numeric(scores),
             confidence_category = ram_plddt_category(scores),
             stringsAsFactors = FALSE)
}

ram_read_confidence_json <- function(path, max_bytes = 32000000) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Install the jsonlite package to read prediction confidence files.")
  if (length(path) != 1L || !file.exists(path) ||
      !is.finite(file.info(path)$size) || file.info(path)$size > max_bytes)
    stop("Confidence JSON is missing or larger than the 32 MB upload limit.")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

ram_square_pae <- function(value, max_tokens = 1600L) {
  if (!is.list(value) || !length(value) || length(value) > max_tokens)
    stop("PAE must be a nonempty square matrix with at most 1,600 tokens.")
  n <- length(value)
  if (!all(vapply(value, function(row) length(row) == n &&
                  is.atomic(unlist(row)), logical(1))))
    stop("The PAE matrix must have equally sized numeric rows.")
  flat <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if (length(flat) != n * n || any(!is.finite(flat) | flat < 0 | flat > 100))
    stop("PAE contains nonnumeric, missing or out-of-range values.")
  matrix(flat, nrow = n, byrow = TRUE)
}

ram_prediction_json <- function(json, torsions, atoms = NULL,
                                source = "alphafold2") {
  if (is.list(json) && is.null(names(json)) && length(json) == 1L)
    json <- json[[1L]] # AlphaFold DB PAE is an array containing one object.
  if (!is.list(json) || is.null(names(json)))
    stop("Unrecognised prediction-confidence JSON schema.")
  n <- nrow(torsions)
  result <- list(pae = NULL, pae_rows = integer(),
                 plddt = rep(NA_real_, n), notes = character(),
                 ptm = NA_real_, iptm = NA_real_)
  for (field in c("ptm", "iptm")) {
    value <- suppressWarnings(as.numeric(unlist(json[[field]])))
    if (length(value) == 1L && is.finite(value) && value >= 0 && value <= 1)
      result[[field]] <- value
  }
  # AFDB confidenceScore contains per-residue pLDDT. Avoid matching arrays
  # to multimeric structures unless all residue rows align unambiguously.
  score <- json$confidenceScore
  if (!is.null(score)) {
    score <- suppressWarnings(as.numeric(unlist(score)))
    if (length(score) != n || length(unique(torsions$chain)) != 1L)
      stop("Confidence scores cannot safely be aligned with this structure.")
    if (any(!is.finite(score) | score < 0 | score > 100))
      stop("Invalid per-residue pLDDT score.")
    result$plddt <- score
  }
  # AF3 atom_plddts are per ATOM, not per token or residue. Only match if
  # Bio3D preserved full atom order and AF3 atom_chain_ids agrees.
  if (!is.null(json$atom_plddts)) {
    values <- suppressWarnings(as.numeric(unlist(json$atom_plddts)))
    if (is.null(atoms) || length(values) != nrow(atoms)) {
      result$notes <- c(result$notes,
        "AF3 atom confidence was not applied: atom order/count could not be verified.")
    } else {
      chain_ids <- unlist(json$atom_chain_ids, use.names = FALSE)
      if (length(chain_ids) && !identical(as.character(chain_ids),
                                          as.character(atoms$chain))) {
        result$notes <- c(result$notes,
          "AF3 atom confidence was not applied: chain order differs.")
      } else if (any(!is.finite(values) | values < 0 | values > 100)) {
        result$notes <- c(result$notes, "Invalid AF3 atom confidence values.")
      } else {
        if (!"insert" %in% names(atoms)) atoms$insert <- ""
        atoms$insert[is.na(atoms$insert)] <- ""
        keys <- ram_prediction_key(atoms$chain, atoms$resno, atoms$insert)
        residues <- ram_prediction_key(torsions$chain, torsions$resi,
                                       torsions$insertion_code)
        ca <- which(atoms$elety == "CA")
        ca <- ca[!duplicated(keys[ca])]
        result$plddt <- values[ca][match(residues, keys[ca])]
      }
    }
  }
  grid <- if (!is.null(json$predicted_aligned_error))
    json$predicted_aligned_error else json$pae
  if (is.null(grid)) return(result)
  pae <- ram_square_pae(grid)
  if (!is.null(json$token_chain_ids) || !is.null(json$token_res_ids)) {
    # AF3 tokens can include ligand atoms and nucleic acids. Only protein
    # tokens with a unique, verifiable chain/residue match may be visualised.
    chains <- as.character(unlist(json$token_chain_ids, use.names = FALSE))
    residue_ids <- suppressWarnings(as.integer(unlist(json$token_res_ids,
                                                        use.names = FALSE)))
    if (length(chains) != nrow(pae) || length(residue_ids) != nrow(pae)) {
      result$notes <- c(result$notes, "PAE omitted: AF3 token IDs are incomplete.")
      return(result)
    }
    keys <- ram_prediction_key(chains, residue_ids)
    torsion_keys <- ram_prediction_key(torsions$chain, torsions$resi)
    if (anyDuplicated(torsion_keys) ||
        any(!is.na(torsions$insertion_code) &
            nzchar(torsions$insertion_code))) {
      result$notes <- c(result$notes,
                        "PAE omitted: insertion codes need explicit token mapping.")
      return(result)
    }
    matches <- match(keys, torsion_keys)
    keep <- which(!is.na(matches) & !duplicated(matches))
    if (!length(keep)) {
      result$notes <- c(result$notes,
                        "PAE omitted: AF3 protein tokens did not match the structure.")
      return(result)
    }
    result$pae <- pae[keep, keep, drop = FALSE]
    result$pae_rows <- matches[keep]
    if (length(keep) != nrow(pae))
      result$notes <- c(result$notes,
                        "Only matched protein tokens are displayed; other tokens omitted.")
  } else if (length(unique(torsions$chain)) == 1L &&
             nrow(pae) == nrow(torsions)) {
    # AFDB and monomeric AF2 arrays follow the sequence row order.
    result$pae <- pae
    result$pae_rows <- seq_len(nrow(torsions))
  } else {
    result$notes <- c(result$notes,
      "PAE omitted: no verifiable token-to-residue mapping for this structure.")
  }
  result
}

ram_apply_prediction <- function(classified, confidence) {
  if (is.null(confidence)) return(classified)
  keys <- ram_prediction_key(classified$chain, classified$resi,
                             classified$insertion_code)
  confidence_keys <- ram_prediction_key(confidence$residues$chain,
    confidence$residues$resi, confidence$residues$insertion_code)
  if (anyDuplicated(confidence_keys))
    stop("Prediction residue identifiers must be unique.")
  classified$plddt <- confidence$residues$plddt[match(keys, confidence_keys)]
  classified$confidence_category <- ram_plddt_category(classified$plddt)
  classified
}

ram_confidence_review <- function(data) {
  if (!"plddt" %in% names(data)) return(rep("Unavailable", nrow(data)))
  category <- ram_plddt_category(data$plddt)
  geometry_issue <- !is.na(data$region) & data$region == "Not allowed"
  result <- rep("Other", nrow(data))
  result[is.finite(data$plddt) & data$plddt < 70] <-
    "Lower-confidence prediction"
  result[is.finite(data$plddt) & data$plddt >= 90 & geometry_issue] <-
    "High confidence · Ramachandran outlier"
  result[is.finite(data$plddt) & data$plddt >= 90 & !geometry_issue] <-
    "High confidence · geometry in range"
  result[category == "Unavailable"] <- "Confidence unavailable"
  result
}

# Optional native AFDB retrieval. Never accept URLs supplied by clients; only
# use the officially returned files URL on the approved AlphaFold host.
ram_afdb_entry <- function(accession,
                            fetch = function(url) jsonlite::fromJSON(url,
                                                     simplifyVector = FALSE)) {
  if (!grepl("^[A-Za-z0-9]{6,12}(?:-[0-9]+)?$", accession))
    stop("Enter a valid UniProt accession for AlphaFold DB.")
  endpoint <- paste0("https://alphafold.ebi.ac.uk/api/prediction/",
                     utils::URLencode(toupper(accession), reserved = TRUE))
  response <- fetch(endpoint)
  if (!is.list(response) || !length(response)) stop("No AlphaFold DB model found.")
  entry <- if (is.null(names(response))) response[[1L]] else response
  if (!is.list(entry)) stop("Unexpected AlphaFold DB API response.")
  valid_url <- function(url) {
    is.character(url) && length(url) == 1L && !is.na(url) &&
      grepl("^https://alphafold[.]ebi[.]ac[.]uk/files/[A-Za-z0-9._/%-]+$", url)
  }
  pdb <- entry$pdbUrl
  cif <- entry$cifUrl
  url <- if (valid_url(pdb)) pdb else if (valid_url(cif)) cif else NULL
  if (is.null(url)) stop("AlphaFold DB did not return a supported structure URL.")
  pae <- entry$paeDocUrl
  list(structure_url = url,
       structure_format = if (identical(url, pdb)) "pdb" else "cif",
       pae_url = if (valid_url(pae)) pae else NULL,
       accession = toupper(accession),
       model_id = if (!is.null(entry$modelEntityId)) entry$modelEntityId else NA_character_)
}

ram_download_afdb <- function(entry,
                              downloader = function(url, file)
                                utils::download.file(url, file, quiet = TRUE,
                                                     mode = "wb")) {
  ext <- if (identical(entry$structure_format, "pdb")) ".pdb" else ".cif"
  coord <- tempfile("ram-afdb-", fileext = ext)
  pae <- if (!is.null(entry$pae_url)) tempfile("ram-afdb-", fileext = ".json") else NULL
  on_error <- TRUE
  on.exit(if (on_error) unlink(c(coord, pae)), add = TRUE)
  downloader(entry$structure_url, coord)
  if (!file.exists(coord) || file.info(coord)$size == 0L)
    stop("AlphaFold DB structure download failed.")
  notes <- character()
  if (!is.null(pae)) {
    tryCatch({
      downloader(entry$pae_url, pae)
      if (!file.exists(pae) || !is.finite(file.info(pae)$size) ||
          file.info(pae)$size > 32000000)
        stop("Missing/oversized PAE file.")
    }, error = function(e) {
      notes <<- paste("PAE unavailable:", conditionMessage(e))
      unlink(pae)
      pae <<- NULL
    })
  }
  on_error <- FALSE
  list(structure = coord, pae = pae, notes = notes,
       original_name = paste0("AF-", entry$accession, ext))
}


# Combine B-factor confidence and independently verified optional JSON without
# altering any of the underlying torsion or Ramachandran calculations.
ram_prepare_prediction <- function(pdb, torsions, source, sidecar = NULL,
                                   summary_file = NULL, notes = character(),
                                   model_id = "") {
  baseline <- ram_prediction_from_atoms(pdb, torsions, source)
  mapped <- if (!is.null(sidecar) && nzchar(sidecar))
    ram_prediction_json(ram_read_confidence_json(sidecar), torsions,
                        atoms = pdb$atom, source = source) else NULL
  if (!is.null(mapped) && any(is.finite(mapped$plddt))) {
    present <- is.finite(mapped$plddt)
    baseline$plddt[present] <- mapped$plddt[present]
    baseline$confidence_category <- ram_plddt_category(baseline$plddt)
  }
  summary <- if (!is.null(summary_file) && nzchar(summary_file))
    ram_prediction_json(ram_read_confidence_json(summary_file), torsions,
                        source = source) else NULL
  metric <- function(key) {
    if (!is.null(summary) && is.finite(summary[[key]])) return(summary[[key]])
    if (!is.null(mapped) && is.finite(mapped[[key]])) return(mapped[[key]])
    NA_real_
  }
  list(source = source, residues = baseline,
       pae = if (!is.null(mapped)) mapped$pae else NULL,
       pae_rows = if (!is.null(mapped)) mapped$pae_rows else integer(),
       ptm = metric("ptm"), iptm = metric("iptm"),
       notes = c(notes, if (!is.null(mapped)) mapped$notes),
       model_id = model_id,
       confidence_file = if (!is.null(sidecar)) basename(sidecar) else "")
}

# Cap Plotly payload size. Sampling is explicitly labelled: never report
# sampled PAE values as a complete matrix or hide a chain/domain boundary.
ram_pae_plot_data <- function(prediction, torsions, max_display = 400L) {
  if (is.null(prediction) || is.null(prediction$pae)) return(NULL)
  pae <- prediction$pae
  if (nrow(pae) != ncol(pae) || nrow(pae) != length(prediction$pae_rows))
    stop("PAE values and residue indices are inconsistent.")
  n <- nrow(pae)
  chosen <- if (n <= max_display) seq_len(n) else
    unique(as.integer(round(seq(1, n, length.out = max_display))))
  rows <- prediction$pae_rows[chosen]
  if (anyNA(rows) || any(rows < 1L | rows > nrow(torsions)))
    stop("PAE residue mapping is invalid.")
  labels <- paste0(torsions$chain[rows], ":",
    torsions$resi[rows], torsions$insertion_code[rows])
  list(z = lapply(chosen, function(i) unname(as.numeric(pae[i, chosen]))),
       labels = labels,
       residues = lapply(rows, function(i) list(
         chain = as.character(torsions$chain[[i]]),
         resi = as.integer(torsions$resi[[i]]),
         insertion_code = as.character(torsions$insertion_code[[i]]))),
       total_tokens = n, displayed_tokens = length(chosen),
       downsampled = length(chosen) != n)
}
