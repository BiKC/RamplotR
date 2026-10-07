# Canonical residue coordinates for cross-structure analysis.
#
# RamplotR preserves PDB author numbering/insertion codes as the local
# coordinate system. Canonical UniProt positions are an additional annotation
# used to compare many structures of the same protein.
#
# SIFTS range mappings are deliberately expanded only when the author-residue
# and UniProt ranges are provably one-to-one. Nonlinear/ambiguous mappings stay
# unresolved until exact residue-level SIFTS data are available.

ram_canonical_empty_segments <- function() {
  data.frame(
    pdb_id=character(),uniprot_accession=character(),
    uniprot_identifier=character(),entity_id=integer(),
    chain=character(),struct_asym_id=character(),
    unp_start=integer(),unp_end=integer(),
    pdb_start=integer(),pdb_end=integer(),
    author_start=integer(),author_end=integer(),
    author_start_insertion=character(),author_end_insertion=character(),
    identity=numeric(),coverage=numeric(),
    safe_linear=logical(),stringsAsFactors=FALSE
  )
}

ram_canonical_empty_map <- function() {
  data.frame(
    chain=character(),resi=integer(),insertion_code=character(),
    uniprot_accession=character(),uniprot_resi=integer(),
    entity_id=integer(),struct_asym_id=character(),
    identity=numeric(),coverage=numeric(),
    canonical_source=character(),stringsAsFactors=FALSE
  )
}

ram_uniprot_accession <- function(accession) {
  value <- toupper(trimws(as.character(accession)))
  if (length(value)!=1L || is.na(value) || !nzchar(value) ||
      !grepl("^[A-Z0-9]{6,12}(-[1-9][0-9]*)?$",value))
    stop("Expected a UniProt accession or isoform accession.",call.=FALSE)
  value
}

ram_canonical_string <- function(value, default="") {
  if (is.null(value) || !length(value) || is.na(value[[1L]])) return(default)
  as.character(value[[1L]])
}

ram_canonical_integer <- function(value) {
  if (is.null(value) || !length(value)) return(NA_integer_)
  suppressWarnings(as.integer(value[[1L]]))
}

ram_canonical_number <- function(value) {
  if (is.null(value) || !length(value)) return(NA_real_)
  suppressWarnings(as.numeric(value[[1L]]))
}

ram_sifts_payload_segments <- function(payload,pdb_id) {
  pdb_id <- tolower(trimws(as.character(pdb_id)))
  if (length(pdb_id)!=1L || is.na(pdb_id) ||
      !grepl("^[a-z0-9]{4}$",pdb_id))
    stop("Expected a four-character PDB accession.",call.=FALSE)
  if (!is.list(payload) || !length(payload)) return(ram_canonical_empty_segments())

  names_lower <- tolower(names(payload))
  entry_index <- match(pdb_id,names_lower)
  if (is.na(entry_index) && length(payload)==1L) entry_index <- 1L
  if (is.na(entry_index)) return(ram_canonical_empty_segments())
  entry <- payload[[entry_index]]
  uniprot <- if (is.list(entry)) entry[["UniProt"]] else NULL
  if (!is.list(uniprot) || !length(uniprot))
    return(ram_canonical_empty_segments())

  rows <- list()
  add_mapping <- function(accession,identifier,mapping) {
    start <- mapping[["start"]]
    end <- mapping[["end"]]
    rows[[length(rows)+1L]] <<- data.frame(
      pdb_id=toupper(pdb_id),
      uniprot_accession=ram_uniprot_accession(accession),
      uniprot_identifier=ram_canonical_string(identifier),
      entity_id=ram_canonical_integer(mapping[["entity_id"]]),
      chain=ram_canonical_string(mapping[["chain_id"]]),
      struct_asym_id=ram_canonical_string(mapping[["struct_asym_id"]]),
      unp_start=ram_canonical_integer(mapping[["unp_start"]]),
      unp_end=ram_canonical_integer(mapping[["unp_end"]]),
      pdb_start=ram_canonical_integer(if(is.list(start)) start[["residue_number"]] else NULL),
      pdb_end=ram_canonical_integer(if(is.list(end)) end[["residue_number"]] else NULL),
      author_start=ram_canonical_integer(if(is.list(start)) start[["author_residue_number"]] else NULL),
      author_end=ram_canonical_integer(if(is.list(end)) end[["author_residue_number"]] else NULL),
      author_start_insertion=ram_canonical_string(
        if(is.list(start)) start[["author_insertion_code"]] else NULL),
      author_end_insertion=ram_canonical_string(
        if(is.list(end)) end[["author_insertion_code"]] else NULL),
      identity=ram_canonical_number(mapping[["identity"]]),
      coverage=ram_canonical_number(mapping[["coverage"]]),
      safe_linear=FALSE,
      stringsAsFactors=FALSE
    )
  }

  for (accession in names(uniprot)) {
    info <- uniprot[[accession]]
    if (!is.list(info)) next
    identifier <- info[["identifier"]]
    mappings <- info[["mappings"]]
    if (is.data.frame(mappings)) {
      # jsonlite can simplify a flat response to a data frame, but nested
      # start/end fields are normally kept as lists. Convert row-wise.
      for (i in seq_len(nrow(mappings))) {
        mapping <- lapply(mappings,function(column) column[[i]])
        add_mapping(accession,identifier,mapping)
      }
    } else if (is.list(mappings)) {
      for (mapping in mappings)
        if (is.list(mapping)) add_mapping(accession,identifier,mapping)
    }
  }

  if (!length(rows)) return(ram_canonical_empty_segments())
  out <- do.call(rbind,rows)
  ram_sifts_mark_safe_segments(out)
}

ram_sifts_normalize_segments <- function(segments,pdb_id="") {
  if (is.null(segments) || !length(segments))
    return(ram_canonical_empty_segments())
  if (is.data.frame(segments)) {
    rows <- lapply(seq_len(nrow(segments)),function(i)
      as.list(segments[i,,drop=FALSE]))
  } else if (is.list(segments) &&
             !is.null(names(segments)) &&
             "uniprot_accession" %in% names(segments)) {
    rows <- list(segments)
  } else if (is.list(segments) &&
             all(vapply(segments,is.list,logical(1L)))) {
    rows <- segments
  } else stop("SIFTS segments must be a data frame or list of records.",
              call.=FALSE)

  get <- function(row,name,default=NULL) {
    value <- row[[name]]
    if (is.null(value) || !length(value)) return(default)
    value[[1L]]
  }
  out <- lapply(rows,function(row) data.frame(
    pdb_id=toupper(ram_canonical_string(get(row,"pdb_id",pdb_id))),
    uniprot_accession=ram_uniprot_accession(
      get(row,"uniprot_accession")),
    uniprot_identifier=ram_canonical_string(get(row,"uniprot_identifier")),
    entity_id=suppressWarnings(as.integer(get(row,"entity_id",NA))),
    chain=ram_canonical_string(get(row,"chain")),
    struct_asym_id=ram_canonical_string(get(row,"struct_asym_id")),
    unp_start=suppressWarnings(as.integer(get(row,"unp_start",NA))),
    unp_end=suppressWarnings(as.integer(get(row,"unp_end",NA))),
    pdb_start=suppressWarnings(as.integer(get(row,"pdb_start",NA))),
    pdb_end=suppressWarnings(as.integer(get(row,"pdb_end",NA))),
    author_start=suppressWarnings(as.integer(get(row,"author_start",NA))),
    author_end=suppressWarnings(as.integer(get(row,"author_end",NA))),
    author_start_insertion=ram_canonical_string(
      get(row,"author_start_insertion")),
    author_end_insertion=ram_canonical_string(
      get(row,"author_end_insertion")),
    identity=suppressWarnings(as.numeric(get(row,"identity",NA))),
    coverage=suppressWarnings(as.numeric(get(row,"coverage",NA))),
    safe_linear=FALSE,stringsAsFactors=FALSE
  ))
  ram_sifts_mark_safe_segments(do.call(rbind,out))
}

ram_sifts_mark_safe_segments <- function(segments) {
  if (!is.data.frame(segments) || !nrow(segments)) {
    out <- ram_canonical_empty_segments()
    if (is.data.frame(segments) && nrow(segments)==0L) return(out)
    return(out)
  }
  required <- c("author_start","author_end","unp_start","unp_end",
                "author_start_insertion","author_end_insertion")
  if (!all(required %in% names(segments)))
    stop("SIFTS segment table is missing required range fields.",call.=FALSE)
  pdb_span <- segments$author_end-segments$author_start
  unp_span <- segments$unp_end-segments$unp_start
  finite <- is.finite(segments$author_start) & is.finite(segments$author_end) &
    is.finite(segments$unp_start) & is.finite(segments$unp_end)
  no_insert <- ifelse(is.na(segments$author_start_insertion),"",
                      segments$author_start_insertion)=="" &
    ifelse(is.na(segments$author_end_insertion),"",
           segments$author_end_insertion)==""
  segments$safe_linear <- finite & no_insert & pdb_span>=0L & unp_span>=0L &
    pdb_span==unp_span
  segments
}

ram_sifts_expand_safe <- function(segments) {
  if (!is.data.frame(segments) || !nrow(segments))
    return(ram_canonical_empty_map())
  if (!"safe_linear" %in% names(segments))
    segments <- ram_sifts_mark_safe_segments(segments)
  safe <- which(!is.na(segments$safe_linear) & segments$safe_linear)
  if (!length(safe)) return(ram_canonical_empty_map())

  rows <- lapply(safe,function(i) {
    n <- segments$author_end[[i]]-segments$author_start[[i]]+1L
    data.frame(
      chain=rep(as.character(segments$chain[[i]]),n),
      resi=seq.int(segments$author_start[[i]],segments$author_end[[i]]),
      insertion_code=rep("",n),
      uniprot_accession=rep(as.character(segments$uniprot_accession[[i]]),n),
      uniprot_resi=seq.int(segments$unp_start[[i]],segments$unp_end[[i]]),
      entity_id=rep(as.integer(segments$entity_id[[i]]),n),
      struct_asym_id=rep(as.character(segments$struct_asym_id[[i]]),n),
      identity=rep(as.numeric(segments$identity[[i]]),n),
      coverage=rep(as.numeric(segments$coverage[[i]]),n),
      canonical_source=rep("PDBe SIFTS range mapping",n),
      stringsAsFactors=FALSE
    )
  })
  do.call(rbind,rows)
}

ram_afdb_canonical_map <- function(data,accession) {
  accession <- ram_uniprot_accession(accession)
  required <- c("chain","resi","insertion_code")
  if (!is.data.frame(data) || !all(required %in% names(data)))
    stop("AlphaFold DB mapping requires residue identifiers.",call.=FALSE)
  if (!nrow(data)) return(ram_canonical_empty_map())
  chains <- unique(as.character(data$chain))
  if (length(chains)!=1L)
    stop("AlphaFold DB canonical mapping expects one protein chain.",call.=FALSE)
  insertion <- ifelse(is.na(data$insertion_code),"",
                      as.character(data$insertion_code))
  if (any(nzchar(insertion)))
    stop("Unexpected insertion codes in AlphaFold DB coordinates.",call.=FALSE)
  valid <- is.finite(data$resi)
  out <- data.frame(
    chain=as.character(data$chain[valid]),
    resi=as.integer(data$resi[valid]),
    insertion_code="",
    uniprot_accession=accession,
    uniprot_resi=as.integer(data$resi[valid]),
    entity_id=NA_integer_,struct_asym_id="",
    identity=1,coverage=NA_real_,
    canonical_source="AlphaFold DB UniProt numbering",
    stringsAsFactors=FALSE
  )
  unique(out)
}

ram_canonical_join <- function(data,mapping) {
  if (!is.data.frame(data)) stop("Expected a residue data frame.",call.=FALSE)
  required <- c("chain","resi","insertion_code")
  if (!all(required %in% names(data)))
    stop("Residue data are missing chain/residue identifiers.",call.=FALSE)

  out <- data
  out$uniprot_accession <- NA_character_
  out$uniprot_resi <- NA_integer_
  out$canonical_source <- NA_character_
  out$canonical_identity <- NA_real_
  out$canonical_coverage <- NA_real_
  out$canonical_status <- "unmapped"
  if (!is.data.frame(mapping) || !nrow(mapping)) return(out)

  map_required <- c("chain","resi","insertion_code","uniprot_accession",
                    "uniprot_resi","canonical_source","identity","coverage")
  if (!all(map_required %in% names(mapping)))
    stop("Canonical mapping table is missing required fields.",call.=FALSE)

  norm_ins <- function(x) ifelse(is.na(x),"",as.character(x))
  key <- paste(as.character(out$chain),as.integer(out$resi),
               norm_ins(out$insertion_code),sep="\r")
  map_key <- paste(as.character(mapping$chain),as.integer(mapping$resi),
                   norm_ins(mapping$insertion_code),sep="\r")

  groups <- split(seq_len(nrow(mapping)),map_key)
  canonical <- vector("list",length(groups))
  ambiguous <- character()
  names(canonical) <- names(groups)
  for (name in names(groups)) {
    ix <- groups[[name]]
    target <- paste(mapping$uniprot_accession[ix],
                    mapping$uniprot_resi[ix],sep=":")
    target <- unique(target[!is.na(target)])
    if (length(target)!=1L) {
      ambiguous <- c(ambiguous,name)
      next
    }
    canonical[[name]] <- mapping[ix[[1L]],,drop=FALSE]
  }
  canonical <- canonical[!vapply(canonical,is.null,logical(1L))]
  if (length(canonical)) {
    map <- do.call(rbind,canonical)
    mk <- paste(as.character(map$chain),as.integer(map$resi),
                norm_ins(map$insertion_code),sep="\r")
    idx <- match(key,mk)
    good <- !is.na(idx)
    out$uniprot_accession[good] <- as.character(map$uniprot_accession[idx[good]])
    out$uniprot_resi[good] <- as.integer(map$uniprot_resi[idx[good]])
    out$canonical_source[good] <- as.character(map$canonical_source[idx[good]])
    out$canonical_identity[good] <- as.numeric(map$identity[idx[good]])
    out$canonical_coverage[good] <- as.numeric(map$coverage[idx[good]])
    out$canonical_status[good] <- "mapped"
  }
  if (length(ambiguous))
    out$canonical_status[key %in% ambiguous] <- "ambiguous"
  out
}
