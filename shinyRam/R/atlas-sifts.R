# Exact PDBe SIFTS mapping for a selected Atlas experimental entity.
# Input records are parsed from PDBe's updated mmCIF SIFTS and polymer sequence
# scheme categories. Do not interpolate or silently substitute author numbers.

ram_atlas_exact_sifts_map <- function(payload, pdb_id, entity_id, accession) {
  pdb <- toupper(trimws(as.character(pdb_id)))
  accession <- ram_uniprot_accession(accession)
  entity <- suppressWarnings(as.integer(entity_id))
  if(length(pdb)!=1L || is.na(pdb) ||
     !grepl("^[A-Z0-9]{4}$",pdb) ||
     length(entity)!=1L || is.na(entity) || entity<1L)
    stop("Invalid experimental PDB/entity identifier.",call.=FALSE)
  if(!is.list(payload) ||
     !identical(toupper(as.character(payload$pdb_id)),pdb) ||
     !identical(as.character(payload$accession),accession) ||
     !identical(as.character(payload$entity_id),as.character(entity)))
    stop("SIFTS data do not match the selected experimental entity.",
         call.=FALSE)
  raw <- payload$rows
  if(is.null(raw)) raw <- list()
  if(!is.list(raw) || length(raw)>50000L)
    stop("Invalid SIFTS residue mapping records.",call.=FALSE)
  if(!length(raw)) return(ram_canonical_empty_map())

  clean <- function(value) {
    if(is.null(value) || length(value)!=1L || is.na(value))
      return(NA_character_)
    as.character(value)
  }
  records <- lapply(raw,function(row) {
    if(!is.list(row)) stop("Malformed exact SIFTS record.",call.=FALSE)
    chain <- clean(row$chain)
    ins <- clean(row$insertion_code)
    target <- clean(row$uniprot_accession)
    asym <- clean(row$struct_asym_id)
    resi <- suppressWarnings(as.integer(row$resi))
    unp <- suppressWarnings(as.integer(row$uniprot_resi))
    assigned_entity <- suppressWarnings(as.integer(row$entity_id))
    if(is.na(chain) || !nzchar(chain) ||
       is.na(ins) || ins %in% c(".","?") ||
       is.na(target) || !identical(target,accession) ||
       is.na(asym) || !nzchar(asym) ||
       length(resi)!=1L || is.na(resi) ||
       length(unp)!=1L || is.na(unp) || unp<1L ||
       length(assigned_entity)!=1L ||
       is.na(assigned_entity) || assigned_entity!=entity)
      stop("Invalid or mismatched exact SIFTS residue.",call.=FALSE)
    data.frame(
      chain=chain,resi=resi,insertion_code=ins,
      uniprot_accession=target,uniprot_resi=unp,
      entity_id=entity,struct_asym_id=asym,
      identity=NA_real_,coverage=NA_real_,
      canonical_source="PDBe updated mmCIF exact SIFTS",
      observed=isTRUE(row$observed),
      stringsAsFactors=FALSE
    )
  })
  out <- unique(do.call(rbind,records))
  # Conflicts between distinct asym chains or different UniProt targets are
  # not resolved by taking the first row; canonical_join marks ambiguous.
  out
}

ram_atlas_sifts_summary <- function(mapping) {
  if(!is.data.frame(mapping) || !all(c("chain","resi","insertion_code",
                                       "uniprot_resi","observed") %in%
                                     names(mapping)))
    stop("Invalid exact SIFTS mapping.",call.=FALSE)
  if(!nrow(mapping))
    return(list(rows=0L,unique_residues=0L,observed_residues=0L,
                conflicting_residues=0L,chains=character()))
  key <- paste(mapping$chain,mapping$resi,mapping$insertion_code,sep="\r")
  target <- paste(mapping$uniprot_accession,mapping$uniprot_resi,sep=":")
  conflicts <- vapply(split(target,key),function(x)
    length(unique(x))>1L,logical(1L))
  list(rows=nrow(mapping),unique_residues=length(unique(key)),
       observed_residues=length(unique(key[!is.na(mapping$observed) &
                                            mapping$observed])),
       conflicting_residues=sum(conflicts),
       chains=unique(as.character(mapping$chain)))
}
