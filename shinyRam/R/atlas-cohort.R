# Collect explicitly verified Atlas SIFTS mappings in UniProt coordinates.
# This is a structural-coverage inventory, not a conformational-state cluster
# or a count of independent experimental observations.

ram_atlas_cohort_table <- function(verified,accession) {
  accession <- ram_uniprot_accession(accession)
  empty <- data.frame(
    pdb_id=character(),entity_id=integer(),chain=character(),
    resi=integer(),insertion_code=character(),
    uniprot_accession=character(),uniprot_resi=integer(),
    observed=logical(),mapping_status=character(),
    struct_asym_id=character(),mon_id=character(),
    canonical_source=character(),
    stringsAsFactors=FALSE
  )
  if(!is.list(verified) || !length(verified)) return(empty)
  rows <- list()
  for(key in names(verified)) {
    record <- verified[[key]]
    if(!is.list(record) || !identical(record$state,"mapped")) next
    if(!grepl("^[A-Z0-9]{4}_[1-9][0-9]*$",key))
      stop("Atlas verification has an invalid entity identifier.",call.=FALSE)
    map <- record$mapping
    if(!is.data.frame(map) || !all(c("chain","resi","insertion_code",
        "uniprot_accession","uniprot_resi","entity_id","struct_asym_id",
        "observed","canonical_source") %in% names(map)))
      stop("Atlas verified mapping is missing required columns.",call.=FALSE)
    if(!nrow(map)) next
    parts <- strsplit(key,"_",fixed=TRUE)[[1L]]
    if(any(is.na(map$uniprot_accession)) ||
       any(map$uniprot_accession!=accession) ||
       any(map$entity_id!=as.integer(parts[[2L]]),na.rm=FALSE))
      stop("Verified entity does not match the current UniProt cohort.",
           call.=FALSE)
    local <- paste(as.character(map$chain),as.integer(map$resi),
                   as.character(map$insertion_code),sep="\r")
    target <- paste(map$uniprot_accession,map$uniprot_resi,sep=":")
    ambiguous <- vapply(split(target,local),function(values)
      length(unique(values))>1L,logical(1L))
    state <- ifelse(ambiguous[local],"ambiguous","mapped")
    rows[[length(rows)+1L]] <- data.frame(
      pdb_id=parts[[1L]],entity_id=as.integer(parts[[2L]]),
      chain=as.character(map$chain),resi=as.integer(map$resi),
      insertion_code=as.character(map$insertion_code),
      uniprot_accession=as.character(map$uniprot_accession),
      uniprot_resi=as.integer(map$uniprot_resi),
      observed=as.logical(map$observed),
      mapping_status=state,
      struct_asym_id=as.character(map$struct_asym_id),
      mon_id=if("mon_id" %in% names(map)) as.character(map$mon_id)
        else rep(NA_character_,nrow(map)),
      canonical_source=as.character(map$canonical_source),
      stringsAsFactors=FALSE
    )
  }
  if(!length(rows)) return(empty)
  unique(do.call(rbind,rows))
}

ram_atlas_cohort_summary <- function(verified,accession) {
  table <- ram_atlas_cohort_table(verified,accession)
  verified_keys <- names(verified)[vapply(verified,function(x)
    is.list(x) && identical(x$state,"mapped"),logical(1L))]
  good <- !is.na(table$observed) & table$observed &
          table$mapping_status=="mapped"
  points <- unique(table[good,c("uniprot_accession","uniprot_resi"),
                         drop=FALSE])
  local <- unique(table[table$mapping_status=="ambiguous",
    c("pdb_id","entity_id","chain","resi","insertion_code"),drop=FALSE])
  list(verified_entities=length(verified_keys),
       exact_rows=nrow(table),observed_positions=nrow(points),
       ambiguous_local_residues=nrow(local),
       observed_entity_positions=length(unique(paste(
         table$pdb_id[good],table$entity_id[good],
         table$uniprot_resi[good],sep=":"))))
}

ram_atlas_cohort_position_support <- function(verified,accession) {
  data <- ram_atlas_cohort_table(verified,accession)
  ok <- !is.na(data$observed) & data$observed &
        data$mapping_status=="mapped"
  data <- data[ok,,drop=FALSE]
  if(!nrow(data)) return(data.frame(
    uniprot_resi=integer(),verified_entity_count=integer()))
  entity_key <- paste(data$pdb_id,data$entity_id,sep="_")
  pairs <- unique(data.frame(
    uniprot_resi=data$uniprot_resi,entity_key=entity_key))
  tab <- table(pairs$uniprot_resi)
  data.frame(uniprot_resi=as.integer(names(tab)),
    verified_entity_count=as.integer(tab))[order(as.integer(names(tab))),,
      drop=FALSE]
}
