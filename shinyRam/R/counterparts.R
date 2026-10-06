# Experimental counterpart discovery for predicted structures.
#
# 3D-Beacons provides a UniProt-centred index of experimental and predicted
# structure models. RamplotR only treats PDBe records explicitly labelled
# EXPERIMENTALLY DETERMINED as experimental counterparts.

ram_uniprot_accession <- function(accession) {
  value <- toupper(trimws(as.character(accession)))
  if(length(value)!=1L || is.na(value) ||
     !grepl("^[A-Z0-9]{6,12}(?:-[0-9]+)?$",value))
    stop("Enter a valid UniProt accession.")
  value
}

ram_3dbeacons_summary_url <- function(accession) {
  accession <- ram_uniprot_accession(accession)
  paste0(
    "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/",
    utils::URLencode(accession,reserved=TRUE),".json"
  )
}

ram_3dbeacons_fetch <- function(url) {
  if(!requireNamespace("jsonlite",quietly=TRUE))
    stop("Install jsonlite to query experimental counterparts.")
  jsonlite::fromJSON(url,simplifyVector=FALSE)
}

ram_counterpart_scalar <- function(value, default=NA_character_) {
  if(is.null(value)) return(default)
  flat <- unlist(value,use.names=FALSE)
  if(!length(flat) || is.na(flat[[1L]])) return(default)
  as.character(flat[[1L]])
}

ram_counterpart_numeric <- function(value) {
  out <- suppressWarnings(as.numeric(ram_counterpart_scalar(value,NA_character_)))
  if(length(out)!=1L || !is.finite(out)) NA_real_ else out
}

ram_counterpart_pdb_id <- function(identifier) {
  value <- toupper(ram_counterpart_scalar(identifier,""))
  hit <- regexpr("^[0-9][A-Z0-9]{3}",value,perl=TRUE)
  if(hit[[1L]]<0L) return(NA_character_)
  regmatches(value,hit)
}

ram_counterpart_chains <- function(entities) {
  if(is.null(entities) || !is.list(entities)) return("")
  chains <- unique(unlist(lapply(entities,function(entity) {
    type <- toupper(ram_counterpart_scalar(entity$entity_type,""))
    if(!identical(type,"POLYMER")) return(character())
    as.character(unlist(entity$chain_ids,use.names=FALSE))
  }),use.names=FALSE))
  chains <- chains[!is.na(chains) & nzchar(chains)]
  paste(chains,collapse=", ")
}

ram_parse_experimental_counterparts <- function(response, accession,
                                                max_results=50L) {
  accession <- ram_uniprot_accession(accession)
  if(!is.list(response)) stop("Invalid 3D-Beacons response.")
  structures <- response$structures
  if(is.null(structures) || !length(structures))
    return(data.frame(
      accession=character(),pdb_id=character(),model_identifier=character(),
      chains=character(),coverage=numeric(),resolution=numeric(),
      uniprot_start=integer(),uniprot_end=integer(),
      experimental_method=character(),model_url=character(),
      stringsAsFactors=FALSE
    ))
  rows <- lapply(structures,function(item) {
    summary <- if(is.list(item)) item$summary else NULL
    if(is.null(summary) || !is.list(summary)) return(NULL)
    provider <- toupper(ram_counterpart_scalar(summary$provider,""))
    category <- toupper(gsub("[-_]"," ",
      ram_counterpart_scalar(summary$model_category,"")))
    if(provider!="PDBE" ||
       !category %in% c("EXPERIMENTALLY DETERMINED","EXPERIMENTAL"))
      return(NULL)
    identifier <- ram_counterpart_scalar(summary$model_identifier,"")
    pdb_id <- ram_counterpart_pdb_id(identifier)
    if(is.na(pdb_id)) return(NULL)
    data.frame(
      accession=accession,
      pdb_id=pdb_id,
      model_identifier=identifier,
      chains=ram_counterpart_chains(summary$entities),
      coverage=ram_counterpart_numeric(summary$coverage),
      resolution=ram_counterpart_numeric(summary$resolution),
      uniprot_start=suppressWarnings(as.integer(
        ram_counterpart_numeric(summary$uniprot_start))),
      uniprot_end=suppressWarnings(as.integer(
        ram_counterpart_numeric(summary$uniprot_end))),
      experimental_method=ram_counterpart_scalar(
        summary$experimental_method,""),
      model_url=ram_counterpart_scalar(summary$model_url,""),
      stringsAsFactors=FALSE
    )
  })
  rows <- rows[!vapply(rows,is.null,logical(1))]
  if(!length(rows)) return(ram_parse_experimental_counterparts(
    list(structures=list()),accession,max_results))
  result <- do.call(rbind,rows)

  # One UniProt mapping may expose multiple records for the same PDB/chain.
  # Keep the best-covered / highest-resolution representative of duplicates.
  coverage_rank <- replace(result$coverage,!is.finite(result$coverage),-Inf)
  resolution_rank <- replace(result$resolution,!is.finite(result$resolution),Inf)
  result <- result[order(-coverage_rank,resolution_rank,result$pdb_id,
                         result$model_identifier),,drop=FALSE]
  key <- paste(result$pdb_id,result$chains,sep="\r")
  result <- result[!duplicated(key),,drop=FALSE]

  max_results <- suppressWarnings(as.integer(max_results))
  if(length(max_results)!=1L || is.na(max_results) || max_results<1L)
    stop("max_results must be a positive integer.")
  utils::head(result,max_results)
}

ram_lookup_experimental_counterparts <- function(
    accession,max_results=50L,fetch=ram_3dbeacons_fetch) {
  accession <- ram_uniprot_accession(accession)
  response <- fetch(ram_3dbeacons_summary_url(accession))
  ram_parse_experimental_counterparts(response,accession,max_results)
}
