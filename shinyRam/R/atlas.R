# Cohort-oriented experimental Atlas page merging.
#
# A page advances by the number of RCSB search hits returned, not by the
# number of successfully enriched entities. Missing metadata must not cause
# subsequent pages to repeat earlier offsets. The first matching entity
# record wins if the archive repeats an identifier across pages.

ram_atlas_page_integer <- function(value,label) {
  n <- suppressWarnings(as.numeric(value))
  if(length(n)!=1L || !is.finite(n) || n<0 || n != floor(n) ||
     n>1000000)
    stop(paste0("Invalid Atlas ",label,"."),call.=FALSE)
  as.integer(n)
}

ram_atlas_record_key <- function(record) {
  if(!is.list(record)) stop("Atlas entity record is not an object.",call.=FALSE)
  pdb <- toupper(as.character(record$pdb_id))
  entity <- as.character(record$entity_id)
  if(length(pdb)!=1L || is.na(pdb) || !grepl("^[A-Z0-9]{4}$",pdb) ||
     length(entity)!=1L || is.na(entity) ||
     !grepl("^[1-9][0-9]*$",entity))
    stop("Atlas record has an invalid PDB/entity identifier.",call.=FALSE)
  paste0(pdb,"_",entity)
}

ram_atlas_merge_page <- function(previous=NULL,page) {
  if(!is.list(page)) stop("Atlas page must be a list.",call.=FALSE)
  accession <- ram_uniprot_accession(page$accession)
  start <- ram_atlas_page_integer(page$start,"page offset")
  total <- ram_atlas_page_integer(page$total_count,"reported total")
  returned <- ram_atlas_page_integer(page$returned_count,"hit count")
  failed <- ram_atlas_page_integer(page$incomplete_metadata,
                                   "missing-metadata count")
  if(failed>returned) stop("Missing metadata exceeds returned hits.",call.=FALSE)
  records <- if(is.null(page$results)) list() else page$results
  if(!is.list(records) || length(records)>returned)
    stop("Atlas page has inconsistent enriched entity records.",call.=FALSE)

  if(start==0L) {
    previous <- NULL
  } else {
    if(!is.list(previous) ||
       !identical(as.character(previous$accession),as.character(accession)))
      stop("Atlas page does not match the active UniProt query.",call.=FALSE)
    if(!identical(as.integer(previous$next_offset),as.integer(start)))
      stop("Atlas page offset does not match the previous page.",call.=FALSE)
    if(!identical(as.integer(previous$total_count),as.integer(total)))
      stop("RCSB total changed between pages. Restart the Atlas search.",
           call.=FALSE)
    if(isFALSE(previous$has_more))
      stop("Atlas query has no more pages to fetch.",call.=FALSE)
  }

  all <- c(if(is.null(previous)) list() else previous$results,records)
  ids <- vapply(all,ram_atlas_record_key,character(1L))
  keep <- !duplicated(ids)
  unique_records <- all[keep]
  offset <- start+returned
  if(offset > total)
    stop("Atlas page returned more hits than the declared total.",call.=FALSE)
  stalled <- returned==0L && offset<total
  list(
    accession=accession,
    total_count=total,
    next_offset=offset,
    returned_count=offset,
    enriched_count=length(unique_records),
    incomplete_metadata=(if(is.null(previous)) 0L
       else previous$incomplete_metadata)+failed,
    duplicate_count=(if(is.null(previous)) 0L
      else previous$duplicate_count)+sum(!keep),
    pages=(if(is.null(previous)) 0L else previous$pages)+1L,
    has_more=offset<total && !stalled,
    stalled=stalled,
    results=unique_records
  )
}
