# Experimental provenance and confounder context for verified Atlas groups.
# Reuses existing RCSB discovery metadata and the exact observed monomer audit.
# Does not infer ligand occupancy or biological/functional state labels.

ram_atlas_experimental_context <- function(cohort,geometry,audit) {
  if(!is.list(cohort) || !is.list(cohort$results) ||
     !is.list(geometry) || !is.data.frame(geometry$assignment) ||
     !is.list(audit) || !is.data.frame(audit$pairs))
    stop("Experimental context requires a discovered, verified cohort.",
      call.=FALSE)
  ids <- as.character(geometry$selected)
  assignment <- geometry$assignment
  if(length(ids)<2L || length(ids)>12L ||
     anyNA(ids) || anyDuplicated(ids) ||
     !all(c("entity","geometric_group","common_coverage") %in%
          names(assignment)) ||
     !setequal(ids,as.character(assignment$entity)) ||
     !identical(as.character(cohort$accession),
                as.character(geometry$accession)) ||
     !setequal(ids,as.character(audit$selected)))
    stop("Discovery, geometry and construct review must use one cohort.",
      call.=FALSE)
  all_keys <- vapply(cohort$results,ram_atlas_record_key,character(1L))
  if(anyDuplicated(all_keys))
    stop("Duplicate RCSB metadata cannot be matched unambiguously.",
      call.=FALSE)
  index <- match(ids,all_keys)
  # Some verified structures may be missing RCSB enrichment. Do not
  # synthesize resolution, method or release date from other entries.
  scalar_text <- function(record,key) {
    value <- if(is.null(record)) NULL else record[[key]]
    if(is.null(value) || length(value)!=1L || is.na(value) ||
       !is.character(value) || !nzchar(trimws(value)))
      return(NA_character_)
    text <- trimws(value)
    if(key=="method" && identical(text,"Unknown method"))
      return(NA_character_)
    text
  }
  scalar_resolution <- function(record) {
    value <- if(is.null(record)) NULL else record$resolution
    if(length(value)!=1L) return(NA_real_)
    numeric <- suppressWarnings(as.numeric(value))
    if(!is.finite(numeric) || numeric<=0) NA_real_ else numeric
  }
  entry <- do.call(rbind,lapply(seq_along(ids),function(i) {
    id <- ids[[i]]
    record <- if(is.na(index[[i]])) NULL else cohort$results[[index[[i]]]]
    profile <- audit$profiles[[id]]
    known <- if(!is.null(profile)) as.integer(profile$known) else NA_integer_
    observed <- if(!is.null(profile)) as.integer(profile$observed) else NA_integer_
    data.frame(entity=id,group=as.character(
      assignment$geometric_group[match(id,assignment$entity)]),
      method=scalar_text(record,"method"),
      resolution_A=scalar_resolution(record),
      initial_release_date=scalar_text(record,"release_date"),
      description=scalar_text(record,"description"),
      observed_canonical_coverage=as.numeric(
        assignment$common_coverage[match(id,assignment$entity)]),
      monomer_known=known,monomer_observed=observed,
      source=if(is.null(record)) "RCSB enrichment unavailable"
        else "RCSB discovered polymer entity",
      stringsAsFactors=FALSE)
  }))
  rownames(entry) <- NULL

  required <- c("entity_a","entity_b","common_observed",
    "chemistry_differences","chemistry_unknown",
    "unmatched_a","unmatched_b")
  if(!all(required %in% names(audit$pairs)))
    stop("Incomplete exact monomer audit for experimental context.",
      call.=FALSE)
  pairs <- audit$pairs
  if(nrow(pairs)!=choose(length(ids),2L) ||
     any(!pairs$entity_a %in% ids) || any(!pairs$entity_b %in% ids))
    stop("Construct audit does not cover all selected experimental pairs.",
      call.=FALSE)
  a <- match(pairs$entity_a,entry$entity)
  b <- match(pairs$entity_b,entry$entity)
  pairs$same_geometric_group <- entry$group[a]==entry$group[b]
  pairs$same_pdb_entry <- substr(pairs$entity_a,1,4)==
    substr(pairs$entity_b,1,4)
  pairs$method_a <- entry$method[a]
  pairs$method_b <- entry$method[b]
  pairs$method_differs <- ifelse(
    !is.na(pairs$method_a) & !is.na(pairs$method_b),
    pairs$method_a!=pairs$method_b,NA)
  pairs$resolution_a_A <- entry$resolution_A[a]
  pairs$resolution_b_A <- entry$resolution_A[b]
  pairs$resolution_gap_A <- ifelse(
    is.finite(pairs$resolution_a_A) &
      is.finite(pairs$resolution_b_A),
    abs(pairs$resolution_a_A-pairs$resolution_b_A),NA_real_)
  pairs$release_a <- entry$initial_release_date[a]
  pairs$release_b <- entry$initial_release_date[b]
  # Metadata differences are annotations, not numerical quality penalties.
  pairs$context_warning <- vapply(seq_len(nrow(pairs)),function(i) {
    notes <- character()
    if(isTRUE(pairs$same_pdb_entry[[i]]))
      notes <- c(notes,"Same PDB deposition")
    if(isTRUE(pairs$method_differs[[i]]))
      notes <- c(notes,"Different experimental methods")
    if(is.na(pairs$method_differs[[i]]))
      notes <- c(notes,"Method comparison unavailable")
    if(pairs$chemistry_differences[[i]]>0L)
      notes <- c(notes,"Known monomer chemistry differs")
    if(pairs$chemistry_unknown[[i]]>0L)
      notes <- c(notes,"Unresolved monomer chemistry")
    if(pairs$unmatched_a[[i]]>0L || pairs$unmatched_b[[i]]>0L)
      notes <- c(notes,"Different observed spans")
    if(!length(notes)) "No listed metadata difference detected" else
      paste(notes,collapse="; ")
  },character(1L))
  list(entries=entry,pairs=pairs,
    missing_metadata=sum(is.na(index)),
    missing_methods=sum(is.na(entry$method)),
    missing_resolution=sum(!is.finite(entry$resolution_A)),
    shared_pdb_pairs=sum(pairs$same_pdb_entry),
    differing_method_pairs=sum(pairs$method_differs %in% TRUE),
    unresolved_method_pairs=sum(is.na(pairs$method_differs)),
    differing_chemistry_pairs=sum(pairs$chemistry_differences>0L),
    ligand_status="Ligand/cofactor occupancy has not been verified by Atlas; no apo/holo assignment is inferred.",
    method=paste("RCSB discovery method, resolution and initial release date;",
      "exact observed PDBe SIFTS/mmCIF polymer monomer compatibility.",
      "These are descriptive confounders, not biological-state labels",
      "or independent-experiment counts."))
}
