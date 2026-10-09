# Conservative experimental-construct compatibility audit for Atlas grouping.
# The residue chemistry is the actual _pdbx_poly_seq_scheme.mon_id linked to
# EXACT observed PDBe SIFTS UniProt positions. "MET" vs "MSE" is a
# chemistry difference, not automatically a gene mutation. Missing monomer
# metadata or different observed spans cannot establish construct equivalence.
#
# This audit does not label structures as biological states or independent
# experimental replicates. It only surfaces specific confounding evidence.

ram_atlas_construct_profile <- function(record, entity, accession) {
  accession <- ram_uniprot_accession(accession)
  positions <- as.integer(entity$coordinates$uniprot_resi)
  stopifnot(!anyNA(positions), !anyDuplicated(positions))
  chemicals <- rep(NA_character_,length(positions))
  names(chemicals) <- as.character(positions)
  map <- record$mapping
  if(is.data.frame(map) && "mon_id" %in% names(map) &&
     all(c("struct_asym_id","uniprot_accession","uniprot_resi",
           "observed") %in% names(map))) {
    relevant <- !is.na(map$struct_asym_id) &
      as.character(map$struct_asym_id)==entity$struct_asym_id &
      !is.na(map$uniprot_accession) &
      as.character(map$uniprot_accession)==accession &
      !is.na(map$observed) & map$observed &
      is.finite(map$uniprot_resi)
    sub <- map[relevant,,drop=FALSE]
    if(nrow(sub)) {
      for(i in seq_along(positions)) {
        codes <- unique(as.character(sub$mon_id[
          sub$uniprot_resi==positions[[i]]]))
        codes <- codes[!is.na(codes) & nzchar(codes) &
          !codes %in% c("UNK","UNL","X")]
        # Conflicting / multiple identities are unknown, never first-wins.
        if(length(codes)==1L &&
           !any(is.na(sub$mon_id[sub$uniprot_resi==positions[[i]]])))
          chemicals[[i]] <- toupper(codes[[1L]])
      }
    }
  }
  list(positions=positions,chemistry=chemicals,
       known=sum(!is.na(chemicals)),observed=length(positions))
}

ram_atlas_construct_audit <- function(verified,accession,selected,
                                      available=NULL) {
  ids <- unique(as.character(selected))
  if(length(ids)<2L || length(ids)>12L || anyNA(ids) ||
     any(!grepl("^[A-Z0-9]{4}_[1-9][0-9]*$",ids)))
    stop("Select 2 to 12 verified experimental entities for construct review.",
         call.=FALSE)
  if(is.null(available))
    available <- ram_atlas_geometry_entities(verified,accession)
  if(any(!ids %in% names(available)))
    stop("Missing verified geometry for construct review.",call.=FALSE)
  profiles <- lapply(ids,function(id)
    ram_atlas_construct_profile(verified[[id]],available[[id]],accession))
  names(profiles) <- ids
  rows <- list()
  for(i in seq_len(length(ids)-1L)) for(j in seq.int(i+1L,length(ids))) {
    a <- profiles[[i]]
    b <- profiles[[j]]
    shared <- sort(intersect(a$positions,b$positions))
    first <- a$chemistry[match(shared,a$positions)]
    second <- b$chemistry[match(shared,b$positions)]
    known <- !is.na(first) & !is.na(second)
    differs <- known & first!=second
    changes <- which(differs)
    rows[[length(rows)+1L]] <- data.frame(
      entity_a=ids[[i]],entity_b=ids[[j]],
      common_observed=length(shared),
      coverage_a=if(a$observed) round(length(shared)/a$observed,3)
        else NA_real_,
      coverage_b=if(b$observed) round(length(shared)/b$observed,3)
        else NA_real_,
      unmatched_a=a$observed-length(shared),
      unmatched_b=b$observed-length(shared),
      chemistry_checked=sum(known),
      chemistry_unknown=sum(!known),
      chemistry_differences=sum(differs),
      chemistry_agreement=if(any(known))
        round(sum(known & !differs)/sum(known),3) else NA_real_,
      difference_examples=if(length(changes))
        paste(head(sprintf("%d:%s/%s",shared[changes],
          first[changes],second[changes]),12L),collapse=", ") else "",
      stringsAsFactors=FALSE
    )
  }
  pairs <- do.call(rbind,rows)
  rownames(pairs) <- NULL
  list(accession=ram_uniprot_accession(accession),
       selected=ids,profiles=profiles,pairs=pairs,
       pair_count=nrow(pairs),
       differs=any(pairs$chemistry_differences>0L),
       incomplete=any(pairs$chemistry_unknown>0L) ||
         any(pairs$unmatched_a>0L | pairs$unmatched_b>0L),
       all_chemistry_known=all(pairs$chemistry_unknown==0L),
       method=paste("Observed first-model C-alpha residues with exact",
         "PDBe SIFTS UniProt coordinate and mmCIF polymer monomer IDs; ",
         "pairwise shared-core chemistry, not biological-state validation"))
}
