# Bridge between exact-verified experimental Atlas entities and group comparison.
# The Atlas provides first-model N/CA/C coordinates and exact SIFTS UniProt
# positions, so neither external downloads nor guessed sequence offsets are
# necessary. Comparisons remain descriptive, not functional-state labelling.

ram_atlas_group_select <- function(geometry,group_a,group_b) {
  if(!is.list(geometry) || !is.null(geometry$error) ||
     !is.data.frame(geometry$assignment) ||
     !length(geometry$common_positions))
    stop("Run a verified Atlas geometry comparison first.",call.=FALSE)
  a <- as.character(group_a); b <- as.character(group_b)
  ids <- as.character(geometry$selected)
  if(!length(a) || !length(b) || length(a)>12L || length(b)>12L ||
     anyNA(c(a,b)) || anyDuplicated(a) || anyDuplicated(b) ||
     length(intersect(a,b)) || any(!c(a,b) %in% ids))
    stop("Choose two disjoint, non-empty sets of verified Atlas entities.",
         call.=FALSE)
  list(group_a=a,group_b=b,
       selected=unique(c(a,b)),
       accession=geometry$accession,
       core=sort(unique(as.integer(geometry$common_positions))))
}

ram_atlas_group_prepare <- function(verified,geometry,group_a,group_b,
                                    label_a="Group A",label_b="Group B") {
  picked <- ram_atlas_group_select(geometry,group_a,group_b)
  if(length(picked$core)<5L)
    stop("Too few common verified UniProt positions.",call.=FALSE)
  available <- ram_atlas_geometry_entities(verified,picked$accession)
  if(any(!picked$selected %in% names(available)))
    stop("The selected experimental entities are no longer verified.",
         call.=FALSE)

  tables <- lapply(picked$selected,function(id) {
    entity <- available[[id]]
    raw <- ram_atlas_entity_torsions(verified[[id]],entity$struct_asym_id)
    if(!nrow(raw) || anyDuplicated(raw$uniprot_resi))
      stop(paste("No unique exact-mapped backbone torsions for",id),
           call.=FALSE)
    raw <- raw[raw$uniprot_resi %in% picked$core,,drop=FALSE]
    if(sum(is.finite(raw$phi) & is.finite(raw$psi))<5L)
      stop(paste("Insufficient complete N/CA/C torsions for",id),
           call.=FALSE)
    raw
  })
  names(tables) <- picked$selected

  # A residue name is displayed only when the selected models agree on
  # verified mmCIF monomer identity. Otherwise use an explicit unknown.
  # No Rama8000 labels are generated from C-alpha/backbone-only metadata.
  profiles <- lapply(picked$selected,function(id)
    ram_atlas_construct_profile(verified[[id]],available[[id]],
                                picked$accession))
  chemistry <- rep("UNK",length(picked$core))
  for(i in seq_along(picked$core)) {
    position <- picked$core[[i]]
    codes <- vapply(profiles,function(profile) {
      loc <- match(position,profile$positions)
      if(is.na(loc)) NA_character_
      else profile$chemistry[[loc]]
    },character(1L))
    if(all(!is.na(codes)) && length(unique(codes))==1L)
      chemistry[[i]] <- codes[[1L]]
  }
  names(chemistry) <- as.character(picked$core)
  models <- lapply(picked$selected,function(id) {
    tbl <- tables[[id]]
    data.frame(chain="UniProt",resi=as.integer(tbl$uniprot_resi),
      insertion_code="",resn=unname(chemistry[
        as.character(tbl$uniprot_resi)]),
      phi=as.numeric(tbl$phi),psi=as.numeric(tbl$psi),
      source_resn=unname(profiles[[match(id,picked$selected)]]$chemistry[
        match(tbl$uniprot_resi,
          profiles[[match(id,picked$selected)]]$positions)]),
      region=NA_character_,rama8000_region=NA_character_,
      stringsAsFactors=FALSE)
  })
  names(models) <- picked$selected
  names_a <- picked$group_a; names_b <- picked$group_b
  label_a <- trimws(as.character(label_a)); label_b <- trimws(as.character(label_b))
  if(length(label_a)!=1L || length(label_b)!=1L ||
     is.na(label_a) || is.na(label_b) ||
     !nzchar(label_a) || !nzchar(label_b))
    stop("Group labels cannot be empty.",call.=FALSE)
  comparison <- ram_group_conformation_compare(
    models[[names_a[[1L]]]],models[names_a],models[names_b],
    label_a,label_b)
  fingerprint <- ram_group_fingerprint(models[names_a],models[names_b],
                                        label_a,label_b)
  metadata <- do.call(rbind,lapply(picked$selected,function(id) {
    entity <- available[[id]]
    table <- models[[id]]
    data.frame(model=id,chain=entity$chain,
      group=if(id %in% names_a) label_a else label_b,
      identity=NA_real_,
      reference_coverage=nrow(table)/length(picked$core),
      candidate_coverage=nrow(table)/nrow(tables[[id]]),
      aligned=nrow(table),
      complete_torsions=sum(is.finite(table$phi) &
                            is.finite(table$psi)),
      mapping="Exact observed PDBe SIFTS UniProt",
      stringsAsFactors=FALSE)
  }))
  audit <- ram_atlas_construct_audit(verified,picked$accession,
                                    picked$selected,available)
  ligand_context <- ram_atlas_observed_ligand_context(verified,geometry)
  list(comparison=comparison,fingerprint=fingerprint,members=metadata,
    ligand_context=ligand_context,
    label_a=label_a,label_b=label_b,n_a=length(names_a),
    n_b=length(names_b),reference_chain="UniProt",
    source="atlas",accession=picked$accession,
    core_positions=length(picked$core),
    unknown_chemistry_positions=sum(chemistry=="UNK"),
    known_chemistry_differences=audit$differs,
    method="Exact observed PDBe SIFTS UniProt backbone torsions; no Rama8000 classification or claims of biological states")
}
