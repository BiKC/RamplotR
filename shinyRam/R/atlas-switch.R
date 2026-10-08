# Atlas local-backbone comparison on exact SIFTS/UniProt coordinates.
# Reuses ram_dihedral() from backbone.R and the geometry group's selected
# structure chains. No residue-number interpolation, secondary-structure
# inference or statistical/functional state assignment.

ram_atlas_backbone_records <- function(records) {
  empty <- data.frame(struct_asym_id=character(),label_seq_id=integer(),
    atom_name=character(),x=numeric(),y=numeric(),z=numeric())
  if(is.null(records) || !length(records)) return(empty)
  if(is.data.frame(records)) out <- records else {
    if(!is.list(records) || length(records)>15000L)
      stop("Invalid Atlas backbone atom records.",call.=FALSE)
    out <- do.call(rbind,lapply(records,function(row) {
      if(!is.list(row)) stop("Malformed backbone atom.",call.=FALSE)
      data.frame(
        struct_asym_id=as.character(row$struct_asym_id),
        label_seq_id=as.integer(row$label_seq_id),
        atom_name=as.character(row$atom_name),
        x=as.numeric(row$x),y=as.numeric(row$y),z=as.numeric(row$z))
    }))
  }
  required <- names(empty)
  if(!all(required %in% names(out)))
    stop("Backbone atom fields are incomplete.",call.=FALSE)
  out <- out[,required,drop=FALSE]
  if(anyNA(out) || any(!is.finite(as.matrix(out[,c("x","y","z")]))) ||
     any(!out$atom_name %in% c("N","CA","C")) ||
     any(out$label_seq_id<1L))
    stop("Invalid exact-mapped backbone atoms.",call.=FALSE)
  out
}

ram_atlas_entity_torsions <- function(record,asym_id) {
  empty <- data.frame(uniprot_resi=integer(),phi=numeric(),
    psi=numeric(),struct_asym_id=character(),label_seq_id=integer(),
    stringsAsFactors=FALSE)
  if(!is.list(record) || !identical(record$state,"mapped") ||
     !is.data.frame(record$mapping) || !nrow(record$mapping))
    return(empty)
  atoms <- ram_atlas_backbone_records(record$backbone_atoms)
  if(!nrow(atoms)) return(empty)
  asym_id <- as.character(asym_id)
  if(length(asym_id)!=1L || is.na(asym_id) || !nzchar(asym_id))
    stop("Choose one exact-mapped polymer chain.",call.=FALSE)
  atoms <- atoms[atoms$struct_asym_id==asym_id,,drop=FALSE]
  if(!nrow(atoms)) return(empty)
  atom_key <- paste(atoms$label_seq_id,atoms$atom_name,sep=":")
  # Ambiguous atom records cannot be used to construct an angle.
  atoms <- atoms[!duplicated(atom_key) &
    !duplicated(atom_key,fromLast=TRUE),,drop=FALSE]
  if(!nrow(atoms)) return(empty)
  seqs <- sort(unique(atoms$label_seq_id))
  coord <- function(atom) {
    rows <- atoms[atoms$atom_name==atom,,drop=FALSE]
    idx <- match(seqs,rows$label_seq_id)
    as.matrix(rows[idx,c("x","y","z"),drop=FALSE])
  }
  nxyz <- coord("N"); caxyz <- coord("CA"); cxyz <- coord("C")
  good <- stats::complete.cases(cbind(nxyz,caxyz,cxyz))
  connected <- rep(FALSE,length(seqs))
  if(length(seqs)>1L) for(i in seq_len(length(seqs)-1L)) {
    d <- sqrt(sum((cxyz[i,]-nxyz[i+1L,])^2))
    connected[i] <- good[i] && good[i+1L] &&
      seqs[i+1L]==seqs[i]+1L && is.finite(d) && d>=1.0 && d<=1.9
  }
  phi <- psi <- rep(NA_real_,length(seqs))
  if(length(seqs)>1L) for(i in seq_along(seqs)) {
    if(i>1L && connected[i-1L])
      phi[i] <- ram_dihedral(cxyz[i-1L,],nxyz[i,],caxyz[i,],cxyz[i,])
    if(i<length(seqs) && connected[i])
      psi[i] <- ram_dihedral(nxyz[i,],caxyz[i,],cxyz[i,],nxyz[i+1L,])
  }
  mapped <- record$mapping
  req <- c("struct_asym_id","label_seq_id","uniprot_resi","observed")
  if(!all(req %in% names(mapped))) return(empty)
  mapped <- mapped[mapped$struct_asym_id==asym_id &
    !is.na(mapped$observed) & mapped$observed &
    is.finite(mapped$label_seq_id) &
    is.finite(mapped$uniprot_resi),,drop=FALSE]
  if(!nrow(mapped)) return(empty)
  mapped <- unique(mapped[,c("label_seq_id","uniprot_resi","chain",
    "resi","insertion_code"),drop=FALSE])
  # Both label->UniProt and UniProt->label assignments must be one-to-one.
  bad <- duplicated(mapped$label_seq_id) |
    duplicated(mapped$label_seq_id,fromLast=TRUE) |
    duplicated(mapped$uniprot_resi) |
    duplicated(mapped$uniprot_resi,fromLast=TRUE)
  mapped <- mapped[!bad,,drop=FALSE]
  matched <- match(mapped$label_seq_id,seqs)
  mapped <- mapped[!is.na(matched),,drop=FALSE]
  matched <- matched[!is.na(matched)]
  if(!nrow(mapped)) return(empty)
  result <- data.frame(uniprot_resi=as.integer(mapped$uniprot_resi),
    phi=phi[matched],psi=psi[matched],
    struct_asym_id=asym_id,
    label_seq_id=as.integer(mapped$label_seq_id),
    chain=as.character(mapped$chain),resi=as.integer(mapped$resi),
    insertion_code=as.character(mapped$insertion_code),
    stringsAsFactors=FALSE)
  result[order(result$uniprot_resi),,drop=FALSE]
}

ram_atlas_torsion_delta <- function(a,b,id_a,id_b,threshold=30) {
  threshold <- as.numeric(threshold)
  if(length(threshold)!=1L || !is.finite(threshold) ||
     threshold<5 || threshold>180)
    stop("Backbone change threshold must be between 5 and 180 degrees.",
         call.=FALSE)
  if(!all(c("uniprot_resi","phi","psi") %in% names(a)) ||
     !all(c("uniprot_resi","phi","psi") %in% names(b)))
    stop("Invalid canonical backbone torsion tables.",call.=FALSE)
  if(anyDuplicated(a$uniprot_resi) || anyDuplicated(b$uniprot_resi))
    stop("Backbone mapping has ambiguous UniProt positions.",call.=FALSE)
  pos <- sort(unique(c(a$uniprot_resi,b$uniprot_resi)))
  ia <- match(pos,a$uniprot_resi)
  ib <- match(pos,b$uniprot_resi)
  fetch <- function(table,index,column,missing=NA_real_) {
    if(column %in% names(table)) table[[column]][index]
    else rep(missing,length(index))
  }
  pa <- fetch(a,ia,"phi"); qa <- fetch(a,ia,"psi")
  pb <- fetch(b,ib,"phi"); qb <- fetch(b,ib,"psi")
  wrap <- function(a,b) {
    x <- (b-a+180) %% 360-180
    x[!is.finite(a) | !is.finite(b)] <- NA_real_
    x
  }
  dp <- wrap(pa,pb); dq <- wrap(qa,qb)
  change <- sqrt(dp^2+dq^2)
  data.frame(uniprot_resi=pos,entity_a=id_a,entity_b=id_b,
    chain_a=fetch(a,ia,"chain",NA_character_),
    resi_a=as.integer(fetch(a,ia,"resi")),
    insertion_a=fetch(a,ia,"insertion_code",NA_character_),
    label_seq_a=as.integer(fetch(a,ia,"label_seq_id")),
    chain_b=fetch(b,ib,"chain",NA_character_),
    resi_b=as.integer(fetch(b,ib,"resi")),
    insertion_b=fetch(b,ib,"insertion_code",NA_character_),
    label_seq_b=as.integer(fetch(b,ib,"label_seq_id")),
    phi_a=pa,psi_a=qa,phi_b=pb,psi_b=qb,
    delta_phi=dp,delta_psi=dq,angular_shift=change,
    comparable=is.finite(change),
    candidate=is.finite(change) & change>=threshold,
    stringsAsFactors=FALSE)
}

ram_atlas_candidate_regions <- function(differences) {
  empty <- data.frame(start=integer(),end=integer(),residues=integer(),
    mean_shift=numeric(),peak_shift=numeric(),
    stringsAsFactors=FALSE)
  if(!nrow(differences)) return(empty)
  if(!all(c("uniprot_resi","candidate","angular_shift") %in%
          names(differences)))
    stop("Invalid Atlas torsion differences.",call.=FALSE)
  data <- differences[order(differences$uniprot_resi),,drop=FALSE]
  flags <- which(!is.na(data$candidate) & data$candidate &
                   is.finite(data$angular_shift))
  if(!length(flags))return(empty)
  segments <- split(flags,cumsum(c(TRUE,diff(flags)!=1L |
      diff(data$uniprot_resi[flags])!=1L)))
  do.call(rbind,lapply(segments,function(idx) data.frame(
    start=min(data$uniprot_resi[idx]),end=max(data$uniprot_resi[idx]),
    residues=length(idx),mean_shift=round(mean(data$angular_shift[idx]),1),
    peak_shift=round(max(data$angular_shift[idx]),1))))
}

ram_atlas_group_switches <- function(verified,geometry,threshold=30,
  representative_ids=NULL) {
  if(!is.list(geometry) || is.null(geometry$representatives) ||
     is.null(geometry$assignment))
    stop("Run an Atlas experimental geometry comparison first.",
         call.=FALSE)
  representatives <- unname(geometry$representatives)
  if(length(representatives)<2L)
    stop("Only one geometric group; no between-group representative comparison.",
         call.=FALSE)
  # Any two distinct geometric-group medoids may be compared, without
  # silently treating them as independent biological states.
  ids <- if(is.null(representative_ids)) representatives[1:2] else
    as.character(representative_ids)
  if(length(ids)!=2L || anyNA(ids) || any(!ids %in% representatives) ||
     identical(ids[[1L]],ids[[2L]]))
    stop("Choose two distinct experimental geometry-group representatives.",
         call.=FALSE)
  available <- ram_atlas_geometry_entities(verified,geometry$accession)
  tables <- lapply(ids,function(id) {
    structure <- available[[id]]
    if(is.null(structure)) stop("Missing verified representative.",call.=FALSE)
    ram_atlas_entity_torsions(verified[[id]],structure$struct_asym_id)
  })
  diff <- ram_atlas_torsion_delta(tables[[1L]],tables[[2L]],
    ids[[1L]],ids[[2L]],threshold)
  list(representatives=ids,threshold=as.numeric(threshold),
    residues=diff,regions=ram_atlas_candidate_regions(diff),
    comparable=sum(diff$comparable),total_positions=nrow(diff),
    method="Paired exact UniProt residue phi/psi, circular wrapped angle differences; peptide bond 1.0–1.9 Å; first model, two geometric group representatives")
}
