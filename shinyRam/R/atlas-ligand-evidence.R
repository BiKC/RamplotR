# Explicitly observed non-water nonpolymer components near exact-mapped
# first-model protein heavy atoms from the SAME PDBe updated mmCIF.
# Proximity does not establish occupancy, biochemical binding or apo/holo states.
ram_atlas_observed_ligand_context <- function(verified,geometry) {
  ids <- as.character(geometry$selected)
  if(!is.list(geometry) || length(ids)<2L || length(ids)>12L ||
     anyNA(ids) || anyDuplicated(ids) || !is.list(verified) ||
     any(!ids %in% names(verified)))
    stop("Ligand proximity requires verified, selected Atlas structures.",
      call.=FALSE)
  rows <- list()
  contacts <- list()
  for(id in ids) {
    item <- verified[[id]]
    raw <- item$ligand_contacts
    status <- if(is.list(raw) && identical(raw$status,"measured"))
      "measured" else "unavailable"
    warning <- if(is.list(raw) && is.character(raw$warning) &&
      length(raw$warning)==1L) raw$warning else
      "Ligand extraction was not completed; no absence can be inferred."
    sites <- if(identical(status,"measured") && is.list(raw$sites))
      raw$sites else list()
    if(length(sites)>250L) {
      status <- "unavailable";sites <- list()
      warning <- "Ligand contact results exceeded the report limit."
    }
    parsed <- list()
    if(identical(status,"measured") && length(sites)) {
      parsed <- tryCatch(lapply(sites,function(site) {
        if(!is.list(site)) stop("Invalid deposited component entry.")
        comp <- toupper(as.character(site$comp_id))
        asym <- as.character(site$asym_id)
        auth <- as.character(site$auth_seq_id)
        distance <- suppressWarnings(as.numeric(site$min_distance_A))
        pos <- suppressWarnings(as.integer(site$nearest_uniprot_resi))
        count <- suppressWarnings(as.integer(site$heavy_atoms))
        if(length(comp)!=1L || is.na(comp) ||
           !grepl("^[A-Z0-9-]{1,8}$",comp) ||
           length(asym)!=1L || is.na(asym) ||
           length(auth)!=1L || is.na(auth) ||
           length(distance)!=1L || !is.finite(distance) ||
           distance<0 || distance>4.5 ||
           length(pos)!=1L || is.na(pos) || pos<1L ||
           length(count)!=1L || is.na(count) || count<1L)
          stop("Invalid observed ligand proximity fields.")
        data.frame(entity=id,comp_id=comp,asym_id=asym,
          auth_seq_id=auth,nearest_uniprot_resi=pos,
          min_distance_A=round(distance,3),heavy_atoms=count,
          stringsAsFactors=FALSE)
      }),error=function(e) e)
      if(inherits(parsed,"error")) {
        warning <- conditionMessage(parsed)
        status <- "unavailable";parsed <- list()
      }
    }
    if(length(parsed)) contacts[[length(contacts)+1L]] <- do.call(rbind,parsed)
    reported <- if(status=="measured") length(parsed) else NA_integer_
    total <- if(status=="measured" && !is.null(raw$total_nonwater_sites))
      suppressWarnings(as.integer(raw$total_nonwater_sites)) else NA_integer_
    if(length(total)!=1L || is.na(total) || total<reported) {
      if(status=="measured") {
        status <- "unavailable";reported <- NA_integer_
        warning <- "Invalid deposited nonpolymer component count."
        if(length(parsed)) contacts[[length(contacts)]]<-NULL
      }
      total <- NA_integer_
    }
    rows[[length(rows)+1L]] <- data.frame(
      entity=id,geometry_group=as.character(
        geometry$assignment$geometric_group[
          match(id,geometry$assignment$entity)]),
      status=status,nearby_nonwater_sites=reported,
      all_nonwater_sites=total,radius_A=4.5,
      warning=warning,stringsAsFactors=FALSE)
  }
  summary <- do.call(rbind,rows)
  hits <- if(length(contacts)) do.call(rbind,contacts) else data.frame(
    entity=character(),comp_id=character(),asym_id=character(),
    auth_seq_id=character(),nearest_uniprot_resi=integer(),
    min_distance_A=numeric(),heavy_atoms=integer(),
    stringsAsFactors=FALSE)
  rownames(summary) <- NULL; rownames(hits) <- NULL
  list(entries=summary,contacts=hits,
    measured=sum(summary$status=="measured"),
    unavailable=sum(summary$status!="measured"),
    with_proximity=sum(summary$nearby_nonwater_sites>0L,na.rm=TRUE),
    method=paste(
      "PDBe updated mmCIF first-model nonwater nonpolymer HETATM heavy atoms",
      "within 4.5 Angstrom of the selected chain's exact SIFTS-mapped",
      "protein heavy atoms. This is geometric proximity only, not",
      "ligand binding, verified occupancy, apo/holo or functional state.",
      "No local contact does not prove no ligand or an apo structure."))
}
