# Fragment-level local backbone fingerprints using Protein Blocks (PBs).
#
# Protein Blocks are 16 pentapeptide prototypes (a-p) defined by eight
# consecutive backbone dihedral angles. The reference values below are from:
# de Brevern AG, Etchebest C, Hazout S. Proteins. 2000;41:271-288.
#
# The assignment matches the public PBxplore convention:
# psi(n-2), phi(n-1), psi(n-1), phi(n), psi(n), phi(n+1),
# psi(n+1), phi(n+2), with periodic angular differences.
#
# PBs are a structural alphabet, not a validation score or a secondary-
# structure assignment. See static/protein-blocks/SOURCE.md for provenance.

ram_protein_block_references <- rbind(
  a=c( 41.14,  75.53,  13.92, -99.80, 131.88, -96.27, 122.08, -99.68),
  b=c(108.24, -90.12, 119.54, -92.21, -18.06,-128.93, 147.04, -99.90),
  c=c(-11.61,-105.66,  94.81,-106.09, 133.56,-106.93, 135.97,-100.63),
  d=c(141.98,-112.79, 132.20,-114.79, 140.11,-111.05, 139.54,-103.16),
  e=c(133.25,-112.37, 137.64,-108.13, 133.00, -87.30, 120.54,  77.40),
  f=c(116.40,-105.53, 129.32, -96.68, 140.72, -74.19, -26.65, -94.51),
  g=c(  0.40, -81.83,   4.91,-100.59,  85.50, -71.65, 130.78,  84.98),
  h=c(119.14,-102.58, 130.83, -67.91, 121.55,  76.25,  -2.95, -90.88),
  i=c(130.68, -56.92, 119.26,  77.85,  10.42, -99.43, 141.40, -98.01),
  j=c(114.32,-121.47, 118.14,  82.88,-150.05, -83.81,  23.35, -85.82),
  k=c(117.16, -95.41, 140.40, -59.35, -29.23, -72.39, -25.08, -76.16),
  l=c(139.20, -55.96, -32.70, -68.51, -26.09, -74.44, -22.60, -71.74),
  m=c(-39.62, -64.73, -39.52, -65.54, -38.88, -66.89, -37.76, -70.19),
  n=c(-35.34, -65.03, -38.12, -66.34, -29.51, -89.10,  -2.91,  77.90),
  o=c(-45.29, -67.44, -27.72, -87.27,   5.13,  77.49,  30.71, -93.23),
  p=c(-27.09, -86.14,   0.30,  59.85,  21.51, -96.30, 132.67, -92.91)
)
colnames(ram_protein_block_references) <- c(
  "psi_m2","phi_m1","psi_m1","phi_0",
  "psi_0","phi_p1","psi_p1","phi_p2"
)

ram_pb_angle_difference <- function(reference,observed) {
  ((as.numeric(reference)-as.numeric(observed)+180) %% 360)-180
}

ram_pb_distance <- function(angles,
                            references=ram_protein_block_references) {
  angles <- as.numeric(angles)
  if(length(angles)!=8L || any(!is.finite(angles)))
    return(rep(NA_real_,nrow(references)))
  if(!is.matrix(references) || ncol(references)!=8L)
    stop("Protein Block reference matrix must contain eight angles.",
         call.=FALSE)
  sqrt(rowMeans(vapply(seq_len(nrow(references)),function(i)
    ram_pb_angle_difference(references[i,],angles)^2,
    numeric(8L))))
}

ram_pb_assign_angles <- function(angles,
                                 references=ram_protein_block_references) {
  distance <- ram_pb_distance(angles,references)
  if(!length(distance) || all(!is.finite(distance)))
    return(c(block=NA_character_,rmsda=NA_character_))
  i <- which.min(distance)
  c(block=rownames(references)[[i]],rmsda=as.character(distance[[i]]))
}

ram_protein_blocks <- function(data,
                               references=ram_protein_block_references) {
  required <- c("chain","resi","insertion_code","phi","psi")
  if(!is.data.frame(data) || !all(required %in% names(data)))
    stop("Protein Blocks require chain/residue identifiers and phi/psi.",
         call.=FALSE)
  n <- nrow(data)
  out <- data.frame(
    protein_block=rep(NA_character_,n),
    protein_block_rmsda=rep(NA_real_,n),
    protein_block_complete=rep(FALSE,n),
    stringsAsFactors=FALSE
  )
  if(n<5L) return(cbind(data,out))

  chains <- as.character(data$chain)
  has_bond <- if("bonded_to_next" %in% names(data))
    !is.na(data$bonded_to_next) & data$bonded_to_next
  else {
    # Without explicit peptide-continuity bookkeeping, only use adjacent rows
    # in the same chain. This is acceptable for already validated aligned
    # tables but raw structures should carry bonded_to_next.
    c(chains[-n]==chains[-1L],FALSE)
  }

  for(i in 3L:(n-2L)) {
    window <- (i-2L):(i+2L)
    if(length(unique(chains[window]))!=1L) next
    if(!all(has_bond[(i-2L):(i+1L)])) next
    angles <- c(
      data$psi[[i-2L]],
      data$phi[[i-1L]],data$psi[[i-1L]],
      data$phi[[i]],data$psi[[i]],
      data$phi[[i+1L]],data$psi[[i+1L]],
      data$phi[[i+2L]]
    )
    if(any(!is.finite(angles))) next
    distances <- ram_pb_distance(angles,references)
    j <- which.min(distances)
    out$protein_block[[i]] <- rownames(references)[[j]]
    out$protein_block_rmsda[[i]] <- distances[[j]]
    out$protein_block_complete[[i]] <- TRUE
  }
  cbind(data,out)
}

ram_pb_frequencies <- function(blocks) {
  values <- as.character(blocks)
  values <- values[!is.na(values) & values %in%
                   rownames(ram_protein_block_references)]
  counts <- table(factor(values,
    levels=rownames(ram_protein_block_references)))
  if(!sum(counts)) return(stats::setNames(rep(0,16L),rownames(
    ram_protein_block_references)))
  as.numeric(counts)/sum(counts) |>
    stats::setNames(names(counts))
}

ram_pb_equivalent_number <- function(blocks) {
  p <- ram_pb_frequencies(blocks)
  p <- p[p>0]
  if(!length(p)) return(NA_real_)
  exp(-sum(p*log(p)))
}

ram_pb_prototype_coarse_benchmark <- function(
    references=ram_protein_block_references) {
  if(!exists("ram_backbone_basin",mode="function"))
    stop("Source conformation.R before benchmarking coarse states.",
         call.=FALSE)
  if(!exists("ram_backbone_basin_centers"))
    stop("Backbone basin centres are unavailable.",call.=FALSE)

  central_phi <- references[,"phi_0"]
  central_psi <- references[,"psi_0"]
  wrapped <- function(value,center)
    ((value-center+180) %% 360)-180
  distances <- vapply(seq_len(nrow(ram_backbone_basin_centers)),function(i) {
    dphi <- wrapped(central_phi,ram_backbone_basin_centers$phi[[i]])
    dpsi <- wrapped(central_psi,ram_backbone_basin_centers$psi[[i]])
    sqrt(dphi^2+dpsi^2)
  },numeric(length(central_phi)))
  if(length(central_phi)==1L) distances <- matrix(distances,nrow=1L)
  ordered <- t(apply(distances,1L,sort))
  data.frame(
    protein_block=rownames(references),
    central_phi=as.numeric(central_phi),
    central_psi=as.numeric(central_psi),
    coarse_state=ram_backbone_basin(central_phi,central_psi),
    nearest_coarse_distance=ordered[,1L],
    coarse_margin=ordered[,2L]-ordered[,1L],
    stringsAsFactors=FALSE
  )
}
