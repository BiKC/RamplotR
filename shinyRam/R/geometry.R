# Additional *diagnostic* geometry for experimental and predicted proteins.
# Reference-grade clashes and rotamer outliers come from an independent
# validation report, NOT these descriptive geometry measurements.

ram_geom_omega_status <- function(degrees, tolerance = 30) {
  out <- rep("Missing", length(degrees))
  valid <- is.finite(degrees)
  out[valid] <- "Twisted"
  out[valid & abs(degrees) <= tolerance] <- "Cis"
  out[valid & abs(abs(degrees) - 180) <= tolerance] <- "Trans"
  out
}

ram_geom_chi1_atom <- c(
  ARG="CG", ASN="CG", ASP="CG", CYS="SG", GLN="CG", GLU="CG",
  HIS="CG", ILE="CG1", LEU="CG", LYS="CG", MET="CG", PHE="CG",
  PRO="CG", SER="OG", THR="OG1", TRP="CG", TYR="CG", VAL="CG1"
)

ram_extra_geometry <- function(pdb, torsions, peptide_min = 1.0,
                               peptide_max = 1.9) {
  required <- c("chain","resno","resid","elety","x","y","z")
  atoms <- pdb$atom
  if (!is.data.frame(atoms) || !all(required %in% names(atoms)))
    stop("A Bio3D structure with named atom coordinates is required.")
  n <- nrow(torsions)
  result <- data.frame(
    omega=rep(NA_real_,n), peptide_bond_length=rep(NA_real_,n),
    omega_status=rep("Missing",n), chi1=rep(NA_real_,n),
    chi1_available=rep(FALSE,n),
    cb_ca_distance=rep(NA_real_,n), cb_signed_volume=rep(NA_real_,n),
    stringsAsFactors=FALSE)
  if (!n) return(result)
  if (!"insert" %in% names(atoms)) atoms$insert <- ""
  if (!"alt" %in% names(atoms)) atoms$alt <- ""
  atoms$insert[is.na(atoms$insert)] <- ""
  atoms$alt[is.na(atoms$alt)] <- ""
  atom_keys <- paste(atoms$chain, atoms$resno, atoms$insert, sep="\r")
  residue_keys <- paste(torsions$chain, torsions$resi,
                        torsions$insertion_code, sep="\r")
  # The original atom table may contain repeated atom names for alternate
  # conformers. Unlabelled atoms take priority, then alternate A.
  valid <- atoms$alt %in% c("", "A") & is.finite(atoms$x) &
           is.finite(atoms$y) & is.finite(atoms$z)
  atoms <- atoms[valid, , drop=FALSE]
  atom_keys <- atom_keys[valid]
  atom_keys_full <- paste(atom_keys, atoms$elety, sep="\r")
  ord <- order(atoms$alt != "", seq_len(nrow(atoms)))
  first <- ord[!duplicated(atom_keys_full[ord])]
  atoms <- atoms[first, , drop=FALSE]
  atom_keys_full <- atom_keys_full[first]
  get_atom <- function(residue, atom) {
    keys <- paste(residue_keys, atom, sep="\r")
    idx <- match(keys, atom_keys_full)
    coordinates <- matrix(NA_real_, nrow=n, ncol=3L)
    ok <- which(!is.na(idx))
    if (length(ok))
      coordinates[ok, ] <- as.matrix(atoms[idx[ok],c("x","y","z"),drop=FALSE])
    coordinates
  }
  ca <- get_atom(seq_len(n), "CA")
  c_atom <- get_atom(seq_len(n), "C")
  n_atom <- get_atom(seq_len(n), "N")
  cb <- get_atom(seq_len(n), "CB")
  cb_rows <- which(complete.cases(cbind(ca,cb)))
  if(length(cb_rows)) result$cb_ca_distance[cb_rows] <-
    sqrt(rowSums((cb[cb_rows,,drop=FALSE]-ca[cb_rows,,drop=FALSE])^2))
  # Signed tetrahedral volume describes N–CA–C–CB chirality. It is a
  # measured geometry value, NOT a validated Cβ-deviation/outlier score.
  cb_volume_rows <- which(complete.cases(cbind(n_atom,ca,c_atom,cb)))
  if(length(cb_volume_rows)) {
    v1 <- n_atom[cb_volume_rows,,drop=FALSE]-ca[cb_volume_rows,,drop=FALSE]
    v2 <- c_atom[cb_volume_rows,,drop=FALSE]-ca[cb_volume_rows,,drop=FALSE]
    v3 <- cb[cb_volume_rows,,drop=FALSE]-ca[cb_volume_rows,,drop=FALSE]
    cross <- cbind(v1[,2]*v2[,3]-v1[,3]*v2[,2],
                   v1[,3]*v2[,1]-v1[,1]*v2[,3],
                   v1[,1]*v2[,2]-v1[,2]*v2[,1])
    result$cb_signed_volume[cb_volume_rows] <- rowSums(cross*v3)
  }
  chi_target <- unname(ram_geom_chi1_atom[as.character(torsions$resn)])
  target <- matrix(NA_real_, nrow=n, ncol=3L)
  for (atom_name in unique(chi_target[!is.na(chi_target)])) {
    rows <- which(chi_target == atom_name & !is.na(chi_target))
    target[rows, ] <- get_atom(seq_len(n),atom_name)[rows,,drop=FALSE]
  }
  chi_rows <- which(complete.cases(cbind(n_atom,ca,cb,target)))
  if (length(chi_rows)) {
    result$chi1[chi_rows] <- ram_dihedral_batch(
      n_atom[chi_rows,,drop=FALSE],ca[chi_rows,,drop=FALSE],
      cb[chi_rows,,drop=FALSE],target[chi_rows,,drop=FALSE])
    result$chi1_available[chi_rows] <- is.finite(result$chi1[chi_rows])
  }
  if (n>1L && "bonded_to_next" %in% names(torsions)) {
    rows <- which(torsions$bonded_to_next[seq_len(n-1L)] %in% TRUE)
    if (length(rows)) {
      j <- rows+1L
      distance <- sqrt(rowSums((c_atom[rows,,drop=FALSE]-
                                n_atom[j,,drop=FALSE])^2))
      valid_pair <- is.finite(distance) & distance>=peptide_min &
                    distance<=peptide_max &
                    complete.cases(cbind(ca[rows,,drop=FALSE],
                      c_atom[rows,,drop=FALSE], n_atom[j,,drop=FALSE],
                      ca[j,,drop=FALSE]))
      rows <- rows[valid_pair]; j <- j[valid_pair]
      if(length(rows)) {
        result$peptide_bond_length[rows] <-
          sqrt(rowSums((c_atom[rows,,drop=FALSE]-
                        n_atom[j,,drop=FALSE])^2))
        result$omega[rows] <- ram_dihedral_batch(
          ca[rows,,drop=FALSE],c_atom[rows,,drop=FALSE],
          n_atom[j,,drop=FALSE],ca[j,,drop=FALSE])
      }
    }
  }
  result$omega_status <- ram_geom_omega_status(result$omega)
  result
}

ram_join_geometry <- function(classified, geometry) {
  if (nrow(classified)!=nrow(geometry))
    stop("Extra-geometry rows must match the extracted backbone rows.")
  cbind(classified,geometry)
}
