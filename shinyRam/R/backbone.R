# Atom-based peptide bookkeeping. Every chain and insertion code is a
# separate residue; torsions cross boundaries only when a C-N bond is present.
ram_cross <- function(a, b) {
  c(a[2L] * b[3L] - a[3L] * b[2L],
    a[3L] * b[1L] - a[1L] * b[3L],
    a[1L] * b[2L] - a[2L] * b[1L])
}

ram_dihedral <- function(p0, p1, p2, p3) {
  if (any(!is.finite(c(p0, p1, p2, p3)))) return(NA_real_)
  b0 <- p0 - p1
  b1 <- p2 - p1
  b2 <- p3 - p2
  norm <- sqrt(sum(b1 * b1))
  if (!is.finite(norm) || norm < 1e-10) return(NA_real_)
  b1 <- b1 / norm
  v <- b0 - sum(b0 * b1) * b1
  w <- b2 - sum(b2 * b1) * b1
  if (sqrt(sum(v * v)) < 1e-10 || sqrt(sum(w * w)) < 1e-10) {
    return(NA_real_)
  }
  atan2(sum(ram_cross(b1, v) * w), sum(v * w)) * 180 / pi
}


# Compute many dihedral angles at once. Rows of each argument are XYZ points.
ram_dihedral_batch <- function(p0, p1, p2, p3) {
  n <- nrow(p0)
  stopifnot(is.matrix(p0), all(dim(p0) == c(n, 3L)),
            identical(dim(p0), dim(p1)), identical(dim(p0), dim(p2)),
            identical(dim(p0), dim(p3)))
  if (n == 0L) return(numeric())
  b0 <- p0 - p1
  b1 <- p2 - p1
  b2 <- p3 - p2
  len <- sqrt(rowSums(b1^2))
  u <- b1 / pmax(len, 1e-10)
  v <- b0 - u * rowSums(b0 * u)
  w <- b2 - u * rowSums(b2 * u)
  good <- complete.cases(cbind(p0, p1, p2, p3)) &
    len > 1e-10 & sqrt(rowSums(v^2)) > 1e-10 &
    sqrt(rowSums(w^2)) > 1e-10
  out <- rep(NA_real_, n)
  if (!any(good)) return(out)
  uv <- cbind(
    u[, 2L] * v[, 3L] - u[, 3L] * v[, 2L],
    u[, 3L] * v[, 1L] - u[, 1L] * v[, 3L],
    u[, 1L] * v[, 2L] - u[, 2L] * v[, 1L]
  )
  out[good] <- atan2(rowSums((uv * w)[good, , drop = FALSE]),
                     rowSums((v * w)[good, , drop = FALSE])) * 180 / pi
  out
}

ram_extract_torsions <- function(pdb, amino_acids = c(
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
  "TRP", "TYR", "VAL"
), max_peptide_bond = 1.9, min_peptide_bond = 1.0) {
  atoms <- pdb$atom
  required <- c("chain", "resno", "resid", "elety", "x", "y", "z")
  if (!is.data.frame(atoms) || !all(required %in% names(atoms))) {
    stop("Expected PDB/mmCIF atom records with residue and backbone coordinates.")
  }
  if ("model" %in% names(atoms) && length(unique(atoms$model)) > 1L) {
    stop("Multiple structural models found; choose a single model before analysis.")
  }
  empty <- data.frame(
    resi = integer(), insertion_code = character(), chain = character(),
    resn = character(), phi = numeric(), psi = numeric(),
    next_resn = character(), bonded_to_next = logical(),
    stringsAsFactors = FALSE
  )
  if (nrow(atoms) == 0L) return(empty)
  atoms$chain[is.na(atoms$chain)] <- ""
  atoms$resid <- trimws(as.character(atoms$resid))
  atoms$elety <- trimws(as.character(atoms$elety))
  atoms$insert <- if ("insert" %in% names(atoms)) {
    trimws(as.character(atoms$insert))
  } else rep("", nrow(atoms))
  atoms$insert[is.na(atoms$insert)] <- ""
  atoms$alt <- if ("alt" %in% names(atoms)) {
    trimws(as.character(atoms$alt))
  } else rep("", nrow(atoms))
  atoms$alt[is.na(atoms$alt)] <- ""
  atoms <- atoms[!is.na(atoms$resid) & atoms$resid %in% amino_acids, , drop = FALSE]
  if (!nrow(atoms)) return(empty)

  # The ordered group index prevents accidental joins between insertion codes
  # or chains that reuse the same residue numbering.
  key <- paste(atoms$chain, atoms$resno, atoms$insert, sep = "\r")
  indices <- split(seq_len(nrow(atoms)), factor(key, levels = unique(key)))
  pick_atom <- function(record, atom_name) {
    options <- record[record$elety == atom_name &
                      record$alt %in% c("", "A") &
                      is.finite(record$x) & is.finite(record$y) &
                      is.finite(record$z), , drop = FALSE]
    if (!nrow(options)) return(rep(NA_real_, 3L))
    # Prefer the unlabelled conformation, followed by alternate A.
    chosen <- if (any(options$alt == "")) {
      options[which(options$alt == "")[1L], , drop = FALSE]
    } else options[1L, , drop = FALSE]
    as.numeric(unlist(chosen[1L, c("x", "y", "z")], use.names = FALSE))
  }
  records <- lapply(indices, function(idx) {
    rec <- atoms[idx, , drop = FALSE]
    list(
      resi = rec$resno[1L], insertion_code = rec$insert[1L],
      chain = rec$chain[1L], resn = rec$resid[1L],
      N = pick_atom(rec, "N"), CA = pick_atom(rec, "CA"),
      C = pick_atom(rec, "C")
    )
  })
  n <- length(records)
  # Matrix operations avoid thousands of tiny R calls for large assemblies.
  coords <- function(atom_name) {
    do.call(rbind, lapply(records, function(rec) rec[[atom_name]]))
  }
  nxyz <- coords("N")
  caxyz <- coords("CA")
  cxyz <- coords("C")
  chains <- vapply(records, function(rec) as.character(rec$chain), character(1))
  residue_names <- vapply(records, function(rec) as.character(rec$resn),
                          character(1))

  bonded_to_next <- rep(FALSE, n)
  next_resn <- rep(NA_character_, n)
  if (n > 1L) {
    i <- seq_len(n - 1L)
    j <- i + 1L
    offset <- cxyz[i, , drop = FALSE] - nxyz[j, , drop = FALSE]
    dist <- sqrt(rowSums(offset^2))
    bonded_to_next[i] <- chains[i] == chains[j] & is.finite(dist) &
      dist >= min_peptide_bond & dist <= max_peptide_bond
    next_resn[i[bonded_to_next[i]]] <- residue_names[j[bonded_to_next[i]]]
  }
  valid_backbone <- complete.cases(cbind(nxyz, caxyz, cxyz))
  phi <- rep(NA_real_, n)
  psi <- rep(NA_real_, n)
  if (n > 1L) {
    i_phi <- which(bonded_to_next[seq_len(n - 1L)] &
                     valid_backbone[2L:n]) + 1L
    if (length(i_phi)) {
      phi[i_phi] <- ram_dihedral_batch(
        cxyz[i_phi - 1L, , drop = FALSE],
        nxyz[i_phi, , drop = FALSE],
        caxyz[i_phi, , drop = FALSE],
        cxyz[i_phi, , drop = FALSE])
    }
    i_psi <- which(bonded_to_next[seq_len(n - 1L)] &
                     valid_backbone[seq_len(n - 1L)])
    if (length(i_psi)) {
      psi[i_psi] <- ram_dihedral_batch(
        nxyz[i_psi, , drop = FALSE],
        caxyz[i_psi, , drop = FALSE],
        cxyz[i_psi, , drop = FALSE],
        nxyz[i_psi + 1L, , drop = FALSE])
    }
  }
  data.frame(
    resi = vapply(records, function(rec) as.integer(rec$resi), integer(1)),
    insertion_code = vapply(records, function(rec) as.character(rec$insertion_code),
                            character(1)),
    chain = chains, resn = residue_names,
    phi = phi, psi = psi, next_resn = next_resn,
    bonded_to_next = bonded_to_next, stringsAsFactors = FALSE
  )
}
