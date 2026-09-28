# Atom-based peptide bookkeeping. Every chain and insertion code is a
# separate residue; torsions cross boundaries only when a C-N bond is present.
ram_cross <- function(a, b) {
  c(a[2L] * b[3L] - a[3L] * b[2L],
    a[3L] * b[1L] - a[1L] * b[3L],
    a[1L] * b[2L] - a[2L] * b[1L])
}

ram_dihedral <- function(p0, p1, p2, p3) {
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
  connected <- function(i, j) {
    if (i < 1L || j > n || records[[i]]$chain != records[[j]]$chain) {
      return(FALSE)
    }
    c_atom <- records[[i]]$C
    n_atom <- records[[j]]$N
    if (anyNA(c_atom) || anyNA(n_atom)) return(FALSE)
    distance <- sqrt(sum((c_atom - n_atom)^2))
    is.finite(distance) && distance >= min_peptide_bond &&
      distance <= max_peptide_bond
  }
  # Allocate output once. Appending rows to a data frame for each residue
  # repeatedly copies the growing table on large structures.
  phi <- rep(NA_real_, n)
  psi <- rep(NA_real_, n)
  bonded_to_next <- rep(FALSE, n)
  next_resn <- rep(NA_character_, n)
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      bonded_to_next[[i]] <- connected(i, i + 1L)
      if (bonded_to_next[[i]]) next_resn[[i]] <- records[[i + 1L]]$resn
    }
  }
  for (i in seq_len(n)) {
    rec <- records[[i]]
    if (anyNA(c(rec$N, rec$CA, rec$C))) next
    if (i > 1L && bonded_to_next[[i - 1L]]) {
      phi[[i]] <- ram_dihedral(records[[i - 1L]]$C, rec$N, rec$CA, rec$C)
    }
    if (bonded_to_next[[i]]) {
      psi[[i]] <- ram_dihedral(rec$N, rec$CA, rec$C, records[[i + 1L]]$N)
    }
  }
  data.frame(
    resi = vapply(records, function(rec) as.integer(rec$resi), integer(1)),
    insertion_code = vapply(records, function(rec) as.character(rec$insertion_code),
                            character(1)),
    chain = vapply(records, function(rec) as.character(rec$chain), character(1)),
    resn = vapply(records, function(rec) as.character(rec$resn), character(1)),
    phi = phi, psi = psi, next_resn = next_resn,
    bonded_to_next = bonded_to_next,
    stringsAsFactors = FALSE
  )
}
