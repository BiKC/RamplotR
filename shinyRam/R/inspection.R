# Pure residue inspection and pairwise comparison helpers.
# Nothing in this file changes torsion extraction or scientific classification.

ram_heavy_atoms <- function(atoms) {
  if (!is.data.frame(atoms) || !nrow(atoms)) return(atoms)
  element <- if ("elesy" %in% names(atoms))
    toupper(trimws(as.character(atoms$elesy))) else rep("",nrow(atoms))
  atom_name <- if ("elety" %in% names(atoms))
    toupper(trimws(as.character(atoms$elety))) else rep("",nrow(atoms))
  hydrogen <- element %in% c("H","D") |
    (!nzchar(element) & grepl("^[0-9]*[HD]",atom_name))
  atoms[!hydrogen,,drop=FALSE]
}

ram_nearby_hetero_context <- function(pdb, row, max_distance = 6,
                                      max_hits = 3L) {
  empty <- data.frame(
    resn=character(),chain=character(),resi=integer(),
    insertion_code=character(),distance=numeric(),
    target_atom=character(),hetero_atom=character(),
    stringsAsFactors=FALSE
  )
  if (!is.list(pdb) || !is.data.frame(pdb$atom) ||
      !is.data.frame(row) || nrow(row)!=1L ||
      is.null(pdb$hetero_atom) || !is.data.frame(pdb$hetero_atom) ||
      !nrow(pdb$hetero_atom)) return(empty)
  required <- c("chain","resno","resid","elety","x","y","z")
  if (!all(required %in% names(pdb$atom)) ||
      !all(required %in% names(pdb$hetero_atom))) return(empty)

  normalize_atoms <- function(atoms) {
    if (!"insert" %in% names(atoms)) atoms$insert <- ""
    atoms$chain <- as.character(atoms$chain)
    atoms$chain[is.na(atoms$chain)] <- ""
    atoms$insert <- as.character(atoms$insert)
    atoms$insert[is.na(atoms$insert)] <- ""
    atoms$resid <- toupper(trimws(as.character(atoms$resid)))
    atoms$elety <- trimws(as.character(atoms$elety))
    atoms
  }
  atoms <- normalize_atoms(pdb$atom)
  hetero <- normalize_atoms(pdb$hetero_atom)
  chain <- as.character(row$chain[[1L]])
  if (is.na(chain)) chain <- ""
  insertion <- as.character(row$insertion_code[[1L]])
  if (is.na(insertion)) insertion <- ""
  resi <- suppressWarnings(as.integer(row$resi[[1L]]))
  resn <- toupper(as.character(row$resn[[1L]]))
  target <- atoms[
    atoms$chain==chain & suppressWarnings(as.integer(atoms$resno))==resi &
      atoms$insert==insertion & atoms$resid==resn,
    ,drop=FALSE
  ]
  target <- ram_heavy_atoms(target)
  hetero <- ram_heavy_atoms(hetero)
  water_names <- c("HOH","WAT","DOD","H2O","SOL","TIP","TIP3","TIP3P")
  hetero <- hetero[!hetero$resid %in% water_names,,drop=FALSE]
  target <- target[is.finite(target$x)&is.finite(target$y)&is.finite(target$z),
                   ,drop=FALSE]
  hetero <- hetero[is.finite(hetero$x)&is.finite(hetero$y)&is.finite(hetero$z),
                   ,drop=FALSE]
  if (!nrow(target) || !nrow(hetero)) return(empty)

  target_xyz <- as.matrix(target[,c("x","y","z"),drop=FALSE])
  nearest <- vapply(seq_len(nrow(hetero)),function(i) {
    delta <- sweep(target_xyz,2L,
      as.numeric(hetero[i,c("x","y","z")]),FUN="-")
    min(sqrt(rowSums(delta^2)))
  },numeric(1))
  hetero$.ram_distance <- nearest
  hetero <- hetero[is.finite(hetero$.ram_distance) &
                   hetero$.ram_distance<=max_distance,,drop=FALSE]
  if (!nrow(hetero)) return(empty)

  keys <- paste(hetero$chain,hetero$resno,hetero$insert,hetero$resid,sep="\r")
  groups <- split(seq_len(nrow(hetero)),keys)
  results <- lapply(groups,function(ix) {
    local <- ix[which.min(hetero$.ram_distance[ix])]
    atom_delta <- sweep(target_xyz,2L,
      as.numeric(hetero[local,c("x","y","z")]),FUN="-")
    target_index <- which.min(sqrt(rowSums(atom_delta^2)))
    data.frame(
      resn=hetero$resid[[local]],
      chain=hetero$chain[[local]],
      resi=suppressWarnings(as.integer(hetero$resno[[local]])),
      insertion_code=hetero$insert[[local]],
      distance=hetero$.ram_distance[[local]],
      target_atom=target$elety[[target_index]],
      hetero_atom=hetero$elety[[local]],
      stringsAsFactors=FALSE
    )
  })
  result <- do.call(rbind,results)
  result <- result[order(result$distance,result$resn,result$chain,result$resi),
                   ,drop=FALSE]
  head(result,max(1L,as.integer(max_hits)))
}

ram_residue_evidence <- function(row, boundary_margin = 2, local_context = NULL,
                                 ensemble_context = NULL) {
  if (!is.data.frame(row) || nrow(row) != 1L)
    stop("Residue evidence expects exactly one residue row.")
  evidence <- list()
  add <- function(level, title, detail, source) {
    evidence[[length(evidence)+1L]] <<- data.frame(
      level=level,title=title,detail=detail,source=source,
      stringsAsFactors=FALSE
    )
  }
  value <- function(name, default=NA) {
    if (!name %in% names(row) || !length(row[[name]])) return(default)
    row[[name]][[1L]]
  }
  phi <- suppressWarnings(as.numeric(value("phi",NA_real_)))
  psi <- suppressWarnings(as.numeric(value("psi",NA_real_)))
  native <- as.character(value("region",NA_character_))
  standard <- as.character(value("rama8000_region",NA_character_))
  plddt <- suppressWarnings(as.numeric(value("plddt",NA_real_)))

  if (!is.finite(phi) || !is.finite(psi)) {
    add("warning","Backbone angles unavailable",
        "Phi/psi cannot be evaluated for this residue, commonly because it is terminal or required backbone atoms are missing.",
        "Coordinates")
  } else {
    if (!is.na(standard) && standard == "Outlier") {
      add("high","Rama8000 backbone outlier",
          "The current six-class standard validation places this residue outside the allowed Rama8000 region.",
          "Rama8000")
    } else if (!is.na(standard) && standard == "Allowed") {
      add("info","Rama8000 allowed region",
          "The residue is outside the favored Rama8000 region but remains within the allowed region.",
          "Rama8000")
    }
    if (!is.na(native) && native == "Not allowed") {
      add("warning","Unusual RamplotR density position",
          "The selected RamplotR reference distribution labels this position Not allowed. This is separate from Rama8000 outlier status.",
          "RamplotR density")
    }
    density <- suppressWarnings(as.numeric(value("density",NA_real_)))
    if (is.finite(density) &&
        any(abs(density-c(85,98,99.95)) <= boundary_margin)) {
      add("info","Near a RamplotR contour boundary",
          sprintf("Density percentile %.1f lies within %.1f points of a display-classification contour.",
                  density,boundary_margin),
          "RamplotR density")
    }
  }

  if (is.finite(plddt)) {
    if (plddt < 50) {
      add("warning","Very low prediction confidence",
          sprintf("pLDDT %.1f indicates very low local model confidence; local geometry should be interpreted cautiously.",plddt),
          "Prediction confidence")
    } else if (plddt < 70) {
      add("info","Low prediction confidence",
          sprintf("pLDDT %.1f indicates low local model confidence.",plddt),
          "Prediction confidence")
    }
    if (plddt >= 90 && !is.na(standard) && standard == "Outlier") {
      add("high","High-confidence prediction with unusual backbone geometry",
          sprintf("pLDDT %.1f is very high while Rama8000 classifies the backbone as an outlier; inspect the local structural context.",plddt),
          "Combined evidence")
    }
  }

  omega_status <- as.character(value("omega_status",NA_character_))
  omega <- suppressWarnings(as.numeric(value("omega",NA_real_)))
  if (!is.na(omega_status) && omega_status == "Twisted") {
    add("high","Twisted peptide bond",
        if (is.finite(omega)) sprintf("Peptide omega is %.1f degrees.",omega)
        else "The preceding peptide bond is classified as twisted.",
        "Local geometry")
  } else if (!is.na(omega_status) && omega_status == "Cis") {
    add("info","Cis peptide bond",
        if (is.finite(omega)) sprintf("Peptide omega is %.1f degrees.",omega)
        else "The preceding peptide bond is cis.",
        "Local geometry")
  }

  official_rama <- tolower(as.character(value("wwpdb_rama",NA_character_)))
  if (!is.na(official_rama) && official_rama == "outlier")
    add("high","Official wwPDB Ramachandran outlier",
        "The attached official validation report identifies this residue as a Ramachandran outlier.",
        "wwPDB")
  rotamer <- tolower(as.character(value("wwpdb_rotamer",NA_character_)))
  if (!is.na(rotamer) && rotamer %in% c("outlier","outliers"))
    add("warning","Official wwPDB rotamer outlier",
        "The attached official validation report identifies the side-chain rotamer as an outlier.",
        "wwPDB")
  clashes <- suppressWarnings(as.numeric(value("wwpdb_clashes",NA_real_)))
  if (is.finite(clashes) && clashes > 0)
    add("warning","Official wwPDB local clash",
        sprintf("The attached report records %d local clash%s for this residue.",
                as.integer(clashes),if (clashes==1) "" else "es"),
        "wwPDB")
  bond <- suppressWarnings(as.numeric(value("wwpdb_bond_outliers",NA_real_)))
  angle <- suppressWarnings(as.numeric(value("wwpdb_angle_outliers",NA_real_)))
  if ((is.finite(bond) && bond > 0) || (is.finite(angle) && angle > 0))
    add("warning","Official covalent-geometry outlier",
        sprintf("The attached report records %d bond-length and %d bond-angle outlier%s.",
          ifelse(is.finite(bond),as.integer(bond),0L),
          ifelse(is.finite(angle),as.integer(angle),0L),
          ifelse((ifelse(is.finite(bond),bond,0)+ifelse(is.finite(angle),angle,0))==1,
                 "","s")),
        "wwPDB")

  if (is.data.frame(ensemble_context) && nrow(ensemble_context)==1L) {
    ensemble_value <- function(name, default=NA) {
      if (!name %in% names(ensemble_context) ||
          !length(ensemble_context[[name]])) return(default)
      ensemble_context[[name]][[1L]]
    }
    phi_sd <- suppressWarnings(as.numeric(ensemble_value("phi_sd",NA_real_)))
    psi_sd <- suppressWarnings(as.numeric(ensemble_value("psi_sd",NA_real_)))
    finite_spread <- c(phi_sd,psi_sd)
    finite_spread <- finite_spread[is.finite(finite_spread)]
    spread <- if(length(finite_spread)) max(finite_spread) else NA_real_
    plddt_mean <- suppressWarnings(as.numeric(
      ensemble_value("plddt_mean",NA_real_)))
    plddt_sd <- suppressWarnings(as.numeric(
      ensemble_value("plddt_sd",NA_real_)))
    models_present <- suppressWarnings(as.integer(
      ensemble_value("models_present",NA_integer_)))
    total_models <- suppressWarnings(as.integer(
      ensemble_value("ensemble_models_total",NA_integer_)))
    standard_changes <- isTRUE(ensemble_value("rama8000_changes",FALSE))
    standard_consistency <- suppressWarnings(as.numeric(
      ensemble_value("rama8000_consistency",NA_real_)))
    basin_changes <- isTRUE(ensemble_value("basin_changes",FALSE))
    basin_consistency <- suppressWarnings(as.numeric(
      ensemble_value("basin_consistency",NA_real_)))
    basin_mode <- as.character(ensemble_value("basin_mode",NA_character_))

    if (is.finite(spread) && spread >= 20 && is.finite(plddt_mean) &&
        plddt_mean >= 90) {
      add("high","High-confidence predictions disagree on local backbone",
          sprintf(paste0(
            "Across the prediction ensemble, the largest circular backbone ",
            "SD is %.1f degrees while mean pLDDT is %.1f. The models are ",
            "individually confident but do not converge on one local ",
            "backbone conformation."),spread,plddt_mean),
          "Prediction ensemble")
    } else if (is.finite(spread) && spread >= 30) {
      add("warning","Strong prediction-ensemble backbone disagreement",
          sprintf(paste0(
            "The largest circular SD across phi/psi is %.1f degrees. This ",
            "describes disagreement between prediction models or seeds, not ",
            "experimental molecular motion."),spread),
          "Prediction ensemble")
    } else if (is.finite(spread) && spread >= 15) {
      add("info","Moderate prediction-ensemble backbone variation",
          sprintf(paste0(
            "The largest circular SD across phi/psi is %.1f degrees across ",
            "the analysed prediction models."),spread),
          "Prediction ensemble")
    }

    if (basin_changes) {
      add("warning","Prediction models choose different backbone states",
          if (is.finite(basin_consistency))
            sprintf(paste0(
              "Only %.1f%% of models with finite phi/psi occupy the modal ",
              "coarse backbone state%s. This state label is a comparison aid, ",
              "not a secondary-structure assignment."),
              100*basin_consistency,
              if(!is.na(basin_mode)) paste0(" (",basin_mode,")") else "")
          else paste0(
            "The prediction models occupy different coarse backbone states. ",
            "These labels are comparison aids, not secondary-structure assignments."),
          "Prediction ensemble")
    }

    if (standard_changes) {
      add("warning","Rama8000 category differs across prediction models",
          if (is.finite(standard_consistency))
            sprintf("Only %.1f%% of classified ensemble models share the modal Rama8000 category at this residue.",
                    100*standard_consistency)
          else "The prediction models do not all share the same Rama8000 category at this residue.",
          "Prediction ensemble")
    }

    if (is.finite(plddt_sd) && plddt_sd >= 10) {
      add("info","Prediction confidence varies across ensemble",
          sprintf("pLDDT has an SD of %.1f across the contributing prediction models.",
                  plddt_sd),
          "Prediction ensemble")
    }

    if (is.finite(models_present) && is.finite(total_models) &&
        total_models > 0L && models_present < total_models) {
      add("info","Residue is absent from some prediction models",
          sprintf("%d of %d analysed models contain this exact chain/residue/insertion/amino-acid identity.",
                  models_present,total_models),
          "Prediction ensemble")
    }
  }

  if (is.data.frame(local_context) && nrow(local_context)) {
    label <- function(i) {
      chain <- as.character(local_context$chain[[i]])
      number <- local_context$resi[[i]]
      insertion <- as.character(local_context$insertion_code[[i]])
      position <- if (is.finite(number))
        paste0(if(nzchar(chain)) paste0(chain,":") else "",
               as.integer(number),insertion)
        else if(nzchar(chain)) chain else "unnumbered"
      sprintf("%s %s at %.1f Å",
        local_context$resn[[i]],position,local_context$distance[[i]])
    }
    nearby <- paste(vapply(seq_len(nrow(local_context)),label,character(1L)),
                    collapse="; ")
    add("info",
        if(nrow(local_context)==1L) "Nearby non-water hetero residue"
        else "Nearby non-water hetero residues",
        paste0(nearby,
          ". Distances are nearest heavy-atom distances and indicate spatial proximity only, not biochemical binding."),
        "Local structure context")
  }

  if (!length(evidence))
    return(data.frame(level=character(),title=character(),
      detail=character(),source=character(),stringsAsFactors=FALSE))
  out <- do.call(rbind,evidence)
  priority <- match(out$level,c("high","warning","info"))
  out[order(priority,out$source,out$title),,drop=FALSE]
}

ram_review_queue <- function(data, boundary_margin = 2) {
  if (!nrow(data)) return(data)
  missing <- is.na(data$phi) | is.na(data$psi) | is.na(data$region)
  not_allowed <- !missing & data$region == "Not allowed"
  standard_outlier <- if ("rama8000_region" %in% names(data))
    !missing & !is.na(data$rama8000_region) & data$rama8000_region == "Outlier"
    else rep(FALSE, nrow(data))
  # The density percentile is cumulative mass above the residue's density.
  # Closeness to any RamplotR contour probability is a review hint only.
  near <- !missing & !is.na(data$density) &
    vapply(data$density, function(value) {
      any(abs(value - c(85, 98, 99.95)) <= boundary_margin)
    }, logical(1))
  data$review_status <- ifelse(missing, "Missing angles",
    ifelse(standard_outlier, "Rama8000 outlier",
      ifelse(not_allowed, "Not allowed",
        ifelse(near, "Near boundary", "Other"))))
  priority <- match(data$review_status,
                    c("Rama8000 outlier", "Not allowed", "Missing angles",
                      "Near boundary", "Other"))
  data[order(priority, data$chain, data$resi, data$insertion_code), ,
       drop = FALSE]
}

ram_amino_acid_letters <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E", GLN="Q",
  GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
  PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V"
)

ram_chain_query_sequence <- function(data, chain) {
  residues <- ram_sequence_data(data)
  if (!nrow(residues)) return(list(sequence="", residues=residues,
                                   known_fraction=NA_real_))
  target <- as.character(chain)
  rows <- residues$chain == target
  residues <- residues[rows, , drop=FALSE]
  if (!nrow(residues)) return(list(sequence="", residues=residues,
                                   known_fraction=NA_real_))
  letters <- toupper(as.character(residues$letter))
  letters[!grepl("^[ACDEFGHIKLMNPQRSTVWY]$", letters)] <- "X"
  known <- letters != "X"
  list(
    sequence=paste0(letters,collapse=""),
    residues=residues,
    known_fraction=mean(known)
  )
}

ram_sequence_data <- function(data) {
  if (!nrow(data)) return(data.frame(
    chain=character(), resi=integer(), insertion_code=character(),
    resn=character(), letter=character(), region=character(),
    rama8000_region=character(), rama8000_group=character(),
    rama8000_score=numeric(), plddt=numeric(),
    stringsAsFactors=FALSE))
  key <- paste(data$chain, data$resi, data$insertion_code, sep="\r")
  data <- data[!duplicated(key), , drop=FALSE]
  out <- data.frame(
    chain=as.character(data$chain),
    resi=as.integer(data$resi),
    insertion_code=as.character(data$insertion_code),
    resn=as.character(data$resn),
    letter=unname(ram_amino_acid_letters[data$resn]),
    region=as.character(data$region),
    rama8000_region=if ("rama8000_region" %in% names(data))
      as.character(data$rama8000_region) else rep(NA_character_,nrow(data)),
    rama8000_group=if ("rama8000_group" %in% names(data))
      as.character(data$rama8000_group) else rep(NA_character_,nrow(data)),
    rama8000_score=if ("rama8000_score" %in% names(data))
      as.numeric(data$rama8000_score) else rep(NA_real_,nrow(data)),
    plddt=if ("plddt" %in% names(data)) as.numeric(data$plddt) else
      rep(NA_real_, nrow(data)),
    stringsAsFactors=FALSE
  )
  out$letter[is.na(out$letter)] <- "X"
  out
}

ram_angular_difference <- function(a, b) {
  delta <- ((as.numeric(b) - as.numeric(a) + 180) %% 360) - 180
  delta[is.na(a) | is.na(b)] <- NA_real_
  delta
}

# A simple local displacement in phi/psi space after each angular component
# has been wrapped independently. This is a navigation/ranking measure, not a
# statistical significance score or a Cartesian structural distance.
ram_backbone_angular_displacement <- function(delta_phi, delta_psi) {
  phi <- as.numeric(delta_phi)
  psi <- as.numeric(delta_psi)
  out <- sqrt(phi^2 + psi^2)
  out[!is.finite(phi) | !is.finite(psi)] <- NA_real_
  out
}

ram_backbone_shift_band <- function(displacement) {
  value <- as.numeric(displacement)
  out <- rep("Unavailable", length(value))
  finite <- is.finite(value)
  out[finite & value < 15] <- "Small"
  out[finite & value >= 15 & value < 30] <- "Moderate"
  out[finite & value >= 30 & value < 60] <- "Large"
  out[finite & value >= 60] <- "Very large"
  out
}

# Needleman-Wunsch global alignment of one chain from each structure. Rows
# containing gaps are retained for display; differences are NA without both
# measured angles. A modest cell limit prevents unbounded Shiny allocations.
ram_align_residues <- function(a, b, max_cells = 4e6) {
  n <- nrow(a); m <- nrow(b)
  if (n * m > max_cells)
    stop("Comparison exceeds the sequence-alignment size limit; select shorter chains.")
  if (!n || !m) return(data.frame(
    index_a = if (n) seq_len(n) else rep(NA_integer_, m),
    index_b = if (m) seq_len(m) else rep(NA_integer_, n)
  ))
  symbols_a <- unname(ram_amino_acid_letters[a$resn])
  symbols_b <- unname(ram_amino_acid_letters[b$resn])
  symbols_a[is.na(symbols_a)] <- "X"
  symbols_b[is.na(symbols_b)] <- "X"
  # Dynamic programming retains only an integer direction for traceback.
  score <- matrix(0, n + 1L, m + 1L)
  direction <- matrix(0L, n + 1L, m + 1L)
  score[, 1L] <- -2 * (0:n)
  score[1L, ] <- -2 * (0:m)
  direction[-1L, 1L] <- 2L
  direction[1L, -1L] <- 3L
  for (i in seq_len(n)) for (j in seq_len(m)) {
    choice <- c(score[i, j] + if (symbols_a[i] == symbols_b[j]) 2 else -1,
                score[i, j + 1L] - 2, score[i + 1L, j] - 2)
    best <- which.max(choice)
    score[i + 1L, j + 1L] <- choice[[best]]
    direction[i + 1L, j + 1L] <- best
  }
  i <- n; j <- m
  ia <- integer(n + m); ib <- integer(n + m); k <- 0L
  while (i > 0L || j > 0L) {
    k <- k + 1L
    d <- direction[i + 1L, j + 1L]
    if (d == 1L) {
      ia[k] <- i; ib[k] <- j; i <- i - 1L; j <- j - 1L
    } else if (d == 2L) {
      ia[k] <- i; ib[k] <- NA_integer_; i <- i - 1L
    } else {
      ia[k] <- NA_integer_; ib[k] <- j; j <- j - 1L
    }
  }
  data.frame(index_a=rev(ia[seq_len(k)]),
             index_b=rev(ib[seq_len(k)]))
}

ram_compare_torsions <- function(a, b, pairing=NULL) {
  if(is.null(pairing)) {
    pairing <- ram_align_residues(a,b)
  } else {
    if(!is.data.frame(pairing) ||
       !all(c("index_a","index_b") %in% names(pairing)) ||
       anyNA(pairing[,c("index_a","index_b"),drop=FALSE]) ||
       any(pairing$index_a < 1L | pairing$index_a > nrow(a) |
           pairing$index_b < 1L | pairing$index_b > nrow(b)) ||
       anyDuplicated(pairing$index_a) || anyDuplicated(pairing$index_b))
      stop("Canonical alignment has invalid or ambiguous source indices.",
           call.=FALSE)
  }
  value <- function(data, indices, key, missing) {
    out <- rep(missing, nrow(pairing))
    valid <- !is.na(indices)
    if (any(valid)) out[valid] <- data[[key]][indices[valid]]
    out
  }
  result <- data.frame(
    chain_a=value(a, pairing$index_a, "chain", NA_character_),
    residue_a=value(a, pairing$index_a, "resi", NA_integer_),
    insertion_a=value(a, pairing$index_a, "insertion_code", NA_character_),
    amino_a=value(a, pairing$index_a, "resn", NA_character_),
    chain_b=value(b, pairing$index_b, "chain", NA_character_),
    residue_b=value(b, pairing$index_b, "resi", NA_integer_),
    insertion_b=value(b, pairing$index_b, "insertion_code", NA_character_),
    amino_b=value(b, pairing$index_b, "resn", NA_character_),
    phi_a=value(a, pairing$index_a, "phi", NA_real_),
    phi_b=value(b, pairing$index_b, "phi", NA_real_),
    psi_a=value(a, pairing$index_a, "psi", NA_real_),
    psi_b=value(b, pairing$index_b, "psi", NA_real_),
    region_a=value(a, pairing$index_a, "region", NA_character_),
    region_b=value(b, pairing$index_b, "region", NA_character_),
    stringsAsFactors=FALSE
  )
  result$delta_phi <- ram_angular_difference(result$phi_a, result$phi_b)
  result$delta_psi <- ram_angular_difference(result$psi_a, result$psi_b)
  result$angular_displacement <- ram_backbone_angular_displacement(
    result$delta_phi, result$delta_psi)
  result$shift_band <- ram_backbone_shift_band(result$angular_displacement)
  if (exists("ram_backbone_basin", mode="function")) {
    result$basin_a <- ram_backbone_basin(result$phi_a,result$psi_a)
    result$basin_b <- ram_backbone_basin(result$phi_b,result$psi_b)
    result$basin_changed <- !is.na(result$basin_a) & !is.na(result$basin_b) &
      result$basin_a != result$basin_b
  }
  result$class_changed <- !is.na(result$region_a) &
    !is.na(result$region_b) & result$region_a != result$region_b
  optional_value <- function(data, indices, key, missing) {
    if (!key %in% names(data)) return(rep(missing,nrow(pairing)))
    value(data,indices,key,missing)
  }
  if ("plddt" %in% names(a) || "plddt" %in% names(b)) {
    result$plddt_a <- optional_value(a,pairing$index_a,"plddt",NA_real_)
    result$plddt_b <- optional_value(b,pairing$index_b,"plddt",NA_real_)
    result$confidence_a <- optional_value(
      a,pairing$index_a,"confidence_category",NA_character_)
    result$confidence_b <- optional_value(
      b,pairing$index_b,"confidence_category",NA_character_)
    result$delta_plddt <- result$plddt_b-result$plddt_a
    result$confidence_changed <- !is.na(result$confidence_a) &
      !is.na(result$confidence_b) &
      result$confidence_a != result$confidence_b
  }
  if (all(c("rama8000_region","rama8000_group","rama8000_score") %in% names(a)) &&
      all(c("rama8000_region","rama8000_group","rama8000_score") %in% names(b))) {
    result$rama8000_region_a <- value(
      a, pairing$index_a, "rama8000_region", NA_character_)
    result$rama8000_region_b <- value(
      b, pairing$index_b, "rama8000_region", NA_character_)
    result$rama8000_group_a <- value(
      a, pairing$index_a, "rama8000_group", NA_character_)
    result$rama8000_group_b <- value(
      b, pairing$index_b, "rama8000_group", NA_character_)
    result$rama8000_score_a <- value(
      a, pairing$index_a, "rama8000_score", NA_real_)
    result$rama8000_score_b <- value(
      b, pairing$index_b, "rama8000_score", NA_real_)
    result$rama8000_changed <- !is.na(result$rama8000_region_a) &
      !is.na(result$rama8000_region_b) &
      result$rama8000_region_a != result$rama8000_region_b
  }
  if("uniprot_resi" %in% names(pairing))
    result$uniprot_resi <- as.integer(pairing$uniprot_resi)
  result$alignment <- ifelse(is.na(pairing$index_a), "Insertion",
                      ifelse(is.na(pairing$index_b), "Deletion",
                      ifelse(result$amino_a == result$amino_b,
                             "Match", "Substitution")))
  result
}

ram_comparison_alignment_quality <- function(data) {
  empty <- list(
    aligned=0L,matches=0L,substitutions=0L,
    residues_a=0L,residues_b=0L,
    identity=NA_real_,coverage_a=NA_real_,coverage_b=NA_real_
  )
  if (!is.data.frame(data) || !nrow(data) ||
      !all(c("residue_a","residue_b","alignment") %in% names(data)))
    return(empty)
  has_a <- !is.na(data$residue_a)
  has_b <- !is.na(data$residue_b)
  aligned <- has_a & has_b
  n_aligned <- sum(aligned)
  matches <- sum(aligned & data$alignment=="Match",na.rm=TRUE)
  substitutions <- sum(aligned & data$alignment=="Substitution",na.rm=TRUE)
  residues_a <- sum(has_a)
  residues_b <- sum(has_b)
  list(
    aligned=as.integer(n_aligned),
    matches=as.integer(matches),
    substitutions=as.integer(substitutions),
    residues_a=as.integer(residues_a),
    residues_b=as.integer(residues_b),
    identity=if(n_aligned) matches/n_aligned else NA_real_,
    coverage_a=if(residues_a) n_aligned/residues_a else NA_real_,
    coverage_b=if(residues_b) n_aligned/residues_b else NA_real_
  )
}


# Keep all selected chains available in one navigator, in first-appearance
# order. This is intentionally independent of any "active chain" input.
ram_sequence_groups <- function(data) {
  residues <- ram_sequence_data(data)
  if (!nrow(residues)) return(list())
  chain <- ifelse(is.na(residues$chain) | !nzchar(residues$chain),
                  "Unassigned", residues$chain)
  split(residues, factor(chain, levels = unique(chain)), drop = TRUE)
}

ram_sequence_status <- function(region) {
  out <- rep("missing", length(region))
  out[!is.na(region) & region == "Favoured"] <- "favoured"
  out[!is.na(region) & region == "Allowed"] <- "allowed"
  out[!is.na(region) & region == "Generously allowed"] <- "generously-allowed"
  out[!is.na(region) & region == "Not allowed"] <- "outlier"
  out
}

# Confidence is an independent per-residue signal: never use its colours to
# replace the Ramachandran classification background.
ram_plddt_color <- function(score) {
  vapply(as.numeric(score), function(value) {
    if (!is.finite(value)) "#cbd7db" else if (value < 50) "#d75e56"
    else if (value < 70) "#d6ac52" else if (value < 90) "#7bbcb1"
    else "#126e74"
  }, character(1))
}

# Permanent position labels mark every tenth *PDB residue number*, not every
# tenth item in the sequence. Keep insertion codes on labelled residues.
ram_sequence_position_labels <- function(resi, insertion_code = rep("", length(resi))) {
  stopifnot(length(resi) == length(insertion_code))
  if (!length(resi)) return(character())
  index <- seq_along(resi)
  label <- index == 1L | index == length(resi) |
    (!is.na(resi) & resi %% 10L == 0L)
  out <- rep("", length(resi))
  out[label] <- paste0(resi[label], ifelse(is.na(insertion_code[label]), "",
                                                insertion_code[label]))
  out
}

# Map an NGL/sequence residue to its pair using chain, PDB numbering and
# insertion code. The displayed amino-acid order may contain alignment gaps.
ram_comparison_find <- function(data, side, chain, resi, insertion_code = "") {
  stopifnot(side %in% c("a", "b"))
  if (!nrow(data) || length(chain) != 1L || length(resi) != 1L ||
      length(insertion_code) != 1L || is.na(chain) || is.na(resi) ||
      is.na(insertion_code)) return(NA_integer_)
  number <- suppressWarnings(as.integer(resi))
  if (is.na(number)) return(NA_integer_)
  ix <- which(!is.na(data[[paste0("residue_", side)]]) &
    data[[paste0("chain_", side)]] == as.character(chain) &
    data[[paste0("residue_", side)]] == number &
    data[[paste0("insertion_", side)]] == as.character(insertion_code))
  if (!length(ix)) NA_integer_ else as.integer(ix[[1L]])
}

# Small overview strips represent *positions*, not aggregate percentages.
# For long chains each strip cell represents a consecutive residue bin.
# A bin takes the highest-priority review state so isolated outliers remain
# visible instead of vanishing when hundreds of residues are compressed.
ram_sequence_overview_bins <- function(region, max_bins = 180L) {
  if (!length(region)) return(character(0))
  stopifnot(is.numeric(max_bins), length(max_bins) == 1L,
            is.finite(max_bins), max_bins >= 1L)
  status <- ram_sequence_status(region)
  bins <- min(length(status), as.integer(max_bins))
  bucket <- pmin(bins, floor((seq_along(status) - 1) * bins /
                              length(status)) + 1L)
  priority <- c("outlier", "missing", "generously-allowed",
                "allowed", "favoured")
  vapply(seq_len(bins), function(i) {
    candidates <- status[bucket == i]
    priority[match(TRUE, priority %in% candidates)]
  }, character(1))
}


# Per-chain pLDDT strip uses the lowest known confidence in each consecutive
# sequence bin; low-confidence linkers remain visible on long chains.
ram_plddt_overview_bins <- function(plddt, max_bins = 180L) {
  if (!length(plddt)) return(numeric())
  bins <- min(length(plddt), as.integer(max_bins))
  bucket <- pmin(bins, floor((seq_along(plddt) - 1) * bins /
                              length(plddt)) + 1L)
  vapply(seq_len(bins), function(i) {
    values <- plddt[bucket == i]
    if (!any(is.finite(values))) NA_real_ else min(values[is.finite(values)])
  }, numeric(1))
}
