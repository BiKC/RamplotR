# Import external wwPDB/MolProbity-backed diagnostics without changing
# RamplotR's original reference grids or assigning our own rotamer labels.
# xml2 is optional and used only when a user supplies an official report.

ram_external_validation_read <- function(path, max_bytes=32000000L,
                                         max_residues=100000L) {
  if (!requireNamespace("xml2",quietly=TRUE))
    stop("Install xml2 to read official wwPDB validation XML.")
  if (!file.exists(path) || !is.finite(file.info(path)$size) ||
      file.info(path)$size <= 0 || file.info(path)$size > max_bytes)
    stop("Provide a nonempty wwPDB validation XML (or XML.gz), 32 MB maximum.")
  raw <- readBin(path,what="raw",n=file.info(path)$size)
  if (grepl("\\.gz$",path,ignore.case=TRUE))
    raw <- memDecompress(raw,type="gzip")
  if (length(raw) > max_bytes)
    stop("The uncompressed wwPDB report exceeds the 32 MB limit.")
  doc <- xml2::read_xml(raw)
  nodes <- xml2::xml_find_all(doc,".//*[local-name()='ModelledSubgroup']")
  if (!length(nodes) || length(nodes) > max_residues)
    stop("No supported residue records found, or report exceeds the size limit.")
  attr <- function(name) {
    z <- trimws(xml2::xml_attr(nodes,name))
    z[z %in% c("",".","?")] <- NA_character_
    z
  }
  count <- function(tag) vapply(nodes,function(node)
    length(xml2::xml_find_all(node,
      sprintf("./*[local-name()='%s']",tag))),integer(1))
  rota <- tolower(attr("rota"))
  rama <- tolower(attr("rama"))
  rama[rama=="favoured"] <- "favored"
  # These labels originate from the official wwPDB report, not RamplotR.
  records <- data.frame(
    model=attr("model"),chain=attr("chain"),
    resi=suppressWarnings(as.integer(attr("resnum"))),
    insertion_code=attr("icode"),resn=toupper(attr("resname")),
    altcode=attr("altcode"),wwpdb_rama=rama,wwpdb_rotamer=rota,
    wwpdb_clashes=count("clash"),wwpdb_symmetry_clashes=count("symm-clash"),
    wwpdb_bond_outliers=count("bond-outlier"),
    wwpdb_angle_outliers=count("angle-outlier"),
    wwpdb_rscc=suppressWarnings(as.numeric(attr("rscc"))),
    wwpdb_rsrz=suppressWarnings(as.numeric(attr("rsrz"))),
    stringsAsFactors=FALSE
  )
  records$chain[is.na(records$chain)] <- ""
  records$insertion_code[is.na(records$insertion_code)] <- ""
  records$altcode[is.na(records$altcode)] <- ""
  records
}

ram_external_validation_join <- function(residues, official, model=1L) {
  required <- c("chain","resi","insertion_code","resn")
  if (!all(required %in% names(residues)) ||
      !all(c(required,"model","altcode") %in% names(official)))
    stop("Residue identifiers and a wwPDB model number are required.")
  other <- official[!is.na(official$model) &
    official$model==as.character(model) & !is.na(official$resi),,drop=FALSE]
  # wwPDB may list multiple alternate conformers. Prefer blank over A,
  # and reject same-priority ambiguity rather than fabricating a match.
  if(nrow(other)) {
    pref <- ifelse(other$altcode %in% c("",".","?"),0L,
                   ifelse(other$altcode=="A",1L,2L))
    other <- other[order(pref),,drop=FALSE]
  }
  key <- function(x) paste(x$chain,x$resi,x$insertion_code,
                           toupper(x$resn),sep="\r")
  ids <- key(other)
  other <- other[!duplicated(ids),,drop=FALSE]
  ids <- unique(ids)
  rows <- match(key(residues),ids)
  out <- residues
  fields <- grep("^wwpdb_",names(official),value=TRUE)
  for(field in fields) out[[field]] <- other[[field]][rows]
  out$wwpdb_matched <- !is.na(rows)
  out$wwpdb_report_model <- rep(as.integer(model),nrow(out))
  out
}

ram_external_validation_summary <- function(data) {
  if (!"wwpdb_matched" %in% names(data)) return(NULL)
  matched <- data$wwpdb_matched %in% TRUE
  counts <- function(column, predicate=function(x) x>0) {
    if (!column %in% names(data)) return(NA_integer_)
    sum(matched & !is.na(data[[column]]) & predicate(data[[column]]))
  }
  list(matched=sum(matched),total=nrow(data),
       official_rama_outliers=counts("wwpdb_rama",function(x) x=="outlier"),
       official_rotamer_outliers=counts("wwpdb_rotamer",
                                     function(x) x %in% c("outlier","outliers")),
       residues_with_clashes=counts("wwpdb_clashes"),
       residues_with_bond_outliers=counts("wwpdb_bond_outliers"),
       residues_with_angle_outliers=counts("wwpdb_angle_outliers"))
}
