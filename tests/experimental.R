source("shinyRam/R/experimental.R")
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
if(!requireNamespace("xml2",quietly=TRUE))
  stop("Install xml2 to run independent-validation import tests.")
source_file <- "validation/fixtures/synthetic-geometry.xml"
official <- ram_external_validation_read(source_file)
assert(nrow(official)==4L, "Official XML residue extraction lost entries.")
assert(official$wwpdb_clashes[[1]]==2L &&
       official$wwpdb_bond_outliers[[1]]==1L &&
       official$wwpdb_angle_outliers[[1]]==1L &&
       official$wwpdb_rotamer[[1]]=="outlier",
       "Independent wwPDB subelements and rotamer labels must be preserved.")
ours <- data.frame(chain=c("A","A","B"),resi=c(1L,2L,1L),
                   insertion_code=c("","A",""),resn=c("SER","THR","ALA"),
                   phi=c(-70,-110,100),stringsAsFactors=FALSE)
joined <- ram_external_validation_join(ours,official,model=1L)
assert(identical(joined$wwpdb_matched,c(TRUE,TRUE,FALSE)),
       "Independent report must match chain, insertion code and resname.")
assert(identical(joined$wwpdb_rama[1:2],c("favored","allowed")) &&
       joined$wwpdb_rotamer[[2]]=="outlier" &&
       is.na(joined$wwpdb_rama[[3]]),
       "Unlabelled conformers must outrank alt A and missing IDs stay missing.")
assert(joined$wwpdb_symmetry_clashes[[2]]==1L,
       "Symmetry clashes must be distinguished from local contacts.")
m2 <- ram_external_validation_join(ours,official,model=2L)
assert(m2$wwpdb_rama[[1]]=="outlier" &&
       !any(m2$wwpdb_matched[2:3]),
       "Validation of one model must never silently mix NMR models.")
summary <- ram_external_validation_summary(joined)
assert(summary$matched==2L && summary$official_rotamer_outliers==2L &&
       summary$residues_with_clashes==1L,
       "Independent report counts must use only matched residues.")
bad <- try(ram_external_validation_read(source_file,max_bytes=8L),silent=TRUE)
assert(inherits(bad,"try-error"),"Oversized XML must be rejected.")
gzpath <- tempfile(fileext=".xml.gz")
con <- gzfile(gzpath,open="wb")
writeBin(readBin(source_file,what="raw",n=file.info(source_file)$size),con)
close(con)
on.exit_gz <- function() unlink(gzpath)
compressed <- ram_external_validation_read(gzpath)
assert(nrow(compressed)==nrow(official),
       "The streamed gzip parser must preserve every official record.")
# A tiny compressed file that expands beyond the configured limit must
# fail while streaming, not after allocating its entire inflated contents.
bombpath <- tempfile(fileext=".xml.gz")
bombcon <- gzfile(bombpath,open="wb")
writeBin(charToRaw(paste0("<root>",strrep(" ",50000),"</root>")),bombcon)
close(bombcon)
assert(file.info(bombpath)$size < 1024L,
       "Synthetic gzip fixture must be smaller than the upload limit.")
blocked <- try(ram_external_validation_read(bombpath,max_bytes=1024L),
               silent=TRUE)
assert(inherits(blocked,"try-error"),
       "Expanded XML limits must also apply to compressed reports.")
unlink(bombpath)
on.exit_gz()
message("External wwPDB geometry/rotamer/clash XML import tests passed.")
