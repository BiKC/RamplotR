# Pure unit coverage; no Shiny session or network required.
source(file.path("shinyRam","R","inspection.R"))
source(file.path("shinyRam","R","io.R"))
assert <- function(test, message) if (!isTRUE(test)) stop(message)
a <- data.frame(
  chain=c("A","A","A","B"), resi=c(1L,2L,3L,1L),
  insertion_code=c("","A","",""),resn=c("ALA","THR","GLY","PRO"),
  phi=c(-60,NA,179, -78),psi=c(-42,NA, -179,65),
  region=c("Favoured",NA,"Not allowed","Allowed"),
  density=c(42,NA,99,86), stringsAsFactors=FALSE
)
q <- ram_review_queue(a)
assert(identical(q$review_status, c("Outlier","Missing angles","Near boundary","Other")),
       "Review queue priority and missing angles")
assert(identical(q$resn, c("GLY","THR","PRO","ALA")),
       "Review queue ordering")
seq <- ram_sequence_data(a)
assert(identical(seq$letter, c("A","T","G","P")),
       "Sequence must retain amino-acid identities")
assert(seq$insertion_code[2L]=="A", "Insertion code must be preserved")
assert(isTRUE(all.equal(ram_angular_difference(179,-179), 2)),
       "Angular difference must wrap across 180 degrees")
assert(is.na(ram_angular_difference(NA,-20)), "Missing angles must stay missing")
ref <- data.frame(
  chain="A",resi=1:4,insertion_code="",resn=c("ALA","SER","THR","GLY"),
  phi=c(-60,-75,179,-40),psi=c(-45,-30,-179,140),
  region=c("Favoured","Favoured","Allowed","Favoured")
)
other <- ref[c(1,2,3,4), ]
other <- other[-2, ,drop=FALSE]
other <- rbind(other[1L,,drop=FALSE],transform(other[1L,,drop=FALSE],
  resn="ASP",resi=88L,phi=-110,psi=85),
  other[2:3,,drop=FALSE])
pair <- ram_compare_torsions(ref,other)
assert(nrow(pair)==5L, "Alignment should preserve insertion and deletion")
assert(sum(pair$alignment=="Insertion")==1L, "Expected insertion")
assert(sum(pair$alignment=="Deletion")==1L, "Expected deletion")
assert(any(!is.na(pair$delta_phi)), "Aligned angles should be comparable")
assert(inherits(try(ram_align_residues(ref,other,max_cells=2),
                    silent=TRUE),"try-error"), "Bounded alignments")
empty <- ram_sequence_data(a[0,,drop=FALSE])
assert(nrow(empty)==0L, "Empty sequence should work")
pdb <- list(atom=data.frame(x=c(1,2),y=c(3,4),z=c(5,6)),
            xyz=matrix(c(1,3,5,2,4,6,
                         11,13,15,12,14,16),nrow=2,byrow=TRUE))
assert(ram_model_count(pdb)==2L, "Two-model Bio3D matrix")
m <- ram_model_at(pdb,2L)
assert(identical(m$atom$x,c(11,12)) &&
       identical(m$atom$z,c(15,16)), "Select model coordinates")
assert(ram_model_count(m)==1L, "Selected model must not retain ensemble flag")
assert(inherits(try(ram_model_at(pdb,3L),silent=TRUE),"try-error"),
       "Invalid model must fail")
message("Inspection, alignment and structural model unit tests passed")
