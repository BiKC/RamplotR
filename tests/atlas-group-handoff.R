# Atlas -> Group Compare direct-input scientific regression.
# No simulated files, external downloads or unverified sequence alignment.
source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","backbone.R"))
source(file.path("shinyRam","R","conformation.R"))
source(file.path("shinyRam","R","inspection.R"))
source(file.path("shinyRam","R","ensemble.R"))
source(file.path("shinyRam","R","group-comparison.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
source(file.path("shinyRam","R","atlas-switch.R"))
source(file.path("shinyRam","R","atlas-construct.R"))
source(file.path("shinyRam","R","atlas-group-handoff.R"))
assert <- function(x,msg) if(!isTRUE(x))stop(msg,call.=FALSE)

positions <- seq_len(45L)
make_record <- function(offset=100L,shift=0,modification=FALSE) {
  atoms <- list()
  ca <- cbind(1.5*positions,0.45*sin(positions),0.45*cos(positions))
  for(i in positions) {
    point <- ca[i,]
    n <- point+c(-0.2,0.3+if(i%in%15:22)shift else 0,0.1)
    c_atom <- point+c(0.2,-0.1-if(i%in%15:22)shift else 0,0.25)
    for(atom in c("N","CA","C")) {
      coord <- switch(atom,N=n,CA=point,C=c_atom)
      atoms[[length(atoms)+1L]] <- list(
        struct_asym_id="X",label_seq_id=i,atom_name=atom,
        x=coord[[1L]],y=coord[[2L]],z=coord[[3L]])
    }
  }
  mon <- rep("ALA",45L)
  if(modification) mon[[20L]] <- "MSE"
  list(state="mapped",mapping=data.frame(
    chain="A",resi=positions+offset,insertion_code="",
    uniprot_accession="P12345",uniprot_resi=positions,
    entity_id=1L,struct_asym_id="X",label_seq_id=positions,
    observed=TRUE,mon_id=mon,canonical_source="exact"),
    backbone_atoms=atoms,
    ca_points=lapply(Filter(function(x) identical(x$atom_name,"CA"),atoms),
      function(x){x$atom_name<-NULL;x}))
}
ids <- c("1ABC_1","2ABC_1","3ABC_1","4ABC_1")
verified <- list(
  "1ABC_1"=make_record(100,0),
  "2ABC_1"=make_record(950,0.02),
  "3ABC_1"=make_record(300,0.25),
  "4ABC_1"=make_record(1800,0.27,TRUE))
geometry <- list(accession="P12345",selected=ids,
  common_positions=positions,
  assignment=data.frame(entity=ids,geometric_group=c(1L,1L,2L,2L)),
  error=NULL)
picked <- ram_atlas_group_select(geometry,ids[1:2],ids[3:4])
assert(identical(picked$group_a,ids[1:2]) &&
       identical(picked$group_b,ids[3:4]),
  "Group selection must preserve researcher assignments.")
for(groups in list(list(ids[1:2],ids[2:3]),
                   list(ids[1],ids[1]),list(ids[1],character()),
                   list(ids[1],"9ABC_1"))) {
  assert(inherits(try(ram_atlas_group_select(geometry,
    groups[[1L]],groups[[2L]]),silent=TRUE),"try-error"),
    "Overlapping, empty or unverified groups must be rejected.")
}
out <- ram_atlas_group_prepare(verified,geometry,ids[1:2],ids[3:4],
  "Enzyme set A","Enzyme set B")
assert(identical(out$source,"atlas") &&
       out$n_a==2L && out$n_b==2L &&
       identical(out$label_a,"Enzyme set A") &&
       out$core_positions==45L,
  "Direct Atlas source should retain group memberships and labels.")
assert(nrow(out$comparison)==45L &&
       all(out$comparison$chain=="UniProt") &&
       identical(sort(out$comparison$resi),positions),
  "Group output must use exact UniProt positions rather than PDB author numbering.")
assert(any(is.finite(out$comparison$angular_displacement)) &&
       out$comparison$a_phi_models[[which(out$comparison$resi==19L)]]==2L &&
       out$comparison$b_phi_models[[which(out$comparison$resi==19L)]]==2L,
  "Paired experimental torsions must enter circular group statistics.")
assert(all(is.na(out$comparison$a_rama8000_mode)) &&
       all(is.na(out$comparison$b_rama8000_mode)),
  "Do not invent Rama8000 classifications from Atlas backbone-only records.")
assert(out$unknown_chemistry_positions==1L &&
       out$known_chemistry_differences,
  "Known MSE-versus-ALA chemistry conflict needs visible uncertainty.")
assert(out$comparison$resn[out$comparison$resi==20L]=="UNK" &&
       out$comparison$resn[out$comparison$resi==19L]=="ALA",
  "Only unanimous observed experimental monomers can label a canonical position.")
assert(nrow(out$members)==4L &&
       all(out$members$mapping=="Exact observed PDBe SIFTS UniProt") &&
       all(out$members$complete_torsions>=40L),
  "Exact canonical mapping and coverage must be preserved in member export.")

# A missing torsion on one experimental member reduces support, but cannot
# silently add interpolated UniProt coordinates.
missing <- verified
missing[["4ABC_1"]]$backbone_atoms <-
  Filter(function(x) x$label_seq_id!=22L,
    missing[["4ABC_1"]]$backbone_atoms)
sparse <- ram_atlas_group_prepare(missing,geometry,ids[1:2],ids[3:4])
row <- sparse$comparison[sparse$comparison$resi==22L,,drop=FALSE]
assert(nrow(row)==1L && row$b_phi_models[[1L]]<=1L,
  "Missing exact backbone atoms must decrease group support.")

# No loaded/reference structure is needed for a verified Atlas cohort.
single <- ram_atlas_group_prepare(verified,geometry,ids[[1L]],ids[[3L]])
assert(single$n_a==1L && single$n_b==1L &&
       !any(single$comparison$high_support_shift,na.rm=TRUE),
  "A one-vs-one exploratory comparison must not claim high-support shifts.")

cat("Exact-SIFTS Atlas -> Compare Groups scientific tests passed.\n")
