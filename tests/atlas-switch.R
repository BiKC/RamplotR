source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","backbone.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
source(file.path("shinyRam","R","atlas-switch.R"))
assert <- function(cond,msg) if(!isTRUE(cond)) stop(msg,call.=FALSE)

# The -180/180° seam is circular; unlike a naive subtraction, this is 2°.
p1 <- data.frame(uniprot_resi=101:104,
  phi=c(179,-179,-65,NA),psi=c(-179,179,125,NA))
p2 <- data.frame(uniprot_resi=101:104,
  phi=c(-179,179,5,NA),psi=c(179,-179,120,NA))
delta <- ram_atlas_torsion_delta(p1,p2,"1ABC_1","2XYZ_1",30)
assert(abs(delta$delta_phi[[1L]]-2)<1e-8 &&
       abs(delta$delta_psi[[1L]]+2)<1e-8 &&
       !delta$candidate[[1L]],"Circular seam interpreted as large shift.")
assert(delta$candidate[[3L]] &&
       is.na(delta$angular_shift[[4L]]),"Shift/missing handling failed.")

differences <- data.frame(uniprot_resi=c(5L,6L,7L,9L,10L,12L),
  candidate=c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE),
  angular_shift=c(35,45,3,70,40,50))
regions <- ram_atlas_candidate_regions(differences)
assert(nrow(regions)==3L &&
       identical(regions$start,c(5L,9L,12L)) &&
       identical(regions$end,c(6L,10L,12L)),
       "Candidate segments must split at missing/nonconsecutive UniProt positions.")

# N-CA-C points create unbroken, non-collinear peptide geometry.
n <- 40L
position <- seq_len(n)
ca <- cbind(1.5*position,0.5*sin(position),0.5*cos(position))
make_record <- function(shift=FALSE,gap=FALSE) {
  atom_rows <- list()
  for(i in position) {
    current <- ca[i,]
    N <- current + c(-0.2,0.3,0.1)
    C <- current + c(0.2,-0.1,0.25)
    if(shift && i %in% 15:19) {
      N <- N+c(0,0.45,-0.1)
      C <- C+c(0,-0.3,0.3)
    }
    if(gap && i==20L) N <- N+c(20,0,0)
    for(atom in c("N","CA","C")) {
      xyz <- switch(atom,N=N,CA=current,C=C)
      atom_rows[[length(atom_rows)+1L]] <- list(struct_asym_id="X",
        label_seq_id=i,atom_name=atom,
        x=xyz[[1L]],y=xyz[[2L]],z=xyz[[3L]])
    }
  }
  list(state="mapped",mapping=data.frame(
    chain="A",resi=position+100L,insertion_code="",
    uniprot_accession="P12345",uniprot_resi=position,
    entity_id=1L,struct_asym_id="X",label_seq_id=position,
    observed=TRUE,canonical_source="exact"),
    backbone_atoms=atom_rows)
}
a <- make_record()
b <- make_record(shift=TRUE)
tors_a <- ram_atlas_entity_torsions(a,"X")
tors_b <- ram_atlas_entity_torsions(b,"X")
assert(nrow(tors_a)==40L && nrow(tors_b)==40L &&
       sum(is.finite(tors_a$phi)&is.finite(tors_a$psi))>=36L,
       "Connected synthetic peptides should give torsion angles.")
same <- ram_atlas_torsion_delta(tors_a,tors_a,"A","A")
assert(all(same$angular_shift[is.finite(same$angular_shift)]<1e-9),
       "Identical backbones must have zero angular displacement.")
changed <- ram_atlas_torsion_delta(tors_a,tors_b,"A","B",30)
assert(any(changed$candidate[12:22]),
       "Changed backbone N/C geometry should produce local torsion shift.")
gap <- make_record(gap=TRUE)
gapped <- ram_atlas_entity_torsions(gap,"X")
assert(is.na(gapped$psi[[19L]]) && is.na(gapped$phi[[20L]]),
       "Broken C-N peptide continuity must invalidate both adjacent angles.")
conflicted <- make_record()
conflicted$mapping <- rbind(conflicted$mapping,
  transform(conflicted$mapping[5L,,drop=FALSE],uniprot_resi=80L))
reduced <- ram_atlas_entity_torsions(conflicted,"X")
assert(!5L %in% reduced$uniprot_resi &&
       !80L %in% reduced$uniprot_resi,
       "Conflicting exact SIFTS targets may not contribute torsion evidence.")

# Mimic two different geometry-cluster medoids without inventing populations.
verified <- list("1ABC_1"=a,"2XYZ_1"=b)
geometry <- list(accession="P12345",
  representatives=c("1"="1ABC_1","2"="2XYZ_1"),
  assignment=data.frame(entity=c("1ABC_1","2XYZ_1")))
# Geometry data helper requires C-alpha records in the verified payloads.
for(id in names(verified)) {
  source_atoms <- verified[[id]]$backbone_atoms
  verified[[id]]$ca_points <- lapply(Filter(function(x)
    identical(x$atom_name,"CA"),source_atoms),function(x) {
    x$atom_name <- NULL
    x
  })
}
report <- ram_atlas_group_switches(verified,geometry,30)
assert(identical(report$representatives,c("1ABC_1","2XYZ_1")) &&
       report$comparable>=35L &&
       any(report$residues$candidate),
       "Representative comparison must report local torsion evidence.")
assert(all(c("chain_a","resi_a","insertion_a","label_seq_a",
  "chain_b","resi_b","insertion_b","label_seq_b") %in% names(report$residues)),
  "Atlas residue inspection must retain exact local residue identifiers.")
assert(all(report$residues$resi_a[!is.na(report$residues$resi_a)] ==
  report$residues$uniprot_resi[!is.na(report$residues$resi_a)]+100L),
  "Exact SIFTS author numbering should be preserved.")
reversed <- ram_atlas_group_switches(verified,geometry,30,
  representative_ids=c("2XYZ_1","1ABC_1"))
assert(identical(reversed$representatives,c("2XYZ_1","1ABC_1")) &&
  isTRUE(all.equal(reversed$residues$delta_phi,-report$residues$delta_phi,
    check.attributes=FALSE)),
  "Reversing representatives should reverse wrapped phi change.")
invalid <- try(ram_atlas_group_switches(verified,geometry,30,
  representative_ids=c("1ABC_1","1ABC_1")),silent=TRUE)
assert(inherits(invalid,"try-error"),
  "The same group representative must not be selected twice.")
# Atlas-linked Compare focus must use both exact PDB author positions,
# not the UniProt number or the first sequence-alignment candidate.
aligned <- data.frame(chain_a=c("A","A","A"),residue_a=c(102L,103L,103L),
  insertion_a=c("","","A"),chain_b=c("B","B","B"),
  residue_b=c(202L,203L,203L),insertion_b=c("","",""),
  row_id=1:3)
target <- data.frame(chain_a="A",resi_a=103L,insertion_a="A",
  chain_b="B",resi_b=203L,insertion_b="")
assert(identical(ram_atlas_comparison_pair_index(aligned,target),3L),
  "Atlas focus must match the exact insertion-bearing author pair.")
not_both <- transform(target,resi_b=204L)
assert(is.na(ram_atlas_comparison_pair_index(aligned,not_both)),
  "Wrong second residue must not be selected by coincidence.")
no_mapping <- transform(target,insertion_a=NA_character_)
assert(is.na(ram_atlas_comparison_pair_index(aligned,no_mapping)),
  "Missing insertion-code mapping must reject automatic focus.")
duplicate <- rbind(aligned,aligned[3L,,drop=FALSE])
assert(is.na(ram_atlas_comparison_pair_index(duplicate,target)),
  "Ambiguous alignment pairs must not be auto-selected.")
# Same verified UniProt positions, but substantially different author numbering.
# Sequence offsets must never be used to build this alignment.
second_map <- transform(b$mapping,chain="B",resi=resi+900L)
second_record <- b
second_record$mapping <- second_map
first_torsions <- ram_atlas_entity_torsions(a,"X")
second_torsions <- ram_atlas_entity_torsions(second_record,"X")
canonical <- ram_atlas_canonical_pairing(first_torsions,second_torsions,
  a$mapping,second_map,"P12345","A","X","B","X")
assert(canonical$matched==40L &&
       identical(canonical$pairing$uniprot_resi,1:40) &&
       identical(canonical$pairing$index_a,1:40) &&
       identical(canonical$pairing$index_b,1:40),
  "Canonical alignment should preserve exact matched UniProt positions across numbering changes.")
assert(first_torsions$resi[[10L]] != second_torsions$resi[[10L]],
  "Test must use different author residue numbering.")

# Non-overlapping accession and ambiguous one-to-many mapping cannot silently
# turn into positional matches.
foreign_map <- transform(second_map,uniprot_accession="Q99999")
foreign <- ram_atlas_canonical_pairing(first_torsions,second_torsions,
  a$mapping,foreign_map,"P12345","A","X","B","X")
assert(foreign$matched==0L,
  "Different accession must not produce canonical matches.")
duplicate_map <- rbind(second_map,
  transform(second_map[10L,,drop=FALSE],uniprot_resi=900L))
ambiguous <- ram_atlas_canonical_pairing(first_torsions,second_torsions,
  a$mapping,duplicate_map,"P12345","A","X","B","X")
assert(ambiguous$matched==39L &&
       !10L %in% ambiguous$pairing$uniprot_resi,
  "Conflicting canonical SIFTS residue must be excluded entirely.")
label_conflict <- second_map
label_conflict$label_seq_id[[11L]] <- label_conflict$label_seq_id[[10L]]
label_ambiguous <- ram_atlas_canonical_pairing(first_torsions,second_torsions,
  a$mapping,label_conflict,"P12345","A","X","B","X")
assert(label_ambiguous$matched==38L &&
       !any(c(10L,11L) %in% label_ambiguous$pairing$uniprot_resi),
  "Duplicate label_seq_id must exclude both conflicting canonical residues.")

# Insertion codes are part of the local author identifier.
insertion_map <- a$mapping
insertion_map$insertion_code[[10L]] <- "A"
wrong_insertion <- ram_atlas_canonical_pairing(first_torsions,second_torsions,
  insertion_map,second_map,"P12345","A","X","B","X")
assert(wrong_insertion$matched==39L,
  "PDB insertion codes must match observed torsion records exactly.")
inserted_torsions <- first_torsions
inserted_torsions$insertion_code[[10L]] <- "A"
correct_insertion <- ram_atlas_canonical_pairing(inserted_torsions,
  second_torsions,insertion_map,second_map,
  "P12345","A","X","B","X")
assert(correct_insertion$matched==40L,
  "A genuinely observed insertion code should still pair by UniProt.")

# No assumed gaps: only residues verified on both sides are compared.
truncated <- second_torsions[-(1:7),,drop=FALSE]
partial <- ram_atlas_canonical_pairing(first_torsions,truncated,
  a$mapping,second_map,"P12345","A","X","B","X")
assert(partial$matched==33L &&
       identical(partial$pairing$uniprot_resi,8:40) &&
       partial$total_a==40L && partial$total_b==33L,
  "Partial constructs must report the common mapped core and true denominators.")

message("Experimental Atlas backbone switch-region tests passed.")
