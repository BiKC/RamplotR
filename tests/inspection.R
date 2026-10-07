# Pure unit coverage; no Shiny session or network required.
source(file.path("shinyRam","R","conformation.R"))
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
assert(identical(q$review_status, c("Not allowed","Missing angles","Near boundary","Other")),
       "Review queue priority and missing angles")
a$rama8000_region <- c("Favored",NA,"Outlier","Allowed")
q_standard <- ram_review_queue(a)
assert(identical(q_standard$review_status,
  c("Rama8000 outlier","Missing angles","Near boundary","Other")),
  "Rama8000 outliers should be reviewed separately from RamplotR regions")
assert(identical(q$resn, c("GLY","THR","PRO","ALA")),
       "Review queue ordering")
seq <- ram_sequence_data(a)
assert(identical(seq$letter, c("A","T","G","P")),
       "Sequence must retain amino-acid identities")
assert(identical(unname(ram_plddt_color(c(NA,49,50,69,70,89,90,100))),
  c("#cbd7db","#d75e56","#d6ac52","#d6ac52",
    "#7bbcb1","#7bbcb1","#126e74","#126e74")),
  "Confidence colour thresholds must be independent of region labels")
assert(identical(ram_sequence_position_labels(c(98L,99L,100L,101L,104L,105L)),
  c("98","","100","","","105")),
  "Sequence index labels must use real residue numbering, not array positions")
assert(identical(ram_sequence_position_labels(c(101L,101L), c("","A")),
  c("101","101A")), "Labels retain insertion codes")
assert(identical(ram_sequence_position_labels(integer(), character()),
  character()), "Empty position labels are supported")
assert(seq$insertion_code[2L]=="A", "Insertion code must be preserved")
groups <- ram_sequence_groups(a)
assert(identical(names(groups), c("A", "B")),
       "Navigator must retain every selected protein chain in source order")
assert(identical(vapply(groups, nrow, integer(1)), c(A=3L, B=1L)),
       "Navigator must preserve all residues of each chain")
assert(length(ram_sequence_groups(a[0,,drop=FALSE]))==0L,
       "Empty chain selection should show no groups")
overview <- ram_sequence_overview_bins(
  c(rep("Favoured", 70L), "Not allowed", rep("Favoured", 70L),
    NA_character_), max_bins=12L)
assert(length(overview)==12L,
       "Large chains should use a bounded number of overview cells")
assert("outlier" %in% overview && "missing" %in% overview,
       "Compressed overviews must not lose isolated issues")

assert(isTRUE(all.equal(ram_angular_difference(179,-179), 2)),
       "Angular difference must wrap across 180 degrees")
assert(is.na(ram_angular_difference(NA,-20)), "Missing angles must stay missing")
shift <- ram_backbone_angular_displacement(c(3,18,30,60,NA),
                                           c(4,0,0,0,2))
assert(isTRUE(all.equal(shift[1:4],c(5,18,30,60))) && is.na(shift[[5]]),
       "Backbone angular displacement changed")
assert(identical(ram_backbone_shift_band(shift),
  c("Small","Moderate","Large","Very large","Unavailable")),
  "Backbone shift display bands changed")
ref <- data.frame(
  chain="A",resi=1:4,insertion_code="",resn=c("ALA","SER","THR","GLY"),
  phi=c(-60,-75,179,-40),psi=c(-45,-30,-179,140),
  region=c("Favoured","Favoured","Allowed","Favoured")
)
other <- rbind(ref[1,,drop=FALSE], transform(ref[1,,drop=FALSE],
  resn="ASP",resi=88L,phi=-110,psi=85), ref[2:4,,drop=FALSE])
pair <- ram_compare_torsions(ref,other)
assert(nrow(pair)==5L, "Alignment should preserve inserted residues")
assert(all(c("basin_a","basin_b","basin_changed") %in% names(pair)),
       "Pairwise comparison must expose coarse backbone-state transitions")
state_ref <- ref[1,,drop=FALSE]
state_other <- state_ref
state_ref$phi <- -63; state_ref$psi <- -43
state_other$phi <- -135; state_other$psi <- 135
state_pair <- ram_compare_torsions(state_ref,state_other)
assert(state_pair$basin_a=="Alpha-R" && state_pair$basin_b=="Beta" &&
       isTRUE(state_pair$basin_changed),
       "Pairwise comparison must flag a broad backbone-state transition")
assert(sum(pair$alignment=="Insertion")==1L, "Expected insertion")
quality <- ram_comparison_alignment_quality(pair)
assert(quality$aligned==4L && quality$matches==4L &&
       isTRUE(all.equal(quality$identity,1)) &&
       isTRUE(all.equal(quality$coverage_a,1)) &&
       isTRUE(all.equal(quality$coverage_b,0.8)),
       "Alignment quality must report identity and both chain coverages.")
substitution <- ref
substitution$resn[[2L]] <- "ASP"
sub_quality <- ram_comparison_alignment_quality(
  ram_compare_torsions(ref,substitution))
assert(sub_quality$substitutions==1L &&
       isTRUE(all.equal(sub_quality$identity,0.75)),
       "Aligned substitutions must reduce sequence identity without reducing coverage.")
assert(identical(ram_comparison_find(pair,"a","A",2L), 3L),
  "Find the primary residue by exact PDB numbering")
assert(identical(ram_comparison_find(pair,"b","A",88L), 2L),
  "Locate inserted comparison-only residue")
assert(is.na(ram_comparison_find(pair,"a","A",88L)),
  "Do not assign the partner's residue number to a primary gap")
assert(is.na(ram_comparison_find(pair,"b","A",999L)),
  "Unknown comparison residue must not select an unrelated pair")
deletion <- ram_compare_torsions(ref, ref[-2,,drop=FALSE])
assert(sum(deletion$alignment=="Deletion")==1L, "Expected deletion")
assert(any(!is.na(pair$delta_phi)), "Aligned angles should be comparable")
pred_a <- ref
pred_b <- ref
pred_a$plddt <- c(95,80,72,60)
pred_b$plddt <- c(75,82,92,60)
pred_a$confidence_category <- c("Very high","Confident","Confident","Low")
pred_b$confidence_category <- c("Confident","Confident","Very high","Low")
confidence_pair <- ram_compare_torsions(pred_a,pred_b)
assert(isTRUE(all.equal(confidence_pair$delta_plddt,c(-20,2,20,0))) &&
       identical(confidence_pair$confidence_changed,
                 c(TRUE,FALSE,TRUE,FALSE)),
       "Pairwise prediction confidence differences changed.")
mixed_pair <- ram_compare_torsions(ref,pred_b)
assert(all(is.na(mixed_pair$plddt_a)) &&
       identical(mixed_pair$plddt_b,pred_b$plddt) &&
       all(is.na(mixed_pair$delta_plddt)),
       "Experimental-vs-predicted comparison must keep missing confidence missing.")
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
seq_fixture <- data.frame(
  chain=c("A","A","A","B"),
  resi=c(1L,2L,3L,1L),
  insertion_code=c("","","",""),
  resn=c("ALA","GLY","MSE","VAL"),
  phi=c(-60,-70,-80,-90),
  psi=c(-40,140,150,160),
  region=c("Favoured","Allowed","Allowed","Favoured"),
  stringsAsFactors=FALSE
)
seq_a <- ram_chain_query_sequence(seq_fixture,"A")
assert(identical(seq_a$sequence,"AGX"),
       "Counterpart search must preserve unknown residues as X.")
assert(isTRUE(all.equal(seq_a$known_fraction,2/3)),
       "Known-residue fraction changed for counterpart sequence search.")
assert(identical(ram_chain_query_sequence(seq_fixture,"B")$sequence,"V"),
       "Chain-specific counterpart sequence extraction changed.")


context_pdb <- list(
  atom=data.frame(
    type="ATOM",chain=c("A","A"),resno=c(10L,10L),resid=c("ALA","ALA"),
    insert=c("",""),elety=c("CA","CB"),elesy=c("C","C"),
    x=c(0,1),y=c(0,0),z=c(0,0),stringsAsFactors=FALSE
  ),
  hetero_atom=data.frame(
    type="HETATM",chain=c("B","B","B","B"),resno=c(401L,401L,402L,501L),
    resid=c("ATP","ATP","MG","HOH"),insert=c("","","",""),
    elety=c("P","H1","MG","O"),elesy=c("P","H","MG","O"),
    x=c(4,2,5,2),y=c(0,0,0,0),z=c(0,0,0,0),stringsAsFactors=FALSE
  )
)
context_row <- data.frame(chain="A",resi=10L,insertion_code="",resn="ALA",
                          stringsAsFactors=FALSE)
nearby <- ram_nearby_hetero_context(context_pdb,context_row,max_distance=6)
assert(nrow(nearby)==2L && identical(as.character(nearby$resn),c("ATP","MG")),
       "Nearby hetero context must exclude water and group atoms by hetero residue.")
assert(isTRUE(all.equal(nearby$distance,c(3,4))) &&
       identical(as.character(nearby$target_atom),c("CB","CB")),
       "Nearest hetero distances must use heavy atoms from the selected residue.")
assert(nrow(ram_nearby_hetero_context(
  context_pdb,context_row,max_distance=2.5))==0L,
  "Distant hetero residues must stay out of the local context.")

evidence_row <- data.frame(
  chain="A",resi=10L,insertion_code="",resn="ALA",
  phi=-120,psi=95,region="Not allowed",density=98.4,
  rama8000_region="Outlier",rama8000_group="general",
  rama8000_score=0.0001,plddt=96,confidence_category="Very high",
  omega=72,omega_status="Twisted",
  wwpdb_rama="outlier",wwpdb_rotamer="outlier",
  wwpdb_clashes=2,wwpdb_bond_outliers=1,wwpdb_angle_outliers=0,
  stringsAsFactors=FALSE
)
ev <- ram_residue_evidence(evidence_row,local_context=nearby)
assert(all(c(
  "Rama8000 backbone outlier",
  "Unusual RamplotR density position",
  "High-confidence prediction with unusual backbone geometry",
  "Twisted peptide bond",
  "Official wwPDB Ramachandran outlier",
  "Official wwPDB rotamer outlier",
  "Official wwPDB local clash",
  "Official covalent-geometry outlier",
  "Nearby non-water hetero residues"
) %in% ev$title), "Residue evidence should preserve independent warning sources.")
assert(ev$title[[1L]] %in% c(
  "High-confidence prediction with unusual backbone geometry",
  "Official wwPDB Ramachandran outlier","Twisted peptide bond"),
  "High-severity residue evidence should sort before informational notes.")
quiet_row <- evidence_row
quiet_row$region <- "Favoured"; quiet_row$density <- 40
quiet_row$rama8000_region <- "Favored"; quiet_row$plddt <- 95
quiet_row$omega <- 179; quiet_row$omega_status <- "Trans"
quiet_row$wwpdb_rama <- "favored"; quiet_row$wwpdb_rotamer <- "favored"
quiet_row$wwpdb_clashes <- 0; quiet_row$wwpdb_bond_outliers <- 0
quiet_row$wwpdb_angle_outliers <- 0
assert(nrow(ram_residue_evidence(quiet_row))==0L,
       "Unremarkable residues should not receive artificial warnings.")

ensemble_context <- data.frame(
  models_present=4L,ensemble_models_total=5L,
  phi_sd=34,psi_sd=12,
  rama8000_changes=TRUE,rama8000_consistency=0.75,
  basin_changes=TRUE,basin_consistency=0.75,basin_mode="Alpha-R",
  plddt_mean=94,plddt_sd=12,
  stringsAsFactors=FALSE
)
ensemble_ev <- ram_residue_evidence(
  quiet_row,ensemble_context=ensemble_context)
assert(all(c(
  "High-confidence predictions disagree on local backbone",
  "Rama8000 category differs across prediction models",
  "Prediction models choose different backbone states",
  "Prediction confidence varies across ensemble",
  "Residue is absent from some prediction models"
) %in% ensemble_ev$title),
  "Prediction ensemble evidence should explain geometry, standard-category, confidence and coverage disagreement.")
assert(ensemble_ev$title[[1L]]==
         "High-confidence predictions disagree on local backbone",
       "High-confidence backbone disagreement should be prioritised.")
assert(all(ensemble_ev$source=="Prediction ensemble"),
       "Ensemble-derived evidence must retain its source.")

missing_row <- quiet_row
missing_row$phi <- NA_real_
missing_ev <- ram_residue_evidence(missing_row)
assert("Backbone angles unavailable" %in% missing_ev$title,
       "Missing torsions should be explained explicitly.")

message("Inspection, alignment and structural model unit tests passed")
