source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas-sifts.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
assert <- function(x,message) if(!isTRUE(x)) stop(message,call.=FALSE)

positions <- 1:40
coords <- cbind(positions*.35,sin(positions*.5)*3,cos(positions*.5)*3)
rotation <- matrix(c(0,-1,0, 1,0,0, 0,0,1),3,3)
rotated <- coords %*% rotation + matrix(c(25,-9,8),40,3,byrow=TRUE)
moved <- coords
moved[22:40,1] <- moved[22:40,1] + 8

fixture <- function(pdb,arr) {
  mapping <- data.frame(
    chain="A",resi=positions+100,insertion_code="",
    uniprot_accession="P12345",uniprot_resi=positions,
    entity_id=1L,struct_asym_id="X",
    label_seq_id=positions,observed=TRUE,
    canonical_source="PDBe updated mmCIF exact SIFTS",
    stringsAsFactors=FALSE)
  points <- lapply(seq_along(positions),function(i) list(
    struct_asym_id="X",label_seq_id=i,
    x=arr[i,1],y=arr[i,2],z=arr[i,3]))
  list(state="mapped",mapping=mapping,ca_points=points)
}
verified <- list("1ABC_1"=fixture("1ABC",coords),
  "2XYZ_1"=fixture("2XYZ",rotated),
  "3ABC_1"=fixture("3ABC",moved))

groups <- ram_atlas_geometry_groups(verified,"P12345",names(verified),
  min_core=30,min_fraction=.6,cutoff=1.5)
assert(length(groups$common_positions)==40L &&
       length(groups$sampled_positions)==40L,
       "Shared canonical coverage should be exactly forty residues.")
assert(groups$distance_matrix[1,2]<1e-9,
       "Rigid translation/rotation must not generate conformational distance.")
assert(groups$distance_matrix[1,3]>1.5,
       "An induced domain rearrangement must remain detectable.")
assert(groups$assignment$geometric_group[1] ==
       groups$assignment$geometric_group[2] &&
       groups$assignment$geometric_group[1] !=
       groups$assignment$geometric_group[3],
       "Exploratory clustering must separate shifted internal geometries.")
assert(length(groups$representatives)==2L,
       "Each geometric group should have one observed representative.")
assert(inherits(try(ram_atlas_geometry_groups(verified,"P12345",
  names(verified),min_core=41),silent=TRUE),"try-error"),
  "Insufficient common positions must block clustering.")
short <- verified
short[["3ABC_1"]]$mapping <- short[["3ABC_1"]]$mapping[1:20,]
assert(inherits(try(ram_atlas_geometry_groups(short,"P12345",
  names(short)),silent=TRUE),"try-error"),
  "A short construct must not distort common-core clustering.")
bad <- verified
bad[["2XYZ_1"]]$mapping$uniprot_accession <- "Q99999"
assert(inherits(try(ram_atlas_geometry_groups(bad,"P12345",
  names(bad)),silent=TRUE),"try-error"),
  "Mixed UniProt accessions must be rejected.")
# An ambiguous mapping is excluded from the common core.
conflicted <- verified
conflicted[["1ABC_1"]]$mapping <- rbind(
  conflicted[["1ABC_1"]]$mapping,
  transform(conflicted[["1ABC_1"]]$mapping[1,,drop=FALSE],
    uniprot_resi=80L))
clean <- ram_atlas_geometry_entities(conflicted,"P12345")
assert(nrow(clean[["1ABC_1"]]$coordinates)==39L &&
       !1L %in% clean[["1ABC_1"]]$coordinates$uniprot_resi,
       "Conflicting exact-SIFTS mapping must remove that coordinate.")
# Exercise the exact plotting function the Shiny view calls. Two entities
# need a direct distance chart, whereas three or more use a dendrogram.
pdf_path <- tempfile(fileext=".pdf")
grDevices::pdf(pdf_path)
two <- ram_atlas_geometry_groups(verified,"P12345",
  c("1ABC_1","3ABC_1"),cutoff=1.5)
ram_atlas_geometry_plot(two)
ram_atlas_geometry_plot(groups)
grDevices::dev.off()
assert(file.exists(pdf_path) && file.info(pdf_path)$size>100L,
       "Geometry charts must render for pairs and multi-entity cohorts.")
unlink(pdf_path)
message("Experimental Atlas distance-map grouping tests passed.")
