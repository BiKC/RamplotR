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
auto <- ram_atlas_geometry_groups(verified,"P12345",names(verified),
  min_core=30,min_fraction=.6,cluster_mode="automatic")
assert(identical(auto$cluster_mode,"automatic") &&
       auto$auto_cluster$k==2L &&
       length(auto$representatives)==2L,
  "Automatic clustering should separate the induced domain-motion example.")
assert(auto$assignment$geometric_group[[1L]]==
       auto$assignment$geometric_group[[2L]] &&
       auto$assignment$geometric_group[[1L]]!=
       auto$assignment$geometric_group[[3L]],
  "Automatic grouping must be rigid-motion invariant on verified UniProt coordinates.")
assert(grepl("Suggested",auto$auto_cluster$reason) &&
       nrow(auto$auto_cluster$candidates)==1L &&
       auto$auto_cluster$candidates$singleton_groups[[1L]]==1L,
  "Automatic cluster summary must disclose its lone singleton.")
assert(all(auto$distance_matrix==groups$distance_matrix),
  "Automatic mode must reuse the same fixed-core structural distances.")
flat <- matrix(1,5L,5L,dimnames=list(LETTERS[1:5],LETTERS[1:5]))
diag(flat) <- 0
no_split <- ram_atlas_auto_clusters(flat,
  stats::hclust(stats::as.dist(flat),method="average"))
assert(no_split$k==1L && all(no_split$assignment==1L) &&
       grepl("No sufficiently separated",no_split$reason),
  "Equidistant experimental structures must not be forced into clusters.")
close <- flat*0.05
close_split <- ram_atlas_auto_clusters(close,
  stats::hclust(stats::as.dist(close),method="average"))
assert(close_split$k==1L,
  "Small numerical coordinate differences cannot force biological-looking groups.")
two_auto <- ram_atlas_geometry_groups(verified,"P12345",
  c("1ABC_1","3ABC_1"),cluster_mode="automatic")
assert(two_auto$auto_cluster$k==1L &&
       grepl("Fewer than three",two_auto$auto_cluster$reason),
  "Two PDBs cannot provide within-cluster support for an automatic two-group recommendation.")
# Two tight pairs should produce two suggested groups rather than over-
# partitioning into singletons, with deterministic results on repeated calls.
paired <- matrix(c(0,0.1,2,2,0.1,0,2,2,
                   2,2,0,0.1,2,2,0.1,0),4,4,byrow=TRUE,
                 dimnames=list(LETTERS[1:4],LETTERS[1:4]))
auto_paired <- ram_atlas_auto_clusters(paired,
  stats::hclust(stats::as.dist(paired),method="average"))
assert(auto_paired$k==2L &&
       auto_paired$assignment[[1L]]==auto_paired$assignment[[2L]] &&
       auto_paired$assignment[[3L]]==auto_paired$assignment[[4L]] &&
       auto_paired$assignment[[1L]]!=auto_paired$assignment[[3L]] &&
       all(auto_paired$candidates$singleton_groups>=0),
  "Automatic grouping should find two coherent replicated structural groups.")
assert(identical(auto_paired$assignment,
  ram_atlas_auto_clusters(paired,
    stats::hclust(stats::as.dist(paired),method="average"))$assignment),
  "Recommended groups must be deterministic for the same input distances.")
assert(inherits(try(ram_atlas_auto_clusters(matrix(c(0,2,1,0),2),
  stats::hclust(stats::as.dist(matrix(c(0,2,1,0),2)))),
  silent=TRUE),"try-error"),
  "Asymmetric pairwise distance matrices must not be accepted.")
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
