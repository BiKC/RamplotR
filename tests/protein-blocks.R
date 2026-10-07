# Pure Protein Blocks and coarse-state benchmark tests.
source(file.path("shinyRam","R","conformation.R"))
source(file.path("shinyRam","R","protein-blocks.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

make_pb_window <- function(reference) {
  data.frame(
    chain="A",resi=1:5,insertion_code="",resn="ALA",
    phi=c(NA,reference[[2L]],reference[[4L]],reference[[6L]],reference[[8L]]),
    psi=c(reference[[1L]],reference[[3L]],reference[[5L]],reference[[7L]],NA),
    bonded_to_next=c(TRUE,TRUE,TRUE,TRUE,FALSE),
    stringsAsFactors=FALSE
  )
}

for(block in rownames(ram_protein_block_references)) {
  window <- make_pb_window(ram_protein_block_references[block,])
  assigned <- ram_protein_blocks(window)
  assert(assigned$protein_block[[3L]]==block &&
         assigned$protein_block_rmsda[[3L]]<1e-10,
         paste("Exact Protein Block prototype did not assign to",block))
  assert(all(is.na(assigned$protein_block[c(1L,2L,4L,5L)])),
         "Terminal pentapeptide positions must remain unassigned.")
}

# Wrapped angles must preserve the nearest prototype.
m <- make_pb_window(ram_protein_block_references["m",])
m$phi[is.finite(m$phi)] <- m$phi[is.finite(m$phi)]+360
m$psi[is.finite(m$psi)] <- m$psi[is.finite(m$psi)]-360
wrapped <- ram_protein_blocks(m)
assert(wrapped$protein_block[[3L]]=="m" &&
       wrapped$protein_block_rmsda[[3L]]<1e-10,
       "Protein Block assignment must be periodic.")

# A broken peptide anywhere in the five-residue fragment invalidates it.
broken <- make_pb_window(ram_protein_block_references["d",])
broken$bonded_to_next[[2L]] <- FALSE
broken <- ram_protein_blocks(broken)
assert(is.na(broken$protein_block[[3L]]) &&
       !broken$protein_block_complete[[3L]],
       "Protein Blocks must not bridge peptide discontinuities.")

missing <- make_pb_window(ram_protein_block_references["d",])
missing$psi[[2L]] <- NA_real_
missing <- ram_protein_blocks(missing)
assert(is.na(missing$protein_block[[3L]]),
       "Missing required torsions must leave a Protein Block unassigned.")

# Equivalent-number calculation is 1 for one state and 2 for two equally
# frequent states.
assert(abs(ram_pb_equivalent_number(rep("m",10))-1)<1e-12,
       "Protein Block equivalent number changed for one state.")
assert(abs(ram_pb_equivalent_number(rep(c("m","d"),5))-2)<1e-12,
       "Protein Block equivalent number changed for two equal states.")

benchmark <- ram_pb_prototype_coarse_benchmark()
expected <- c(
  a="PPII",b="Alpha-R",c="Beta",d="Beta",e="Beta",f="PPII",
  g="Beta",h="PPII",i="Alpha-L",j="Other",k="Alpha-R",
  l="Alpha-R",m="Alpha-R",n="Alpha-R",o="Alpha-R",p="Alpha-L"
)
observed <- stats::setNames(benchmark$coarse_state,benchmark$protein_block)
assert(identical(unname(observed[names(expected)]),unname(expected)),
       "Coarse-state relationship to Protein Block prototypes changed.")
assert(benchmark$coarse_margin[benchmark$protein_block=="g"]<5 &&
       benchmark$coarse_margin[benchmark$protein_block=="m"]>100,
       paste0("Benchmark should reveal a coarse Beta/PPII boundary near PB g ",
              "while central alpha PB m remains stable."))

message("Protein Blocks fingerprint and coarse-state benchmark tests passed.")
