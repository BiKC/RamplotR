# Scientific regression for per-member group fingerprints.
source(file.path("shinyRam","R","conformation.R"))
source(file.path("shinyRam","R","ensemble.R"))
source(file.path("shinyRam","R","group-fingerprint.R"))
assert <- function(x,message) if(!isTRUE(x)) stop(message,call.=FALSE)
member <- function(phi,psi,aa="ALA") data.frame(
  chain="A",resi=10L,insertion_code="",resn="ALA",
  source_resn=aa,phi=phi,psi=psi,stringsAsFactors=FALSE)
A <- list("apo-1"=member(179,-179),
          "apo-2"=member(-179,179,aa="VAL"))
B <- list("holo-1"=member(-63,-43),
          "holo-2"=member(NA_real_,-40))
observations <- ram_group_fingerprint(A,B,"apo","holo")
assert(nrow(observations)==4L && sum(observations$paired)==3L,
  "Each member must contribute one row, but only complete pairs count.")
residue <- ram_group_fingerprint_at(observations,"A",10L,"")
assert(nrow(residue)==4L && identical(sort(unique(residue$amino_acid)),
  c("ALA","VAL")),"The actual candidate residue identity must be retained.")
summary <- ram_group_fingerprint_summary(residue,c("apo","holo"),c(2L,2L))
assert(summary$complete_pairs[[1L]]==2L &&
       summary$complete_pairs[[2L]]==1L &&
       summary$members[[2L]]==2L,
  "Incomplete pairs must not contribute to within-group support.")
assert(abs(abs(summary$phi_mean[[1L]])-180)<1e-8,
  "Circular group mean must handle the ±180-degree seam.")
assert(is.na(summary$phi_sd[[2L]]),
  "A group with one paired observation has no measured dispersion.")
assert(nrow(ram_group_fingerprint_at(observations,"A",11L,""))==0L,
  "Unobserved residue positions cannot be fabricated.")
bad <- A
bad[[1L]] <- rbind(bad[[1L]],bad[[1L]])
assert(inherits(try(ram_group_fingerprint(bad,B),silent=TRUE),
               "try-error"),
  "Ambiguous residue identifiers must not silently inflate fingerprints.")
fig <- tempfile(fileext=".pdf")
grDevices::pdf(fig)
ram_group_fingerprint_plot(residue,c("apo","holo"))
grDevices::dev.off()
assert(file.exists(fig) && file.info(fig)$size>100L,
  "Fingerprint plot must render from measured member angles.")
unlink(fig)
message("Group conformational fingerprint tests passed.")
