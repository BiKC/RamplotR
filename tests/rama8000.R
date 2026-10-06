# Run from repository root: Rscript tests/rama8000.R
source(file.path("shinyRam","R","rama8000.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

groups <- ram_rama8000_group(
  resn=c("ALA","GLY","PRO","PRO","VAL","ILE","SER"),
  next_resn=c("PRO","ALA","ALA","ALA","ALA","PRO","ALA"),
  bonded_to_next=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
  omega_prev=c(180,180,0,179,180,180,180)
)
assert(identical(groups,c(
  "pre-proline","glycine","cis-proline","trans-proline",
  "isoleucine or valine","pre-proline","general"
)), "Six-class Rama8000 residue grouping differs from cctbx priority")

dir <- file.path("shinyRam","static","rama8000")

# Regression points are taken from current cctbx tst_ramalyze.py output.
cases <- data.frame(
  group=c("general","general","pre-proline","trans-proline",
          "isoleucine or valine","general","general","general"),
  phi=c(-83.26,-111.53,-42.39,-39.12,-60.38,-61.13,60.09,-37.21),
  psi=c(131.88,71.36,121.87,-31.84,-51.85,-170.23,-80.26,-36.12),
  score_pct=c(35.07,0.74,2.66,0.31,68.24,0.02,0.02,0.13),
  region=c("Favored","Allowed","Favored","Allowed",
           "Favored","Outlier","Outlier","Allowed"),
  stringsAsFactors=FALSE
)
observed <- numeric(nrow(cases))
for(i in seq_len(nrow(cases))) {
  table <- ram_rama8000_table(dir,cases$group[[i]])
  observed[[i]] <- 100 * ram_rama8000_score_one(
    table,cases$phi[[i]],cases$psi[[i]])
}
assert(max(abs(observed-cases$score_pct)) < 0.015,
       paste("Rama8000 interpolation differs from cctbx:",
             paste(round(observed,3),collapse=", ")))
regions <- ram_rama8000_region(observed/100,cases$group)
assert(identical(regions,cases$region),
       "Rama8000 score thresholds differ from cctbx")

# Periodic angle handling must match at the -180/180 boundary.
general <- ram_rama8000_table(dir,"general")
a <- ram_rama8000_score_one(general,181,-179)
b <- ram_rama8000_score_one(general,-179,-179)
assert(isTRUE(all.equal(a,b,tolerance=1e-12)),
       "Rama8000 periodic interpolation changed across phi boundary")

torsions <- data.frame(
  resn=c("ALA","PRO","VAL"), next_resn=c("PRO","VAL",NA),
  bonded_to_next=c(TRUE,TRUE,FALSE), omega_prev=c(NA,5,175),
  phi=c(-60,-65,-60), psi=c(-45,145,-52),
  stringsAsFactors=FALSE
)
classified <- ram_rama8000_classify(torsions,dir)
assert(all(c("rama8000_group","rama8000_score","rama8000_region") %in%
             names(classified)), "Rama8000 output schema changed")
assert(identical(classified$rama8000_group,
  c("pre-proline","cis-proline","isoleucine or valine")),
  "Integrated Rama8000 grouping changed")

message("Rama8000 six-class validation regression tests passed")
