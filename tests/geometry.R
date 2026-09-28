# Reproducible Phase C native diagnostic geometry.
source("shinyRam/R/backbone.R")
source("shinyRam/R/geometry.R")
assert <- function(x, message) if(!isTRUE(x)) stop(message,call.=FALSE)
stopifnot(identical(ram_geom_omega_status(c(NA,-179,0,180,74)),
                    c("Missing","Trans","Cis","Trans","Twisted")))
# Synthetic first residue has a complete chi1, peptide omega to the second.
residue <- function(chain, no, name, insert, atoms, coords) {
  data.frame(chain=chain,resno=no,resid=name,insert=insert,
             elety=atoms,alt="",x=coords[,1],y=coords[,2],z=coords[,3],
             stringsAsFactors=FALSE)
}
first <- residue("A",1,"SER","",c("N","CA","C","CB","OG"),
  rbind(c(0,0,0),c(1,0,0),c(1,1,0),c(1,-1,0),c(2,-1,1)))
second <- residue("A",2,"PRO","",c("N","CA","C","CB","CG"),
  rbind(c(1,2.4,0),c(2,2.4,1),c(3,2.4,1),c(2,1.4,1),c(3,1.4,2)))
pdb <- list(atom=rbind(first,second))
torsion <- ram_extract_torsions(pdb)
assert(nrow(torsion)==2L && torsion$bonded_to_next[[1]],
       "Peptide bond in synthetic geometry was lost.")
geom <- ram_extra_geometry(pdb,torsion)
assert(nrow(geom)==2L && geom$chi1_available[[1]] &&
       geom$chi1_available[[2]] && is.finite(geom$omega[[1]]) &&
       is.finite(geom$peptide_bond_length[[1]]),
       "Complete peptide/side chain must yield finite diagnostic angles.")
assert(is.na(geom$omega[[2]]) && geom$omega_status[[2]]=="Missing",
       "Termini must not be assigned peptide omega.")
# Break between chains even if atom coordinates remain chemically close.
other <- torsion
other$bonded_to_next[] <- FALSE
disconnected <- ram_extra_geometry(pdb,other)
assert(all(is.na(disconnected$omega)),
       "Peptide geometry must not cross an unconnected residue boundary.")
# Missing sidechain atoms produce unknown chi1, not a fake rotamer outlier.
deleted <- pdb
deleted$atom <- deleted$atom[deleted$atom$elety!="OG",,drop=FALSE]
missing <- ram_extra_geometry(deleted,torsion)
assert(is.na(missing$chi1[[1]]) && !missing$chi1_available[[1]],
       "Missing side-chain atoms must remain explicitly unclassified.")
# Insertion-code and chain identifiers are preserved by the backbone join.
joined <- ram_join_geometry(torsion,geom)
assert(nrow(joined)==nrow(torsion) &&
         identical(joined$chain,torsion$chain),
       "Geometry joins must preserve backbone row order.")
assert(inherits(try(ram_join_geometry(torsion,geom[1,,drop=FALSE]),
                    silent=TRUE),"try-error"),
       "Misaligned geometry tables must be rejected.")
message("Descriptive peptide/side-chain geometry checks passed.")
