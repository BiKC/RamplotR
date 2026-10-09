# Exact-residue Atlas SIFTS tests; source location is repository root.
source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas-sifts.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
row <- function(resi,ins,unp,observed=TRUE,mon="ALA") {
  list(chain="A",resi=resi,insertion_code=ins,
       uniprot_accession="P12345",uniprot_resi=unp,
       entity_id=1L,struct_asym_id="X",observed=observed,mon_id=mon)
}
payload <- list(pdb_id="1ABC",entity_id="1",accession="P12345",
  rows=list(row(101L,"",20L,mon="MET"),row(101L,"A",21L,mon="MSE"),
            row(103L,"",25L,FALSE)))
mapping <- ram_atlas_exact_sifts_map(payload,"1abc",1L,"p12345")
assert(identical(mapping$mon_id,c("MET","MSE","ALA")),
  "Experimental sequence scheme must retain modified-residue chemistry.")
assert(nrow(mapping)==3L &&
       identical(mapping$uniprot_resi,c(20L,21L,25L)) &&
       identical(mapping$insertion_code,c("","A","")),
       "Exact SIFTS must preserve non-linear numbering and insertion codes.")
info <- ram_atlas_sifts_summary(mapping)
assert(info$unique_residues==3L && info$observed_residues==2L &&
       info$conflicting_residues==0L,"Exact mapping coverage wrong.")
torsions <- data.frame(chain="A",resi=c(101L,101L,103L,102L),
  insertion_code=c("","A","",""),stringsAsFactors=FALSE)
joined <- ram_canonical_join(torsions,mapping)
assert(identical(joined$canonical_status,c("mapped","mapped","mapped","unmapped")) &&
       identical(joined$uniprot_resi[1:3],c(20L,21L,25L)),
       "Canonical join must use exact chain/residue/insertion identifiers.")
duplicate <- payload
duplicate$rows[[4L]] <- row(101L,"",22L)
different <- ram_atlas_exact_sifts_map(duplicate,"1ABC",1L,"P12345")
assert(ram_atlas_sifts_summary(different)$conflicting_residues==1L,
       "Conflicts must be counted.")
assert(ram_canonical_join(torsions,different)$canonical_status[[1L]]=="ambiguous",
       "Conflicting UniProt residue numbers must remain ambiguous.")
wrong <- payload;wrong$entity_id <- "2"
assert(inherits(try(ram_atlas_exact_sifts_map(wrong,"1ABC",1L,"P12345"),
                    silent=TRUE),"try-error"),"Wrong entity accepted.")
wrong <- payload;wrong$rows[[1L]]$uniprot_accession <- "Q99999"
assert(inherits(try(ram_atlas_exact_sifts_map(wrong,"1ABC",1L,"P12345"),
                    silent=TRUE),"try-error"),"Wrong accession accepted.")
message("Exact PDBe SIFTS residue mapping tests passed.")
