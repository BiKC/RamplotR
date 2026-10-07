# Pure canonical-coordinate tests; no network access required.
source(file.path("shinyRam","R","canonical.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

payload <- list(
  "1abc"=list(
    UniProt=list(
      P12345=list(
        identifier="TEST_HUMAN",
        mappings=list(
          list(
            entity_id=1L,chain_id="A",struct_asym_id="A",
            unp_start=10L,unp_end=12L,identity=1,coverage=0.5,
            start=list(residue_number=1L,author_residue_number=101L,
                       author_insertion_code=""),
            end=list(residue_number=3L,author_residue_number=103L,
                     author_insertion_code="")
          ),
          list(
            entity_id=1L,chain_id="A",struct_asym_id="A",
            unp_start=20L,unp_end=21L,identity=1,coverage=0.5,
            start=list(residue_number=10L,author_residue_number=200L,
                       author_insertion_code=""),
            end=list(residue_number=11L,author_residue_number=201L,
                     author_insertion_code="A")
          )
        )
      )
    )
  )
)

segments <- ram_sifts_payload_segments(payload,"1ABC")
assert(nrow(segments)==2L &&
       identical(as.character(segments$uniprot_accession),
                 c("P12345","P12345")),
       "SIFTS payload parsing lost accession mappings.")
assert(isTRUE(segments$safe_linear[[1L]]) &&
       !isTRUE(segments$safe_linear[[2L]]),
       "Only provably linear author-number ranges may be expanded.")

mapping <- ram_sifts_expand_safe(segments)
assert(nrow(mapping)==3L &&
       identical(mapping$resi,c(101L,102L,103L)) &&
       identical(mapping$uniprot_resi,c(10L,11L,12L)),
       "Safe SIFTS range expansion changed.")
assert(!any(mapping$resi %in% c(200L,201L)),
       "Insertion-code range must not be fabricated into canonical residues.")

residues <- data.frame(
  chain="A",resi=c(101L,102L,103L,200L),insertion_code="",
  resn=c("ALA","SER","GLY","THR"),stringsAsFactors=FALSE
)
joined <- ram_canonical_join(residues,mapping)
assert(identical(as.character(joined$canonical_status),
                 c("mapped","mapped","mapped","unmapped")) &&
       identical(joined$uniprot_resi[1:3],c(10L,11L,12L)),
       "Canonical residue join changed.")

conflict <- rbind(mapping[1,,drop=FALSE],mapping[1,,drop=FALSE])
conflict$uniprot_accession[[2L]] <- "Q99999"
ambiguous <- ram_canonical_join(residues[1,,drop=FALSE],conflict)
assert(ambiguous$canonical_status[[1L]]=="ambiguous" &&
       is.na(ambiguous$uniprot_resi[[1L]]),
       "Conflicting canonical targets must remain explicitly ambiguous.")

af <- data.frame(
  chain="A",resi=1:3,insertion_code="",resn=c("ALA","GLY","SER"),
  stringsAsFactors=FALSE
)
af_map <- ram_afdb_canonical_map(af,"p12345")
assert(nrow(af_map)==3L &&
       all(af_map$uniprot_accession=="P12345") &&
       identical(af_map$uniprot_resi,1:3),
       "AlphaFold DB canonical numbering changed.")
assert(inherits(try(ram_afdb_canonical_map(
  rbind(af,transform(af,chain="B")),"P12345"),silent=TRUE),"try-error"),
  "AlphaFold DB canonical mapping must reject unexpected multichain input.")

bad <- segments[1,,drop=FALSE]
bad$author_end <- 105L
bad <- ram_sifts_mark_safe_segments(bad)
assert(!bad$safe_linear[[1L]],
       "Unequal PDB/UniProt spans must not be treated as one-to-one mappings.")

message("Canonical SIFTS/UniProt mapping tests passed.")
