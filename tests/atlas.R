# Pure Atlas pagination and cohort integrity regressions.
source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas.R"))
assert <- function(x,message) if(!isTRUE(x)) stop(message,call.=FALSE)
rec <- function(pdb,entity) list(pdb_id=pdb,entity_id=as.character(entity),
  description="Test experimental entity")
base_records <- lapply(seq_len(48L),function(i) rec("1ABC",i))
first <- ram_atlas_merge_page(NULL,list(
  accession="P12345",start=0L,total_count=75L,
  returned_count=50L,incomplete_metadata=2L,
  failed_entity_ids=c("1ABC_49","1ABC_50"),results=base_records))
assert(first$next_offset==50L && first$total_count==75L &&
       first$enriched_count==48L && first$incomplete_metadata==2L &&
       first$pages==1L && isTRUE(first$has_more),
       "First page must retain raw offset and partial metadata.")
second_records <- c(list(rec("1ABC",1L)),
                    lapply(seq_len(23L),function(i) rec("2XYZ",i)))
second <- ram_atlas_merge_page(first,list(
  accession="P12345",start=50L,total_count=75L,
  returned_count=25L,incomplete_metadata=1L,
  failed_entity_ids="2XYZ_24",results=second_records))
assert(identical(second$failed_entity_ids,c("1ABC_49","1ABC_50","2XYZ_24")),
       "Failed metadata entity IDs should remain auditable across pages.")
assert(second$next_offset==75L && second$enriched_count==71L &&
       second$incomplete_metadata==3L && second$duplicate_count==1L &&
       second$pages==2L && !isTRUE(second$has_more),
       "Second page must deduplicate without altering RCSB hit offsets.")
assert(identical(second$results[[1L]]$pdb_id,"1ABC"),
       "Earlier entity metadata must win when duplicate keys recur.")
should_fail <- function(page,msg) {
  assert(inherits(try(ram_atlas_merge_page(first,page),silent=TRUE),
                  "try-error"),msg)
}
second_page <- list(accession="P12345",start=50L,total_count=75L,
  returned_count=25L,incomplete_metadata=1L,results=second_records)
should_fail(within(second_page,start<-49L),
            "Out-of-order pages must fail.")
should_fail(within(second_page,total_count<-76L),
            "An archive change must require query restart.")
should_fail(within(second_page,accession<-"Q99999"),
            "Pages from another UniProt protein must fail.")
should_fail(within(second_page,incomplete_metadata<-26L),
            "Inconsistent missing metadata must fail.")
should_fail(within(second_page,returned_count<-1L),
            "Entity records may not exceed returned hit count.")
empty <- ram_atlas_merge_page(NULL,list(
  accession="P12345",start=0L,total_count=0L,
  returned_count=0L,incomplete_metadata=0L,results=list()))
assert(!empty$has_more && empty$enriched_count==0L,
       "Empty searches must terminate.")
stalled <- ram_atlas_merge_page(NULL,list(
  accession="P12345",start=0L,total_count=3L,
  returned_count=0L,incomplete_metadata=0L,results=list()))
assert(stalled$stalled && !stalled$has_more,
       "Missing pages must not allow infinite next-page loops.")
message("Atlas cohort pagination integrity tests passed.")
