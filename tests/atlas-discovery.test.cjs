const assert = require("node:assert/strict");
const atlas = require("../shinyRam/www/atlas-discovery.js");

assert.equal(atlas.validateAccession("p00533"), "P00533");
assert.throws(() => atlas.validateAccession("../bad"), /UniProt/);
assert.throws(() => atlas.validateAccession(""), /UniProt/);

const q = atlas.searchRequest("P00533", 500);
assert.equal(q.query.parameters.attribute,
  "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession");
assert.equal(q.query.parameters.operator, "exact_match");
assert.equal(q.query.parameters.value, "P00533");
assert.equal(q.return_type, "polymer_entity");
assert.deepEqual(q.request_options.results_content_type, ["experimental"]);
assert.deepEqual(q.request_options.paginate, {start: 0, rows: 50});

assert.deepEqual(atlas.parseEntityId("1abc_2"), {pdb_id:"1ABC",entity_id:"2"});
assert.equal(atlas.parseEntityId("bad/name"), null);
const hit = {identifier:"1ABC_2"};
const entity = {
  rcsb_polymer_entity_container_identifiers: {auth_asym_ids:["A","B","A"]},
  rcsb_polymer_entity: {pdbx_description:"Test protein"}
};
const entry = {
  exptl:[{method:"X-RAY DIFFRACTION"}],
  rcsb_entry_info:{resolution_combined:[2.2,1.9]},
  rcsb_accession_info:{initial_release_date:"2025-04-01"},
  struct:{title:"My experimental structure"}
};
assert.deepEqual(atlas.normalizeEntity(hit,entity,entry), {
  pdb_id:"1ABC",entity_id:"2",chains:["A","B"],chain:"A",
  description:"Test protein",method:"X-RAY DIFFRACTION",
  resolution:1.9,release_date:"2025-04-01",
  title:"My experimental structure"
});
assert.equal(atlas.normalizeEntity(hit,null,entry),null);
assert.equal(atlas.normalizeEntity({identifier:"not-an-id"},entity,entry),null);

(async () => {
  let active = 0, peak = 0;
  const list = await atlas.runBounded(Array.from({length:20},(_,i)=>i),
    async i => {
      active++;
      peak=Math.max(peak,active);
      await Promise.resolve();
      active--;
      return i*2;
    });
  assert.ok(peak<=6,"Do not launch unbounded metadata requests.");
  assert.deepEqual(list,Array.from({length:20},(_,i)=>i*2));
  console.log("Atlas discovery input, metadata and request bounds passed.");
})().catch(e=>{console.error(e);process.exitCode=1;});
