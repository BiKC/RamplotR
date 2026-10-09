const assert=require("node:assert/strict");
const fs=require("node:fs");
const path=require("node:path");
const conf=JSON.parse(fs.readFileSync(
  path.join(__dirname,"../benchmarks/atlas-multiprotein-cases.json"),
  "utf8"));
assert.equal(conf.schema_version,1);
assert.ok(conf.proteins.length>=4);
const proteinIds=new Set();
const foldFamilies=new Set();
const fullIds=new Set();
for(const protein of conf.proteins) {
  assert.match(protein.uniprot,/^[A-Z0-9]{6,12}(?:-[1-9]\d*)?$/);
  assert.ok(!proteinIds.has(protein.id),"Duplicate protein benchmark case.");
  proteinIds.add(protein.id);
  assert.match(protein.fold_family,/^[a-z][a-z0-9_]+$/);
  foldFamilies.add(protein.fold_family);
  const entries=new Set();
  for(const e of protein.entries) {
    assert.match(e.pdb,/^[A-Z0-9]{4}$/);
    assert.ok(Number.isInteger(e.entity)&&e.entity>0);
    assert.ok(!entries.has(e.pdb),"Duplicated PDB in benchmark case.");
    entries.add(e.pdb);
    const key=e.pdb+"_"+e.entity;
    assert.ok(!fullIds.has(key),"Unexpected shared entity across benchmark proteins.");
    fullIds.add(key);
    assert.ok(e.label==="open"||e.label==="closed");
  }
  let positive=0;
  let control=0;
  for(const comparison of protein.comparisons) {
    assert.ok(entries.has(comparison.left)&&entries.has(comparison.right));
    if(comparison.class==="documented_open_closed") {
      positive++;
      assert.notEqual(comparison.left,comparison.right);
      const byId=Object.fromEntries(protein.entries.map(x=>[x.pdb,x.label]));
      assert.equal(byId[comparison.left],"open");
      assert.equal(byId[comparison.right],"closed");
    } else if(comparison.class==="within_crystal") {
      control++;
      assert.equal(comparison.left,comparison.right);
      assert.notEqual(comparison.left_chain_rank,comparison.right_chain_rank);
      assert.ok(comparison.left_chain_rank>0 && comparison.right_chain_rank>0);
    } else if(comparison.class==="separate_crystal_open_control") {
      control++;
      assert.notEqual(comparison.left,comparison.right);
      const byId=Object.fromEntries(protein.entries.map(x=>[x.pdb,x.label]));
      assert.equal(byId[comparison.left],"open");
      assert.equal(byId[comparison.right],"open");
    } else {
      throw new Error("Unsupported benchmark control or positive type.");
    }
  }
  assert.equal(positive,1,"Exactly one documented contrast per protein.");
  assert.ok(control>=1,"Explicit negative/same-state control required.");
}
assert.ok(foldFamilies.size>=3,
  "Benchmark must span at least three distinct structural fold families.");
const citrate=conf.proteins.find(x=>x.id==="pig_citrate_synthase");
assert.deepEqual(citrate.entries.map(e=>e.pdb),["1CTS","2CTS","3ENJ"]);
assert.equal(citrate.uniprot,"P00889");
assert.ok(citrate.comparisons.some(x=>x.left==="1CTS" && x.right==="3ENJ" &&
  x.class==="separate_crystal_open_control"));
console.log("Curated Atlas benchmark manifest checks passed ("+
  conf.proteins.length+" proteins, "+fullIds.size+" PDB entities).");
