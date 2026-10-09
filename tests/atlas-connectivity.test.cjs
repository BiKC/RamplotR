"use strict";
const assert=require("node:assert/strict");
require("../shinyRam/www/atlas-discovery.js");
require("../shinyRam/www/atlas-sifts-exact.js");
const api=require("../shinyRam/www/atlas-connectivity.js");

const syntheticCif=[
 "data_test","#","loop_",
 "_pdbx_poly_seq_scheme.asym_id",
 "_pdbx_poly_seq_scheme.entity_id",
 "_pdbx_poly_seq_scheme.seq_id",
 "_pdbx_poly_seq_scheme.pdb_strand_id",
 "_pdbx_poly_seq_scheme.pdb_seq_num",
 "_pdbx_poly_seq_scheme.pdb_ins_code",
 "_pdbx_poly_seq_scheme.mon_id",
 "X 1 1 A 101 . MET","X 1 2 A 102 . GLY","#","loop_",
 "_pdbx_sifts_xref_db.entity_id",
 "_pdbx_sifts_xref_db.asym_id",
 "_pdbx_sifts_xref_db.seq_id",
 "_pdbx_sifts_xref_db.unp_acc",
 "_pdbx_sifts_xref_db.unp_num",
 "_pdbx_sifts_xref_db.observed",
 "1 X 1 P69441 1 1","1 X 2 P69441 2 1","#","loop_",
 "_atom_site.group_PDB",
 "_atom_site.label_atom_id",
 "_atom_site.label_asym_id",
 "_atom_site.label_seq_id",
 "_atom_site.Cartn_x",
 "_atom_site.Cartn_y",
 "_atom_site.Cartn_z",
 "_atom_site.pdbx_PDB_model_num",
 "_atom_site.label_alt_id",
 "ATOM CA X 1 1 2 3 1 .","ATOM CA X 2 4 5 6 1 .","#"
].join("\n");
function response(status,body,type="json") {
  return {status,ok:status>=200&&status<300,
    json:async()=>type==="json"?body:JSON.parse(body),
    text:async()=>type==="text"?body:JSON.stringify(body)};
}
const calls=[];
const successful=async (url,options)=>{
  calls.push({url,options});
  if(url.includes("rcsbsearch"))
    return response(200,{total_count:12,result_set:[{identifier:"4AKE_1"}]});
  if(url.includes("polymer_entity"))
    return response(200,{rcsb_polymer_entity_container_identifiers:
      {entry_id:"4AKE",entity_id:"1"}});
  if(url.includes("4ake_updated.cif"))
    return response(200,syntheticCif,"text");
  throw Error("Unknown endpoint "+url);
};

(async()=>{
  assert.equal(api.descriptors.length,3);
  const result=await api.runChecks({fetcher:successful,timeoutMs:1000});
  assert.equal(result.state,"ok");
  assert.equal(result.passed,3);
  assert.equal(result.total,3);
  assert.equal(calls.length,3);
  assert.deepEqual(result.checks.map(x=>x.state),["ok","ok","ok"]);
  assert.match(result.checks[2].summary,/2 observed mapped residue rows/);
  assert.ok(calls.every(c=>c.options.credentials==="omit" &&
    c.options.cache==="no-store"));
  const submitted=JSON.parse(calls[0].options.body);
  assert.equal(submitted.query.parameters.value,"P69441");
  assert.deepEqual(submitted.request_options.paginate,{start:0,rows:1});
  assert.equal(calls[0].options.method,"POST");
  assert.equal(calls[1].options.method,undefined);

  const progress=[];
  const broken=await api.runChecks({
    fetcher:async (url)=>{
      if(url.includes("rcsbsearch"))
        throw new TypeError("Failed to fetch");
      if(url.includes("polymer_entity"))
        return response(503,{error:"maintenance"});
      return response(200,"data_broken\n#","text");
    },
    onProgress:checks=>progress.push(checks.map(c=>c.state))
  });
  assert.equal(broken.passed,0);
  assert.equal(broken.state,"partial_failure");
  assert.deepEqual(broken.checks.map(x=>x.reason),
    ["network_or_cors","invalid_response","invalid_response"]);
  assert.equal(broken.checks[1].http_status,503);
  assert.match(broken.checks[0].detail,/CORS/);
  assert.deepEqual(progress.map(x=>x.length),[1,2,3]);
  assert.deepEqual(progress[2],["failed","failed","failed"]);

  const mixed=await api.runChecks({fetcher:async (url,options)=>
    url.includes("polymer_entity")?
      response(200,{unrelated:"payload"}):
      successful(url,options)});
  assert.equal(mixed.passed,2);
  assert.equal(mixed.checks[1].reason,"invalid_response");

  const abort=new Error("Aborted");
  abort.name="AbortError";
  const timeout=await api.runChecks({fetcher:async (url,options)=>{
    if(url.includes("rcsbsearch"))throw abort;
    return successful(url,options);
  }});
  assert.equal(timeout.checks[0].reason,"timeout");
  assert.equal(timeout.passed,2);
  console.log("Atlas browser-origin diagnostics passed success, HTTP, CORS-like network errors, missing SIFTS and timeout.");
})().catch(e=>{console.error(e);process.exitCode=1;});
