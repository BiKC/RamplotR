#!/usr/bin/env node
"use strict";
// Reproducible, opt-in real-archive Atlas benchmark dataset builder.
// Downloads PDBe *updated* mmCIF files and uses exactly the browser's SIFTS
// and atom-site parsers. No inferred numbering or silent fallbacks.
const https = require("node:https");
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");
const sifts = require("../shinyRam/www/atlas-sifts-exact.js");

const entries = [
  {pdb:"4AKE",entity:"1",accession:"P69441",context:"apo/open"},
  {pdb:"1AKE",entity:"1",accession:"P69441",context:"inhibitor-bound/closed"}
];
const destination = process.argv[2] || "benchmarks/output/atlas-adk";
fs.mkdirSync(destination,{recursive:true});

function fetchText(url,redirects=0) {
  return new Promise((resolve,reject) => {
    if(redirects>3) return reject(new Error("Too many redirects."));
    https.get(url,{headers:{"User-Agent":"RamplotR-research-benchmark/1.0"},
      timeout:45000},response=>{
      if([301,302,303,307,308].includes(response.statusCode)) {
        const next=response.headers.location;
        response.resume();
        if(!next) return reject(new Error("Empty redirect."));
        return resolve(fetchText(new URL(next,url).href,redirects+1));
      }
      if(response.statusCode!==200) {
        response.resume();
        return reject(new Error("PDBe returned HTTP "+response.statusCode));
      }
      const chunks=[];
      let size=0;
      response.on("data",chunk=>{
        size+=chunk.length;
        if(size>45*1024*1024) response.destroy(
          new Error("Updated mmCIF file exceeds 45 MB safety limit."));
        else chunks.push(chunk);
      });
      response.on("end",()=>resolve(Buffer.concat(chunks)));
      response.on("error",reject);
    }).on("timeout",function(){this.destroy(new Error("PDBe request timeout."));})
      .on("error",reject);
  });
}

(async()=>{
  const manifest={
    generated_at:new Date().toISOString(),
    coordinate_source:"PDBe updated mmCIF, model 1",
    mapping_source:"_pdbx_sifts_xref_db joined to _pdbx_poly_seq_scheme",
    protocol:"Direct SIFTS label_asym_id + label_seq_id only",
    known_limitations:[
      "Chains in one crystal are not independent biological replicates",
      "Construct and isoform equivalence must be independently checked",
      "No functional-state classifier is trained or calibrated"
    ],
    entries:[]
  };
  for(const entry of entries) {
    const url="https://www.ebi.ac.uk/pdbe/static/entry/"+
      entry.pdb.toLowerCase()+"_updated.cif";
    const bytes=await fetchText(url);
    const cif=bytes.toString("utf8");
    const exact=sifts.exactRows(cif,entry.pdb,entry.entity,entry.accession);
    const backbone_atoms=sifts.extractMappedBackbone(cif,exact);
    const ca_points=backbone_atoms.filter(x=>x.atom_name==="CA")
      .map(({atom_name,...point})=>point);
    const uniquePositions=new Set(exact.rows.filter(x=>x.observed)
      .map(x=>x.uniprot_resi));
    if(uniquePositions.size<100 || ca_points.length<100)
      throw new Error(entry.pdb+
        ": <100 observed canonical positions/CA atoms; do not infer mapping.");
    const chains=[...new Set(exact.rows.filter(x=>x.observed)
      .map(x=>x.struct_asym_id))];
    if(chains.length<2)
      throw new Error(entry.pdb+
        ": expected two observed chains for within-crystal controls.");
    const name=entry.pdb+"_1";
    const record={state:"mapped",pdb_id:entry.pdb,entity_id:1,
      accession:entry.accession,source_url:url,
      mapping:exact.rows,backbone_atoms,ca_points};
    fs.writeFileSync(path.join(destination,name+".json"),
      JSON.stringify(record));
    const audit={id:name,context:entry.context,url,
      sha256:crypto.createHash("sha256").update(bytes).digest("hex"),
      size_bytes:bytes.length,observed_uniprot_positions:uniquePositions.size,
      mapped_sifts_rows:exact.rows.length,
      ca_atoms:ca_points.length,backbone_atoms:backbone_atoms.length,
      polymer_asym_ids:chains,
      unmatched_sifts_rows:exact.unlinked_sifts_rows};
    manifest.entries.push(audit);
    console.log(JSON.stringify(audit));
  }
  fs.writeFileSync(path.join(destination,"manifest.json"),
    JSON.stringify(manifest,null,2));
})().catch(error=>{
  console.error("Live PDBe benchmark retrieval failed:",error.stack || error);
  process.exitCode=1;
});
