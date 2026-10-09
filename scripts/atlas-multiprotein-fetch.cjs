#!/usr/bin/env node
"use strict";
// Reproducible, opt-in multi-fold live-archive benchmark dataset builder.
// Downloads PDBe *updated* mmCIF files and uses exactly the browser's SIFTS
// and atom-site parsers. No inferred numbering or silent fallbacks.
const https = require("node:https");
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");
const sifts = require("../shinyRam/www/atlas-sifts-exact.js");

const config = JSON.parse(fs.readFileSync(
  "benchmarks/atlas-multiprotein-cases.json","utf8"));
if(config.schema_version!==1 || !Array.isArray(config.proteins) ||
   config.proteins.length<3)
  throw new Error("Curated multi-protein benchmark configuration missing.");
const entries = [];
for(const protein of config.proteins) {
  if(!/^[A-Z0-9]{6,12}(?:-[1-9][0-9]*)?$/.test(protein.uniprot))
    throw new Error("Invalid UniProt accession in case manifest.");
  for(const item of protein.entries) {
    if(!/^[A-Z0-9]{4}$/.test(item.pdb) ||
       !Number.isInteger(item.entity) || item.entity<1)
      throw new Error("Invalid PDB entity in manifest.");
    entries.push({protein:protein.id,pdb:item.pdb,entity:item.entity,
      accession:protein.uniprot,context:item.label+" / "+item.condition});
  }
}
const ids = entries.map(e=>e.pdb+"_"+e.entity);
if(new Set(ids).size!==ids.length)
  throw new Error("Duplicate benchmark PDB entity IDs.");
const destination=process.argv[2] || "benchmarks/output/atlas-multiprotein";
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
      "Same-crystal chains are not independent biological replicates",
      "Separate crystals are not guaranteed identical functional states",
      "Same UniProt accession does not guarantee matching constructs or mutations",
      "Conformational labels originate from literature, not the algorithm",
      "No universal state classifier has been calibrated"
    ],
    curated_cases:config.proteins.map(p=>({id:p.id,uniprot:p.uniprot,
      fold_family:p.fold_family})),
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
    if(!chains.length)
      throw new Error(entry.pdb+": no observed mapped protein chains.");
    // Explicitly check the residue chemistry across structures. Exact SIFTS
    // mapping alone is not proof of construct/sequence compatibility.
    const scheme=sifts.parseLoops(cif,["_pdbx_poly_seq_scheme"])
      ["_pdbx_poly_seq_scheme"];
    const monByKey=new Map();
    for(const row of scheme.rows) {
      const key=row.asym_id+"|"+row.seq_id;
      if(!monByKey.has(key) && row.mon_id && ![".","?"].includes(row.mon_id))
        monByKey.set(key,row.mon_id);
    }
    const sequence_audit=exact.rows.filter(x=>x.observed).map(x=>({
      struct_asym_id:x.struct_asym_id,
      uniprot_resi:x.uniprot_resi,
      residue_name:monByKey.get(x.struct_asym_id+"|"+x.label_seq_id)||null
    }));
    const name=entry.pdb+"_"+entry.entity;
    const record={state:"mapped",pdb_id:entry.pdb,entity_id:1,
      accession:entry.accession,source_url:url,
      mapping:exact.rows,backbone_atoms,ca_points,sequence_audit};
    fs.writeFileSync(path.join(destination,name+".json"),
      JSON.stringify(record));
    const audit={id:name,protein:entry.protein,context:entry.context,url,
      sha256:crypto.createHash("sha256").update(bytes).digest("hex"),
      size_bytes:bytes.length,observed_uniprot_positions:uniquePositions.size,
      mapped_sifts_rows:exact.rows.length,
      ca_atoms:ca_points.length,backbone_atoms:backbone_atoms.length,
      polymer_asym_ids:chains,
      residue_names_available:sequence_audit.filter(x=>x.residue_name).length,
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
