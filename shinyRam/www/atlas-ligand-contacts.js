(function(root) {
  "use strict";
  // Observed non-water HETATM proximity to *exactly mapped* target-polymer
  // heavy atoms. Does not claim biochemical ligand binding or apo/holo state.
  const WATER = new Set(["HOH","DOD","WAT","H2O","OH2"]);
  const clean = x => typeof x === "string" && x !== "." && x !== "?" && x.length ?
    x : null;
  const number = x => x !== "." && x !== "?" && Number.isFinite(Number(x)) ?
    Number(x) : null;
  const rowCell = (row,fields,name) => {
    const i=fields.indexOf("_atom_site."+name);
    return i<0 ? null : row[i];
  };
  const xyzCell = (row,headers) => ["cartn_x","cartn_y","cartn_z"]
    .map(name=>number(rowCell(row,headers,name)));
  const binKey = (coords) => coords.map(n=>Math.floor(n/5)).join("|");

  function extract(cif,exact,radius=4.5) {
    if(typeof cif!=="string" || cif.length>45*1024*1024)
      throw Error("mmCIF missing or beyond browser coordinate limit.");
    if(!exact || !Array.isArray(exact.rows) ||
       !(radius>=2 && radius<=8))
      throw Error("Invalid exact SIFTS mapping or contact radius.");
    const mapped=new Map();
    for(const r of exact.rows.filter(r=>r.observed)) {
      const key=r.struct_asym_id+"|"+r.label_seq_id;
      // No unknown or ambiguous SIFTS positions contribute contacts.
      if(mapped.has(key) && mapped.get(key)!==r.uniprot_resi)
        mapped.set(key,null);
      else if(!mapped.has(key)) mapped.set(key,r.uniprot_resi);
    }
    if(!mapped.size) return {status:"unavailable",
      warning:"No exact observed protein positions for ligand contacts.",
      radius_A:radius,sites:[],total_nonwater_sites:0};

    const lines=cif.split(/\r?\n/);
    const protein=[],hetero=new Map();
    let atomLoop=false;
    for(let i=0;i<lines.length;i++) {
      if(lines[i].trim()!=="loop_")continue;
      let headers=[];
      while(i+1<lines.length && lines[i+1].trim().startsWith("_"))
        headers.push(lines[++i].trim().toLowerCase());
      if(!headers.length || !headers[0].startsWith("_atom_site."))continue;
      if(atomLoop)throw Error("Multiple atom_site coordinate loops.");
      atomLoop=true;
      const required=["group_pdb","label_atom_id","label_asym_id",
        "label_seq_id","label_comp_id","type_symbol","cartn_x","cartn_y","cartn_z"];
      if(required.some(name=>!headers.includes("_atom_site."+name)))
        throw Error("Updated mmCIF lacks fields needed for ligand contacts.");
      const buffer=[];
      for(;i+1<lines.length;) {
        const next=lines[i+1].trim();
        if(next==="#"||next==="loop_"||next.startsWith("data_")||
           next.startsWith("save_")||next.startsWith("_"))break;
        i++;
        if(!next)continue;
        if(lines[i].startsWith(";"))
          throw Error("Unsupported multiline atom_site coordinate.");
        buffer.push(...root.RamplotRExactSifts.tokenize(lines[i]));
        while(buffer.length>=headers.length) {
          const row=buffer.splice(0,headers.length);
          const group=rowCell(row,headers,"group_pdb");
          if(group!=="ATOM" && group!=="HETATM")continue;
          const model=rowCell(row,headers,"pdbx_pdb_model_num");
          if(model!==null && model!=="1")continue;
          const alt=rowCell(row,headers,"label_alt_id");
          if(alt!==null && alt!=="." && alt!=="?" && alt!=="A")
            continue;
          const element=String(rowCell(row,headers,"type_symbol")||"").toUpperCase();
          if(!element || ["H","D"].includes(element))continue;
          const xyz=xyzCell(row,headers);
          if(xyz.some(n=>n===null || Math.abs(n)>100000))continue;
          const asym=rowCell(row,headers,"label_asym_id");
          const seq=rowCell(row,headers,"label_seq_id");
          const mappedKey=asym+"|"+seq;
          if(group==="ATOM" && Number.isInteger(mapped.get(mappedKey))) {
            protein.push({xyz,pos:mapped.get(mappedKey)});
            if(protein.length>60000)
              throw Error("Too many mapped protein heavy atoms for local contacts.");
          } else if(group==="HETATM") {
            // HETATM can also encode a modified polymer residue, e.g. MSE.
            // A defined polymer label_seq_id is not a nonpolymer ligand site.
            if(seq!==null && seq!=="." && seq!=="?")continue;
            const comp=clean(rowCell(row,headers,"label_comp_id"));
            if(!comp || WATER.has(comp.toUpperCase()))continue;
            const author=clean(rowCell(row,headers,"auth_seq_id")) || "?";
            const ins=clean(rowCell(row,headers,"pdbx_pdb_ins_code")) || "";
            const key=[asym||"?",author,ins,comp].join("|");
            if(!hetero.has(key))hetero.set(key,{
              comp_id:comp.toUpperCase(),asym_id:asym||"?",
              auth_seq_id:author,insertion_code:ins,coords:[]});
            hetero.get(key).coords.push(xyz);
            if(hetero.size>2000)
              throw Error("More than 2000 non-water hetero sites in model 1.");
          }
        }
      }
      if(buffer.length)throw Error("Incomplete mmCIF atom_site record.");
    }
    if(!atomLoop||!protein.length)return {status:"unavailable",
      warning:"No complete model-1 heavy-atom coordinates for the exact-mapped protein.",
      radius_A:radius,sites:[],total_nonwater_sites:hetero.size};
    const grid=new Map();
    for(const p of protein) {
      const key=binKey(p.xyz);
      if(!grid.has(key))grid.set(key,[]);
      grid.get(key).push(p);
    }
    const close=[];
    for(const site of hetero.values()) {
      let best=Infinity,pos=null;
      const residueDistances=new Map();
      for(const atom of site.coords) {
        const c=atom.map(x=>Math.floor(x/5));
        // Neighboring 5 Å grid cells contain every point within 4.5 Å.
        for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)
          for(let dz=-1;dz<=1;dz++) {
            const neighbors=grid.get([c[0]+dx,c[1]+dy,c[2]+dz].join("|"))||[];
            for(const p of neighbors) {
              const d2=(atom[0]-p.xyz[0])**2+(atom[1]-p.xyz[1])**2+
                (atom[2]-p.xyz[2])**2;
              if(d2<=radius*radius) {
                const old=residueDistances.get(p.pos);
                if(old===undefined || d2<old)
                  residueDistances.set(p.pos,d2);
                if(d2<best) { best=d2;pos=p.pos; }
              }
            }
          }
      }
      if(pos!==null)close.push({comp_id:site.comp_id,
        asym_id:site.asym_id,auth_seq_id:site.auth_seq_id,
        insertion_code:site.insertion_code,
        heavy_atoms:site.coords.length,min_distance_A:Math.sqrt(best),
        nearest_uniprot_resi:pos,
        residue_contacts:[...residueDistances].sort((a,b)=>a[0]-b[0])
          .map(([uniprot_resi,d2])=>({uniprot_resi,
            min_distance_A:Math.sqrt(d2)}))});
    }
    close.sort((a,b)=>a.min_distance_A-b.min_distance_A||
      a.comp_id.localeCompare(b.comp_id));
    // An overlarge output is unavailable rather than falsely complete.
    if(close.length>250)throw Error("Too many nearby hetero sites to report.");
    return{status:"measured",radius_A:radius,
      total_nonwater_sites:hetero.size,nearby_sites:close.length,
      sites:close,protein_heavy_atoms:protein.length,
      warning:hetero.size===0?"No non-water nonpolymer HETATM sites in model 1; this does not establish an apo state.":""};
  }
  const api={extract};
  root.RamplotRAtlasLigands=api;
  if(typeof module!=="undefined"&&module.exports)module.exports=api;
})(typeof window!=="undefined"?window:globalThis);
