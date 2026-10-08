(function (root) {
  "use strict";

  // Targeted PDBx/mmCIF reader. Only the two sequence-mapping loops are
  // tokenized; large atom_site loops are skipped. No numbering is interpolated.
  const BASE = "https://www.ebi.ac.uk/pdbe/static/entry/";

  function tokenize(line) {
    const out = [];
    let i = 0;
    while (i < line.length) {
      while (i < line.length && /\s/.test(line[i])) i++;
      if (i >= line.length || line[i] === "#") break;
      let word = "";
      const quote = line[i] === "'" || line[i] === '"' ? line[i++] : null;
      if (quote) {
        let closed = false;
        while (i < line.length) {
          if (line[i] === quote && (i + 1 === line.length ||
              /\s|#/.test(line[i + 1]))) {
            i++; closed = true; break;
          }
          word += line[i++];
        }
        if (!closed) throw new Error("Unterminated quoted mmCIF token.");
      } else {
        while (i < line.length && !/\s/.test(line[i]))
          word += line[i++];
      }
      out.push(word);
    }
    return out;
  }

  function parseLoops(text, categories) {
    if (typeof text !== "string") throw new Error("Expected mmCIF text.");
    const lines = text.split(/\r?\n/);
    const requested = new Set(categories);
    const output = Object.fromEntries(categories.map(category =>
      [category, {headers: [], rows: []}]));
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].trim() !== "loop_") continue;
      const headers = [];
      while (i + 1 < lines.length &&
             lines[i + 1].trim().startsWith("_")) {
        i++;
        const words = tokenize(lines[i].trim());
        if (words.length !== 1 || !words[0].startsWith("_"))
          throw new Error("Malformed mmCIF loop header.");
        headers.push(words[0].toLowerCase());
      }
      if (!headers.length) continue;
      const category = headers[0].split(".")[0];
      const relevant = requested.has(category);
      if (relevant && !headers.every(h => h.startsWith(category + ".")))
        throw new Error("Mixed mmCIF category in mapping loop.");
      const buffer = [];
      while (i + 1 < lines.length) {
        const trimmed = lines[i + 1].trim();
        if (trimmed === "#" || trimmed === "loop_" ||
            trimmed.startsWith("data_") || trimmed.startsWith("save_") ||
            trimmed.startsWith("_")) break;
        i++;
        if (!relevant || !trimmed) continue;
        if (lines[i].startsWith(";"))
          throw new Error("Unexpected multiline value in SIFTS mapping loop.");
        buffer.push(...tokenize(lines[i]));
        if (buffer.length >= headers.length) {
          while (buffer.length >= headers.length) {
            const values = buffer.splice(0, headers.length);
            const row = {};
            for (let j = 0; j < headers.length; j++)
              row[headers[j].slice(category.length + 1)] = values[j];
            output[category].rows.push(row);
            if (output[category].rows.length > 250000)
              throw new Error("SIFTS mapping loop exceeds safety limit.");
          }
        }
      }
      if (relevant) {
        if (buffer.length) throw new Error("Incomplete mmCIF mapping row.");
        output[category].headers = headers;
      }
    }
    return output;
  }

  const present = value => value !== undefined && value !== null &&
    value !== "." && value !== "?" && String(value).trim() !== "";
  const integer = value => present(value) && /^-?\d+$/.test(String(value))
    ? Number(value) : null;
  function validate(pdbId, entityId, accession) {
    const pdb = String(pdbId || "").trim().toLowerCase();
    const entity = String(entityId || "").trim();
    const uniprot = String(accession || "").trim().toUpperCase();
    if (!/^[a-z0-9]{4}$/.test(pdb) ||
        !/^[1-9]\d*$/.test(entity) ||
        !/^[A-Z0-9]{6,12}(?:-[1-9]\d*)?$/.test(uniprot))
      throw new Error("Invalid SIFTS PDB/entity/UniProt identifier.");
    return {pdb, entity, uniprot};
  }

  function exactRows(cif, pdbId, entityId, accession) {
    const requested = validate(pdbId, entityId, accession);
    const categories = ["_pdbx_sifts_xref_db", "_pdbx_poly_seq_scheme"];
    const parsed = parseLoops(cif, categories);
    const sifts = parsed[categories[0]];
    const scheme = parsed[categories[1]];
    const requiredSifts = ["asym_id", "entity_id", "seq_id",
                            "unp_acc", "unp_num"];
    const requiredScheme = ["asym_id", "seq_id", "pdb_strand_id",
                            "pdb_seq_num", "pdb_ins_code"];
    for (const column of requiredSifts)
      if (!sifts.headers.includes(categories[0] + "." + column))
        throw new Error("Updated mmCIF lacks SIFTS field " + column + ".");
    for (const column of requiredScheme)
      if (!scheme.headers.includes(categories[1] + "." + column))
        throw new Error("Updated mmCIF lacks sequence scheme field " + column + ".");

    const bySequence = new Map();
    for (const row of scheme.rows) {
      const seq = integer(row.seq_id);
      if (seq === null || !present(row.asym_id)) continue;
      const key = row.asym_id + "|" + seq;
      if (!bySequence.has(key)) bySequence.set(key, []);
      bySequence.get(key).push(row);
    }
    const mapped = [];
    let matched = 0, unlinked = 0, unobserved = 0;
    for (const row of sifts.rows) {
      if (row.entity_id !== requested.entity ||
          String(row.unp_acc).toUpperCase() !== requested.uniprot) continue;
      const seq = integer(row.seq_id), unp = integer(row.unp_num);
      if (seq === null || unp === null || unp < 1 ||
          !present(row.asym_id)) continue;
      matched++;
      const schemes = bySequence.get(row.asym_id + "|" + seq) || [];
      if (!schemes.length) { unlinked++; continue; }
      for (const s of schemes) {
        const resi = integer(s.pdb_seq_num);
        if (resi === null || !present(s.pdb_strand_id) ||
            (s.pdb_ins_code === "?" || s.pdb_ins_code === undefined)) {
          unlinked++; continue;
        }
        const observed = ["1", "Y", "YES"].includes(
          String(row.observed || "").toUpperCase());
        if (!observed) unobserved++;
        mapped.push({
          chain: String(s.pdb_strand_id), resi,
          insertion_code: present(s.pdb_ins_code) ? String(s.pdb_ins_code) : "",
          uniprot_accession: requested.uniprot, uniprot_resi: unp,
          entity_id: Number(requested.entity), struct_asym_id: row.asym_id,
          label_seq_id: seq, observed
        });
      }
    }
    const unique = new Map();
    for (const row of mapped) {
      const key = [row.chain, row.resi, row.insertion_code,
        row.uniprot_accession, row.uniprot_resi, row.struct_asym_id].join("|");
      if (!unique.has(key)) unique.set(key, row);
    }
    return {
      pdb_id: requested.pdb.toUpperCase(), entity_id: requested.entity,
      accession: requested.uniprot,
      source: "PDBe updated mmCIF _pdbx_sifts_xref_db + _pdbx_poly_seq_scheme",
      matched_sifts_rows: matched, unlinked_sifts_rows: unlinked,
      unobserved_rows: unobserved, rows: [...unique.values()]
    };
  }

  // Read C-alpha coordinates only. Never send entire mmCIF atom tables to
  // Shiny. Require exact (label_asym_id,label_seq_id) references from SIFTS;
  // PDB author numbers are not used for coordinate indexing.
  function extractMappedCA(cif, exact) {
    if (typeof cif !== "string") throw new Error("Expected mmCIF text.");
    if (cif.length > 45 * 1024 * 1024)
      throw new Error("Coordinate table exceeds the 45 MB browser limit.");
    const wanted = new Set(exact.rows.filter(r=>r.observed).map(r=>
      r.struct_asym_id + "|" + r.label_seq_id));
    if (!wanted.size) return [];
    const lines=cif.split(/\r?\n/);
    const output=new Map();
    for(let i=0;i<lines.length;i++) {
      if(lines[i].trim()!=="loop_") continue;
      const headers=[];
      while(i+1<lines.length && lines[i+1].trim().startsWith("_")) {
        const name=lines[++i].trim().toLowerCase();
        headers.push(name);
      }
      if(!headers.length || !headers[0].startsWith("_atom_site.")) continue;
      const fields=["group_pdb","label_atom_id","label_asym_id",
        "label_seq_id","cartn_x","cartn_y","cartn_z"];
      const indices=fields.map(x=>headers.indexOf("_atom_site."+x));
      if(indices.some(x=>x<0))
        throw new Error("Updated mmCIF lacks atom_site C-alpha coordinate fields.");
      const model=headers.indexOf("_atom_site.pdbx_pdb_model_num");
      const alt=headers.indexOf("_atom_site.label_alt_id");
      const buffer=[];
      for(;i+1<lines.length;) {
        const next=lines[i+1].trim();
        if(next==="#" || next==="loop_" || next.startsWith("data_") ||
           next.startsWith("save_") || next.startsWith("_")) break;
        i++;
        if(!next) continue;
        if(lines[i].startsWith(";"))
          throw new Error("Unsupported multiline atom-site record.");
        // Most CIF atom rows occupy one physical line. Retain wrapped rows.
        buffer.push(...tokenize(lines[i]));
        while(buffer.length>=headers.length) {
          const row=buffer.splice(0,headers.length);
          if(row[indices[0]]!=="ATOM" || row[indices[1]]!=="CA") continue;
          if(model>=0 && row[model]!=="1") continue;
          const seq=integer(row[indices[3]]);
          if(seq===null) continue;
          const asym=row[indices[2]], key=asym+"|"+seq;
          if(!wanted.has(key)) continue;
          const alternative=alt<0 ? "." : row[alt];
          if(![".","?","A"].includes(alternative)) continue;
          const xyz=indices.slice(4).map(index=>Number(row[index]));
          if(!xyz.every(n=>Number.isFinite(n)&&Math.abs(n)<100000))
            continue;
          const existing=output.get(key);
          // Prefer the primary unlabelled conformer to alternate A.
          const preferred=alternative==="." || alternative==="?";
          if(!existing || (preferred && !existing.preferred))
            output.set(key,{struct_asym_id:asym,label_seq_id:seq,
              x:xyz[0],y:xyz[1],z:xyz[2],preferred});
          else if(existing.preferred===preferred &&
                  xyz.some((n,j)=>Math.abs(n-[existing.x,existing.y,existing.z][j])>0.001))
            output.set(key,{ambiguous:true});
          if(output.size>5000)
            throw new Error("More than 5000 mapped C-alpha positions.");
        }
      }
      if(buffer.length)
        throw new Error("Incomplete atom_site coordinate record.");
      break;
    }
    return [...output.values()].filter(x=>!x.ambiguous).map(x=>({
      struct_asym_id:x.struct_asym_id,label_seq_id:x.label_seq_id,
      x:x.x,y:x.y,z:x.z
    }));
  }

  async function fetchExact(payload) {
    const requestId = String(payload && payload.request_id || "");
    const id = validate(payload && payload.pdb_id,
      payload && payload.entity_id, payload && payload.accession);
    const url = BASE + id.pdb + "_updated.cif";
    try {
      const controller = new AbortController();
      const timer = setTimeout(() => controller.abort(), 45000);
      let response;
      try {
        response = await fetch(url, {credentials: "omit",
          signal: controller.signal});
        if (!response.ok) throw new Error("PDBe mmCIF HTTP " + response.status);
        const cif = await response.text();
        const result = exactRows(cif, id.pdb, id.entity, id.uniprot);
        if (result.rows.length > 50000)
          throw new Error("SIFTS mapping response is too large.");
        // A coordinate-reader problem does not invalidate exact SIFTS mapping.
        try { result.ca_points = extractMappedCA(cif,result); }
        catch (err) {
          result.ca_points = [];
          result.ca_warning = err && err.message ? err.message : String(err);
        }
        if (root.Shiny && root.Shiny.setInputValue)
          root.Shiny.setInputValue("ramAtlasSiftsExact", {
            request_id: requestId, state: "ok", endpoint: url, ...result
          }, {priority: "event"});
      } finally {
        clearTimeout(timer);
      }
    } catch (e) {
      if (root.Shiny && root.Shiny.setInputValue)
        root.Shiny.setInputValue("ramAtlasSiftsExact", {
          request_id: requestId, state: "error",
          pdb_id: id.pdb.toUpperCase(), entity_id: id.entity,
          accession: id.uniprot,
          message: e && e.message ? e.message : String(e)
        }, {priority: "event"});
    }
  }

  function register() {
    if (!root.Shiny || !root.Shiny.addCustomMessageHandler) return false;
    root.Shiny.addCustomMessageHandler("ram-atlas-sifts-exact", fetchExact);
    return true;
  }
  if (!register() && root.document)
    root.document.addEventListener("shiny:connected", register, {once: true});
  if (root.document) root.document.addEventListener("click", e => {
    const button = e.target && e.target.closest &&
      e.target.closest(".ram-atlas-verify");
    if (!button || !root.Shiny || !root.Shiny.setInputValue) return;
    root.Shiny.setInputValue("ramAtlasVerifyPick", {
      pdb_id: button.dataset.pdb,
      entity_id: button.dataset.entity
    }, {priority: "event"});
  });
  const api = {tokenize, parseLoops, exactRows, extractMappedCA, validate};
  root.RamplotRExactSifts = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
