(function (root) {
  "use strict";
  const SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query";
  const DATA = "https://data.rcsb.org/rest/v1/core";
  const FIELD = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession";
  const MAX_RESULTS = 50;
  const MAX_CONCURRENT = 6;

  function clean(value, fallback = "") {
    return String(value == null ? "" : value).replace(/\s+/g, " ").trim() || fallback;
  }
  function validateAccession(accession) {
    const value = clean(accession).toUpperCase();
    if (!/^[A-Z0-9]{6,12}(?:-[1-9][0-9]*)?$/.test(value))
      throw new Error("Enter a valid UniProt accession.");
    return value;
  }
  function searchRequest(accession, rows = MAX_RESULTS, start = 0) {
    const acc = validateAccession(accession);
    const limit = Math.min(MAX_RESULTS, Math.max(1, Math.floor(Number(rows) || MAX_RESULTS)));
    const offset = Number(start);
    if (!Number.isSafeInteger(offset) || offset < 0 || offset > 1000000)
      throw new Error("Invalid Atlas page offset.");
    return {
      query: {
        type: "terminal",
        service: "text",
        parameters: {
          attribute: FIELD,
          operator: "exact_match",
          value: acc
        }
      },
      return_type: "polymer_entity",
      request_options: {
        paginate: {start: offset, rows: limit},
        sort: [{sort_by: "rcsb_id", direction: "asc"}],
        results_content_type: ["experimental"]
      }
    };
  }
  function parseEntityId(identifier) {
    const match = clean(identifier).match(/^([A-Za-z0-9]{4})_(\d+)$/);
    return match ? {pdb_id: match[1].toUpperCase(), entity_id: match[2]} : null;
  }
  function asFinite(value) {
    if (value === null || value === undefined || value === "") return null;
    const n = Number(value);
    return Number.isFinite(n) ? n : null;
  }
  function normalizeEntity(hit, entity, entry, accession) {
    const id = parseEntityId(hit && hit.identifier);
    if (!id || !entity) return null;
    const container = entity.rcsb_polymer_entity_container_identifiers || {};
    const chains = Array.isArray(container.auth_asym_ids)
      ? [...new Set(container.auth_asym_ids.map(x => clean(x)).filter(Boolean))]
      : clean(entity.entity_poly && entity.entity_poly.pdbx_strand_id)
        .split(",").map(x => x.trim()).filter(Boolean);
    const references = Array.isArray(container.reference_sequence_identifiers)
      ? container.reference_sequence_identifiers : [];
    const matching = references.filter(x =>
      clean(x && x.database_name).toLowerCase() === "uniprot" &&
      clean(x && x.database_accession).toUpperCase() ===
        clean(accession).toUpperCase());
    const coverage = matching.length === 1 ? matching[0] : null;
    const refCoverage = asFinite(coverage && coverage.reference_sequence_coverage);
    const entityCoverage = asFinite(coverage && coverage.entity_sequence_coverage);
    const resolutions = entry && entry.rcsb_entry_info &&
      Array.isArray(entry.rcsb_entry_info.resolution_combined)
      ? entry.rcsb_entry_info.resolution_combined.map(asFinite).filter(x => x !== null)
      : [];
    const methods = entry && Array.isArray(entry.exptl) ?
      entry.exptl.map(x => clean(x && x.method)).filter(Boolean) : [];
    return {
      pdb_id: id.pdb_id, entity_id: id.entity_id,
      chains,
      chain: chains[0] || "",
      description: clean(entity.rcsb_polymer_entity &&
        entity.rcsb_polymer_entity.pdbx_description, "Protein entity"),
      reference_sequence_coverage: refCoverage,
      entity_sequence_coverage: entityCoverage,
      matching_uniprot_reference: matching.length === 1,
      method: [...new Set(methods)].join(", ") || "Unknown method",
      resolution: resolutions.length ? Math.min(...resolutions) : null,
      release_date: clean(entry && entry.rcsb_accession_info &&
        entry.rcsb_accession_info.initial_release_date),
      title: clean(entry && entry.struct && entry.struct.title)
    };
  }
  async function fetchJSON(url, options = {}) {
    const controller = typeof AbortController !== "undefined"
      ? new AbortController() : null;
    const timer = controller ? setTimeout(() => controller.abort(), 15000) : null;
    try {
      const response = await fetch(url, {...options,
        credentials: "omit",
        signal: controller ? controller.signal : undefined});
      if (response.status === 204) return null;
      if (!response.ok) throw new Error("HTTP " + response.status);
      return response.json();
    } finally {
      if (timer) clearTimeout(timer);
    }
  }
  async function runBounded(items, callback) {
    const out = Array(items.length);
    let next = 0;
    async function worker() {
      while (next < items.length) {
        const i = next++;
        out[i] = await callback(items[i], i);
      }
    }
    await Promise.all(Array.from({length: Math.min(MAX_CONCURRENT, items.length)},
      () => worker()));
    return out;
  }
  function notify(inputName, value) {
    if (root.Shiny && typeof root.Shiny.setInputValue === "function")
      root.Shiny.setInputValue(inputName, value, {priority: "event"});
  }
  async function searchAtlas(payload) {
    const requestId = clean(payload && payload.request_id);
    const accession = clean(payload && payload.accession).toUpperCase();
    try {
      const query = searchRequest(accession, payload && payload.rows,
        payload && payload.start == null ? 0 : payload.start);
      const offset = query.request_options.paginate.start;
      const response = await fetchJSON(SEARCH, {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(query)
      });
      const hits = response && Array.isArray(response.result_set)
        ? response.result_set : [];
      const cache = new Map();
      const records = await runBounded(hits, async hit => {
        const id = parseEntityId(hit.identifier);
        if (!id) return null;
        if (!cache.has(id.pdb_id))
          cache.set(id.pdb_id, fetchJSON(DATA + "/entry/" + id.pdb_id)
            .catch(() => null));
        const [entity, entry] = await Promise.all([
          fetchJSON(DATA + "/polymer_entity/" + id.pdb_id + "/" + id.entity_id)
            .catch(() => null),
          cache.get(id.pdb_id)
        ]);
        return normalizeEntity(hit, entity, entry, accession);
      });
      const results = records.filter(Boolean);
      const failedEntityIds = hits.filter((hit,i) => !records[i])
        .map(hit => clean(hit && hit.identifier))
        .filter(identifier => /^[A-Za-z0-9]{4}_[1-9][0-9]*$/.test(identifier));
      notify("ramAtlasResults", {
        request_id: requestId,
        accession,
        start: offset,
        total_count: response && Number.isFinite(Number(response.total_count))
          ? Number(response.total_count) : offset + hits.length,
        returned_count: hits.length,
        next_offset: offset + hits.length,
        incomplete_metadata: hits.length - results.length,
        failed_entity_ids: failedEntityIds,
        results
      });
    } catch (err) {
      notify("ramAtlasResults", {
        request_id: requestId, accession,
        start: Number.isSafeInteger(Number(payload && payload.start))
          ? Number(payload.start) : 0,
        state: "error",
        message: err && err.message ? err.message : "Atlas search failed."
      });
    }
  }
  function setup() {
    if (!root.Shiny || typeof root.Shiny.addCustomMessageHandler !== "function")
      return false;
    root.Shiny.addCustomMessageHandler("ram-atlas-discover", searchAtlas);
    root.Shiny.addCustomMessageHandler("ram-open-group-panel", () => {
      const panel=root.document &&
        root.document.getElementById("ram-group-comparison-panel");
      if(panel) {
        panel.open=true;
        if(typeof panel.scrollIntoView==="function")
          panel.scrollIntoView({behavior:"smooth",block:"start"});
      }
    });
    return true;
  }
  if (!setup() && root.document) root.document.addEventListener("shiny:connected",
    setup, {once: true});
  if (root.document) root.document.addEventListener("click", event => {
    const button = event.target && event.target.closest &&
      event.target.closest(".ram-atlas-compare");
    if (!button) return;
    notify("ramAtlasComparePick", {
      pdb_id: clean(button.dataset.pdb).toUpperCase(),
      chain: clean(button.dataset.chain),
      entity_id: clean(button.dataset.entity)
    });
  });
  const api = {validateAccession, searchRequest, parseEntityId, normalizeEntity,
    runBounded};
  root.RamplotRAtlasDiscovery = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
