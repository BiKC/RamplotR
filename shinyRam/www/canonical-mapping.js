(function (root) {
  "use strict";

  const V2_BASE = "https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/";
  const COMPAT_BASE = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/";

  function text(value, fallback = "") {
    return value === null || value === undefined ? fallback : String(value);
  }

  function number(value) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : null;
  }

  function entryFromPayload(payload, pdbId) {
    if (!payload || typeof payload !== "object") return null;
    const wanted = String(pdbId || "").toLowerCase();
    for (const [key, value] of Object.entries(payload)) {
      if (String(key).toLowerCase() === wanted) return value;
    }
    const values = Object.values(payload);
    return values.length === 1 ? values[0] : null;
  }

  function normalizePayload(payload, pdbId) {
    const entry = entryFromPayload(payload, pdbId);
    if (!entry || typeof entry !== "object") return [];
    const uniprot = entry.UniProt || entry.uniprot;
    if (!uniprot || typeof uniprot !== "object") return [];

    const out = [];
    for (const [accession, info] of Object.entries(uniprot)) {
      if (!info || typeof info !== "object") continue;
      const mappings = Array.isArray(info.mappings) ? info.mappings : [];
      for (const mapping of mappings) {
        if (!mapping || typeof mapping !== "object") continue;
        const start = mapping.start && typeof mapping.start === "object"
          ? mapping.start : {};
        const end = mapping.end && typeof mapping.end === "object"
          ? mapping.end : {};
        out.push({
          pdb_id: String(pdbId || "").toUpperCase(),
          uniprot_accession: text(accession).toUpperCase(),
          uniprot_identifier: text(info.identifier),
          entity_id: number(mapping.entity_id),
          chain: text(mapping.chain_id),
          struct_asym_id: text(mapping.struct_asym_id),
          unp_start: number(mapping.unp_start),
          unp_end: number(mapping.unp_end),
          pdb_start: number(start.residue_number),
          pdb_end: number(end.residue_number),
          author_start: number(start.author_residue_number),
          author_end: number(end.author_residue_number),
          author_start_insertion: text(start.author_insertion_code),
          author_end_insertion: text(end.author_insertion_code),
          identity: number(mapping.identity),
          coverage: number(mapping.coverage)
        });
      }
    }
    return out;
  }

  async function fetchJson(url) {
    const controller = typeof AbortController !== "undefined"
      ? new AbortController() : null;
    const timer = controller
      ? setTimeout(() => controller.abort(), 10000) : null;
    try {
      const response = await fetch(url, {
        method: "GET",
        headers: {Accept: "application/json"},
        credentials: "omit",
        signal: controller ? controller.signal : undefined
      });
      if (!response.ok)
        throw new Error("PDBe mapping request returned " + response.status);
      return response.json();
    } finally {
      if (timer) clearTimeout(timer);
    }
  }

  async function fetchMapping(pdbId) {
    const id = String(pdbId || "").trim().toLowerCase();
    if (!/^[a-z0-9]{4}$/.test(id))
      throw new Error("Invalid PDB accession for canonical mapping.");

    const errors = [];
    for (const base of [V2_BASE, COMPAT_BASE]) {
      try {
        const payload = await fetchJson(base + encodeURIComponent(id));
        return {
          endpoint: base,
          segments: normalizePayload(payload, id)
        };
      } catch (error) {
        errors.push(error && error.message ? error.message : String(error));
      }
    }
    throw new Error(errors.join("; ") || "PDBe mapping request failed.");
  }

  function sendInput(value) {
    if (root.Shiny && typeof root.Shiny.setInputValue === "function") {
      root.Shiny.setInputValue("ramCanonicalMapping", value, {priority: "event"});
    }
  }

  function registerShinyHandler() {
    if (!root.Shiny || typeof root.Shiny.addCustomMessageHandler !== "function")
      return false;
    root.Shiny.addCustomMessageHandler("ram-canonical-map", async function (request) {
      const requestId = text(request && request.request_id);
      const pdbId = text(request && request.pdb_id).toUpperCase();
      try {
        const result = await fetchMapping(pdbId);
        sendInput({
          request_id: requestId,
          pdb_id: pdbId,
          state: "ok",
          endpoint: result.endpoint,
          segments: result.segments
        });
      } catch (error) {
        sendInput({
          request_id: requestId,
          pdb_id: pdbId,
          state: "error",
          message: error && error.message ? error.message : String(error),
          segments: []
        });
      }
    });
    return true;
  }

  let registered = registerShinyHandler();
  if (!registered && root.document) {
    root.document.addEventListener("shiny:connected", function () {
      if (!registered) registered = registerShinyHandler();
    }, {once: true});
  }

  const api = {normalizePayload, entryFromPayload};
  root.RamplotRCanonicalMapping = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
