(function () {
  "use strict";

  const SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query";
  const DATA_BASE = "https://data.rcsb.org/rest/v1/core";

  async function fetchJson(url, options, timeoutMs) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), timeoutMs || 20000);
    try {
      const response = await fetch(url, Object.assign({}, options || {}, {
        signal: controller.signal
      }));
      if (response.status === 204) return null;
      if (!response.ok) {
        throw new Error("HTTP " + response.status + " from " + url);
      }
      return await response.json();
    } finally {
      clearTimeout(timer);
    }
  }

  function compactText(value, fallback) {
    const text = String(value == null ? "" : value)
      .replace(/\s+/g, " ").trim();
    return text || fallback || "";
  }

  function parseEntityIdentifier(identifier) {
    const value = compactText(identifier, "");
    const match = value.match(/^([A-Za-z0-9]{4})_(\d+)$/);
    return match ? {
      pdb_id: match[1].toUpperCase(),
      entity_id: match[2]
    } : null;
  }

  function chainIds(entity) {
    const ids = entity && entity.rcsb_polymer_entity_container_identifiers;
    let chains = ids && Array.isArray(ids.auth_asym_ids)
      ? ids.auth_asym_ids.filter(Boolean).map(String) : [];
    if (!chains.length && entity && entity.entity_poly) {
      chains = compactText(entity.entity_poly.pdbx_strand_id, "")
        .split(",").map(x => x.trim()).filter(Boolean);
    }
    return [...new Set(chains)];
  }

  function entityDescription(entity) {
    return compactText(
      entity && entity.rcsb_polymer_entity &&
      entity.rcsb_polymer_entity.pdbx_description,
      "Protein polymer entity"
    );
  }

  function entryMethod(entry) {
    const values = entry && Array.isArray(entry.exptl)
      ? entry.exptl.map(item => compactText(item && item.method, ""))
        .filter(Boolean) : [];
    return [...new Set(values)].join(", ") || "Experimental structure";
  }

  function entryResolution(entry) {
    const values = entry && entry.rcsb_entry_info &&
      Array.isArray(entry.rcsb_entry_info.resolution_combined)
      ? entry.rcsb_entry_info.resolution_combined.map(Number)
        .filter(Number.isFinite) : [];
    return values.length ? Math.min(...values) : null;
  }

  function entryTitle(entry) {
    return compactText(entry && entry.struct && entry.struct.title, "");
  }

  function releaseDate(entry) {
    return compactText(
      entry && entry.rcsb_accession_info &&
      entry.rcsb_accession_info.initial_release_date, ""
    );
  }

  async function enrichHit(hit, entryCache) {
    const parsed = parseEntityIdentifier(hit && hit.identifier);
    if (!parsed) return null;
    const entityUrl = DATA_BASE + "/polymer_entity/" +
      encodeURIComponent(parsed.pdb_id) + "/" +
      encodeURIComponent(parsed.entity_id);
    const entryKey = parsed.pdb_id;
    if (!entryCache.has(entryKey)) {
      entryCache.set(entryKey, fetchJson(
        DATA_BASE + "/entry/" + encodeURIComponent(entryKey),
        { method: "GET" }, 15000
      ).catch(() => null));
    }
    const [entity, entry] = await Promise.all([
      fetchJson(entityUrl, { method: "GET" }, 15000).catch(() => null),
      entryCache.get(entryKey)
    ]);
    if (!entity) return null;
    const chains = chainIds(entity);
    return {
      pdb_id: parsed.pdb_id,
      entity_id: parsed.entity_id,
      chain: chains.length ? chains[0] : "",
      chains,
      description: entityDescription(entity),
      title: entryTitle(entry),
      method: entryMethod(entry),
      resolution: entryResolution(entry),
      release_date: releaseDate(entry),
      search_score: Number.isFinite(Number(hit && hit.score))
        ? Number(hit.score) : null
    };
  }

  async function searchExperimental(payload) {
    if (!window.Shiny || !window.Shiny.setInputValue) return;
    const requestId = compactText(payload && payload.request_id,
      String(Date.now()));
    const sequence = compactText(payload && payload.sequence, "")
      .replace(/\s+/g, "").toUpperCase();
    const identity = Number(payload && payload.identity_cutoff);
    const rows = Math.max(1, Math.min(25,
      Number(payload && payload.rows) || 12));

    window.Shiny.setInputValue("ramExperimentalSearchStatus", {
      request_id: requestId,
      state: "searching",
      message: "Searching experimental PDB polymer entities..."
    }, { priority: "event" });

    if (sequence.length < 20) {
      window.Shiny.setInputValue("ramExperimentalSearchStatus", {
        request_id: requestId,
        state: "error",
        message: "The selected chain is too short for a useful sequence search."
      }, { priority: "event" });
      return;
    }

    const body = {
      query: {
        type: "terminal",
        service: "sequence",
        parameters: {
          target: "pdb_protein_sequence",
          value: sequence,
          evalue_cutoff: 0.1,
          identity_cutoff: Number.isFinite(identity)
            ? Math.max(0, Math.min(1, identity)) : 0.9
        }
      },
      return_type: "polymer_entity",
      request_options: {
        paginate: { start: 0, rows },
        results_content_type: ["experimental"],
        scoring_strategy: "sequence"
      }
    };

    try {
      const response = await fetchJson(SEARCH_URL, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body)
      }, 30000);
      const hits = response && Array.isArray(response.result_set)
        ? response.result_set : [];
      const entryCache = new Map();
      const enriched = (await Promise.all(
        hits.slice(0, rows).map(hit => enrichHit(hit, entryCache))
      )).filter(Boolean);

      window.Shiny.setInputValue("ramExperimentalSearchResults", {
        request_id: requestId,
        total_count: response && Number.isFinite(Number(response.total_count))
          ? Number(response.total_count) : hits.length,
        results: enriched
      }, { priority: "event" });
      window.Shiny.setInputValue("ramExperimentalSearchStatus", {
        request_id: requestId,
        state: "done",
        message: enriched.length
          ? "Experimental counterparts found."
          : "No experimental counterpart met the search threshold."
      }, { priority: "event" });
    } catch (error) {
      const message = error && error.name === "AbortError"
        ? "The experimental-structure search timed out."
        : "Experimental-structure search failed: " +
          compactText(error && error.message, "unknown error");
      window.Shiny.setInputValue("ramExperimentalSearchStatus", {
        request_id: requestId,
        state: "error",
        message
      }, { priority: "event" });
    }
  }

  function install() {
    if (!window.Shiny || !window.Shiny.addCustomMessageHandler) {
      setTimeout(install, 50);
      return;
    }
    window.Shiny.addCustomMessageHandler(
      "ram-experimental-search", searchExperimental
    );
  }
  install();

  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest &&
      event.target.closest(".ram-experimental-compare");
    if (!button || !window.Shiny || !window.Shiny.setInputValue) return;
    window.Shiny.setInputValue("ramExperimentalComparePick", {
      pdb_id: compactText(button.dataset.pdb, "").toUpperCase(),
      entity_id: compactText(button.dataset.entity, ""),
      chain: compactText(button.dataset.chain, ""),
      title: compactText(button.dataset.title, "")
    }, { priority: "event" });
  });
})();