/* RamplotR plot and small interface interactions.
 * The server owns scientific calculations; this file only draws their results.
 */
(function () {
  "use strict";

  const defaultChainColors = [
    "#CE6A4D", "#317E9A", "#8065A3", "#B68A3E",
    "#498777", "#B65B7A", "#5070A0", "#826B4B"
  ];
  const plot = document.getElementById("plotly");
  const empty = document.getElementById("plot-empty");
  const currentStructure = document.getElementById("ram-current-structure");
  const comparePlot = document.getElementById("comparePlot");
  let selectedResidue = null;
  let selectedCoordinates = new Map();
  let selectionTrace = -1;
  let boundStage = null;
  let boundPickHandler = null;
  let focusedKey = "";
  let focusedStage = null;
  let deferredPlot = null;
  let plotLoadPending = false;
  let deferredComparison = null;
  let comparisonLoadPending = false;

  // A residue identity is made of chain, sequence position and insertion
  // code. Do not focus an arbitrary selector received from a browser event.
  function nglSelection(item) {
    if (!item || item.chain == null || item.resi == null) return null;
    const chain = String(item.chain);
    const ins = String(item.insertion_code || "");
    const resi = Number(item.resi);
    if (!Number.isInteger(resi) || !/^[A-Za-z0-9_-]*$/.test(chain) ||
        !/^[A-Za-z0-9]*$/.test(ins)) return null;
    const residue = String(resi) + (ins ? "^" + ins : "") +
      (chain ? ":" + chain : "");
    const model = Number(item.modelIndex);
    return item.multipleModels && Number.isInteger(model) && model > 0
      ? residue + " and /" + (model - 1) : residue;
  }

  function syncNglSelection(zoom) {
    if (typeof window.getNGLStage !== "function" ||
        typeof window.getNGLStructure !== "function") return false;
    const stage = window.getNGLStage("NGL");
    const structures = window.getNGLStructure("NGL");
    if (!stage || !structures || !structures.length ||
        typeof stage.getRepresentationsByName !== "function") return false;
    const highlight = stage.getRepresentationsByName("ram-highlight");
    if (!highlight || !highlight.list || !highlight.list.length) return false;
    if (stage !== focusedStage) {
      focusedStage = stage;
      focusedKey = "";
    }
    const sele = nglSelection(selectedResidue);
    // Update the already-created orange ball-and-stick representation rather
    // than adding a new representation for every residue click.
    highlight.setSelection(sele || "none");

    const key = sele ? selectionKey(selectedResidue) + "::" +
      String(selectedResidue.modelIndex || 1) : "";
    if (sele && zoom && key !== focusedKey) {
      const component = structures.find(function (item) {
        return item && item.structure && typeof item.autoView === "function";
      });
      if (component) {
        // Pause animations before an animated camera move so the selected
        // side chain stays centered while it is being inspected.
        if (typeof stage.setRock === "function") stage.setRock(false);
        if (typeof stage.setSpin === "function") stage.setSpin(false);
        component.autoView(sele, 650);
        focusedKey = key;
      }
    } else if (!sele && focusedKey) {
      // Clearing a selection restores the whole-structure overview.
      if (typeof stage.autoView === "function") stage.autoView(450);
      focusedKey = "";
    }
    return true;
  }

  function selectionKey(item) {
    if (!item || item.chain == null || item.resi == null) return "";
    return [String(item.chain), String(item.resi),
            String(item.insertion_code || "")].join("\\r");
  }

  function selectionTraceData() {
    const point = selectedCoordinates.get(selectionKey(selectedResidue));
    return point && selectedResidue ? {
      x: [point.phi], y: [point.psi]
    } : { x: [], y: [] };
  }

  function refreshSelection() {
    if (!plot || !window.Plotly || selectionTrace < 0 ||
        !plot.classList.contains("js-plotly-plot")) return;
    const point = selectionTraceData();
    window.Plotly.restyle(plot, { x: [point.x], y: [point.y] },
                        [selectionTrace]);
  }

  function onPlotClick(event) {
    if (!event || !event.points || !event.points.length) return;
    const point = event.points[0];
    const detail = point.customdata ||
      (point.data && point.data.customdata && point.data.customdata[point.pointNumber]);
    if (!Array.isArray(detail) || detail.length < 3) return;
    const pick = {
      chain: String(detail[0]), resi: Number(detail[1]),
      insertion_code: String(detail[2] || "")
    };
    if (!Number.isInteger(pick.resi)) return;
    selectedResidue = pick;
    refreshSelection();
    if (window.Shiny && window.Shiny.setInputValue)
      window.Shiny.setInputValue("ramPlotPick", pick, { priority: "event" });
  }

  function bindNglPick() {
    if (typeof window.getNGLStage !== "function") return false;
    const stage = window.getNGLStage("NGL");
    if (!stage || !stage.signals || !stage.signals.clicked) return false;
    if (stage === boundStage) return true;
    if (boundStage && boundPickHandler &&
        boundStage.signals && boundStage.signals.clicked) {
      boundStage.signals.clicked.remove(boundPickHandler);
    }
    boundStage = stage;
    boundPickHandler = function (proxy) {
      const atom = proxy && (proxy.atom || proxy.closestBondAtom);
      if (!atom || atom.resno == null) return;
      const pick = {
        chain: String(atom.chainname || atom.chainid || ""),
        resi: Number(atom.resno),
        insertion_code: String(atom.inscode || "")
      };
      if (!Number.isInteger(pick.resi)) return;
      if (window.Shiny && window.Shiny.setInputValue)
        window.Shiny.setInputValue("ramNglPick", pick, { priority: "event" });
    };
    stage.signals.clicked.add(boundPickHandler);
    // A selection can be made before NGL finishes loading. Once the widget
    // announces it is ready, apply both the stick representation and focus.
    syncNglSelection(true);
    return true;
  }

  function safeColor(value, fallback) {
    return typeof value === "string" && /^#[0-9a-f]{6}$/i.test(value)
      ? value : fallback;
  }

  function escapeText(value) {
    return String(value == null ? "" : value).replace(/[&<>"']/g, function (char) {
      return {
        "&": "&amp;", "<": "&lt;", ">": "&gt;",
        '"': "&quot;", "'": "&#39;"
      }[char];
    });
  }

  function array(value) {
    return Array.isArray(value) ? value : value == null ? [] : [value];
  }

  // A delegated sequence click survives Shiny's HTML re-rendering, including
  // when the user switches chains or changes scientific reference datasets.
  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest &&
      event.target.closest(".ram-group-pick");
    if (!button || !window.Shiny || !window.Shiny.setInputValue) return;
    const resi = Number(button.dataset.resi);
    if (!Number.isInteger(resi)) return;
    window.Shiny.setInputValue("ramGroupComparisonPick", {
      chain: String(button.dataset.chain || ""),
      resi,
      insertion_code: String(button.dataset.insertion || "")
    }, { priority: "event" });
  });

  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest &&
      event.target.closest(".ram-ensemble-cell");
    if (!button) return;
    const resi = Number(button.dataset.resi);
    if (!Number.isInteger(resi) || !window.Shiny || !window.Shiny.setInputValue)
      return;
    window.Shiny.setInputValue("ramPredictionEnsemblePick", {
      chain: String(button.dataset.chain || ""),
      resi,
      insertion_code: String(button.dataset.insertion || "")
    }, { priority: "event" });
  });

  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest &&
      event.target.closest(".ram-seq-res");
    if (!button) return;
    const pick = {
      chain: button.dataset.chain,
      resi: Number(button.dataset.resi),
      insertion_code: button.dataset.insertion || ""
    };
    if (!Number.isInteger(pick.resi)) return;
    if (window.Shiny && window.Shiny.setInputValue)
      window.Shiny.setInputValue("ramSeqPick", pick, { priority: "event" });
  });

  // Jump to an actual residue number without counting hundreds of letters.
  // The alignment may contain gaps and structures may use non-1 numbering.
  document.addEventListener("click", function (event) {
    const trigger = event.target && event.target.closest &&
      event.target.closest(".ram-seq-jump");
    if (!trigger || !trigger.closest) return;
    const chain = trigger.closest(".ram-sequence-chain");
    if (!chain || !chain.querySelector) return;
    const field = chain.querySelector(".ram-seq-jump-input");
    const note = chain.querySelector(".ram-seq-current");
    const number = field && Number(field.value);
    if (!field || field.value.trim() === "" || !Number.isInteger(number)) {
      if (note) note.textContent = "Enter a residue number.";
      return;
    }
    const buttons = Array.from(chain.querySelectorAll(".ram-seq-res"));
    const button = buttons.find(item => Number(item.dataset.resi) === number &&
      item.dataset.chain === trigger.dataset.chain) || null;
    if (!button) {
      if (note) note.textContent = "Residue " + number + " is not in this chain.";
      return;
    }
    const strip = button.closest(".ram-sequence-grid");
    if (strip) strip.scrollLeft = Math.max(0, button.offsetLeft -
      strip.offsetLeft - strip.clientWidth / 2);
    if (button.disabled) {
      if (note) note.textContent = "Residue " + number +
        " is hidden by the current plot filters.";
    } else if (typeof button.click === "function") {
      button.click();
    }
  });
  document.addEventListener("keydown", function (event) {
    const input = event.target;
    if (event.key !== "Enter" || !input || !input.matches ||
        !input.matches(".ram-seq-jump-input")) return;
    const chain = input.closest(".ram-sequence-chain");
    const button = chain && chain.querySelector(".ram-seq-jump");
    if (button && typeof button.click === "function") {
      event.preventDefault();
      button.click();
    }
  });

  let lastSequenceScrollKey = "";
  function markSequenceSelection() {
    if (typeof document.querySelectorAll !== "function") return;
    const selected = selectionKey(selectedResidue);
    let activeButton = null;
    document.querySelectorAll(".ram-seq-res").forEach(function (button) {
      const key = selectionKey({
        chain: button.dataset.chain, resi: Number(button.dataset.resi),
        insertion_code: button.dataset.insertion || ""
      });
      const active = !!selected && selected === key;
      button.setAttribute("aria-pressed", String(active));
      if (active) activeButton = button;
    });
    document.querySelectorAll(".ram-sequence-chain").forEach(function (chain) {
      const note = chain.querySelector && chain.querySelector(".ram-seq-current");
      const active = chain.querySelector &&
        chain.querySelector('.ram-seq-res[aria-pressed="true"]');
      if (!note) return;
      note.textContent = active ?
        "Selected " + active.dataset.resi + (active.dataset.insertion || "") +
        (active.dataset.plddt ? " · pLDDT " + active.dataset.plddt : "") :
        "Select a residue";
    });
    // With all chains present, keep the selected letter in view only within
    // its own horizontal sequence row. Do not scroll the entire page.
    if (activeButton && selected !== lastSequenceScrollKey &&
        typeof activeButton.closest === "function") {
      const strip = activeButton.closest(".ram-sequence-grid");
      if (strip && strip.clientWidth && Number.isFinite(activeButton.offsetLeft)) {
        strip.scrollLeft = Math.max(0, activeButton.offsetLeft -
          strip.offsetLeft - strip.clientWidth / 2);
        lastSequenceScrollKey = selected;
      }
    }
    if (!selected) lastSequenceScrollKey = "";
  }

  document.addEventListener("click", function (event) {
    const target = event.target && event.target.closest &&
      event.target.closest(".ram-change-pick");
    if (!target) return;
    const rowId = Number(target.dataset.rowId);
    if (!Number.isInteger(rowId)) return;
    if (window.Shiny && window.Shiny.setInputValue)
      window.Shiny.setInputValue("ramCompareTrackPick", rowId,
        { priority: "event" });
  });

  function changeCompareSource() {
    const selected = document.querySelector('input[name="compareInputSource"]:checked');
    const upload = selected && selected.value === "upload";
    const pdb = document.getElementById("ram-compare-pdb");
    const file = document.getElementById("ram-compare-upload");
    if (pdb) pdb.classList.toggle("is-hidden", upload);
    if (file) file.classList.toggle("is-hidden", !upload);
  }

  // The expanded all-chain navigator may render after the selection message.
  // Re-apply selection markers whenever Shiny inserts its residue buttons.
  document.addEventListener("shiny:value", function (event) {
    if (event.target && event.target.id === "sequenceView" &&
        typeof window.requestAnimationFrame === "function") {
      window.requestAnimationFrame(markSequenceSelection);
    }
  });

  function changeSource() {
    const selected = document.querySelector('input[name="inputSource"]:checked');
    const upload = selected && selected.value === "upload";
    const afdb = selected && selected.value === "afdb";
    const pdbWrap = document.getElementById("ram-pdb-wrap");
    const uploadWrap = document.getElementById("ram-upload-wrap");
    const afdbWrap = document.getElementById("ram-afdb-wrap");
    const predictionFields = document.getElementById("ram-prediction-upload");
    if (pdbWrap) pdbWrap.classList.toggle("is-hidden", !!(upload || afdb));
    if (uploadWrap) uploadWrap.classList.toggle("is-hidden", !upload);
    if (afdbWrap) afdbWrap.classList.toggle("is-hidden", !afdb);
    if (predictionFields) {
      predictionFields.classList.toggle("is-hidden", !upload);
      const app = document.querySelector(".ram-app");
      if (upload && (!app || !app.classList.contains("ram-has-data"))) {
        predictionFields.open = true;
      }
    }
  }

  function changePredictionSource() {
    const source = document.getElementById("predictionSource");
    const sidecars = document.getElementById("ram-confidence-sidecars");
    if (!source || !sidecars) return;
    sidecars.classList.toggle("is-hidden",
      source.value !== "alphafold2" && source.value !== "alphafold3");
  }

  function initSourceControl() {
    changeSource();
    changePredictionSource();
    changeCompareSource();
    document.addEventListener("change", function (event) {
      if (event.target && event.target.name === "inputSource") changeSource();
      if (event.target && event.target.id === "predictionSource")
        changePredictionSource();
      if (event.target && event.target.name === "compareInputSource")
        changeCompareSource();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initSourceControl);
  } else {
    initSourceControl();
  }

  function plotHeight() {
    // Square angular axes remain enforced in Plotly. On short laptops, use
    // the available vertical space rather than blindly matching card width.
    const width = plot.clientWidth;
    // Reserve enough physical pixels for equal 360° axes and their labels.
    // The previous fixed 385px minimum shrank the angular area to <48%
    // of its width on compact laptop screens with prediction controls.
    const minimum = width < 540 ? 285 :
      Math.max(385, Math.min(690, Math.ceil(width * 0.52 + 124)));
    let height = Math.max(minimum, Math.min(690,
      Math.round(width + (width < 540 ? 35 : 75))));
    const app = document.querySelector && document.querySelector(".ram-app");
    if (app && app.classList.contains("ram-has-data") &&
        (window.innerWidth || 0) >= 900 && (window.innerHeight || 0) < 1000 &&
        typeof plot.getBoundingClientRect === "function") {
      const top = plot.getBoundingClientRect().top;
      if (top > 0 && top < window.innerHeight) {
        height = Math.min(height,
          Math.max(minimum, Math.floor(window.innerHeight - top - 30)));
      }
    }
    return height;
  }

  function contourTrace(matrix, operation, cutoff, color) {
    return {
      type: "contour",
      x: matrix.x, y: matrix.y, z: matrix.z,
      contours: { type: "constraint", operation: operation, value: cutoff },
      fillcolor: color,
      hoverinfo: "skip",
      showscale: false,
      showlegend: false,
      line: { width: 0 }
    };
  }

  function drawPlot(obj) {
    if (!plot) return;
    if (!window.Plotly) {
      // A structure has arrived before the plotting library. Keep only the
      // latest reactive update; palette/filter changes may arrive meanwhile.
      deferredPlot = obj;
      if (empty) {
        empty.hidden = false;
        const message = empty.querySelector("p");
        if (message) message.textContent = "Loading the plot renderer…";
      }
      if (!plotLoadPending) {
        plotLoadPending = true;
        window.ramLoadPlotly().then(function () {
          plotLoadPending = false;
          const latest = deferredPlot;
          deferredPlot = null;
          if (latest) drawPlot(latest);
        }).catch(function () {
          plotLoadPending = false;
          if (empty) {
            const message = empty.querySelector("p");
            if (message) message.textContent =
              "The plotting library could not load. Check your connection and try again.";
          }
        });
      }
      return;
    }

    const frame = obj && obj.matrix;
    const limits = array(obj && obj.limits);
    const shades = array(obj && obj.backgroundColors);
    if (!frame || !frame.x || !frame.y || !frame.z || limits.length < 3) return;

    // The plot container starts hidden while empty. Measure its actual width
    // only after showing it, otherwise the first Plotly layout becomes tiny.
    plot.style.display = "block";
    plot.style.height = plotHeight() + "px";

    const background = [
      safeColor(shades[0], "#FFF8ED"),
      safeColor(shades[1], "#D4ECE7"),
      safeColor(shades[2], "#7DB9B5"),
      safeColor(shades[3], "#126E74")
    ];
    // Preserve the original scientific contour cutoffs and layer order.
    const traces = [
      contourTrace(frame, "<", limits[2], background[1]),
      contourTrace(frame, "<", limits[1], background[2]),
      contourTrace(frame, "<", limits[0], background[3]),
      contourTrace(frame, ">", limits[2], background[0])
    ];

    const df = obj.df || {};
    const chains = array(df.chain);
    const amino = array(df.resn);
    const residueIds = array(df.resi);
    const insertionCodes = array(df.insertion_code);
    selectedCoordinates = new Map();
    const phis = array(df.phi);
    const psis = array(df.psi);
    const uniqueChains = [...new Set(chains.map(String))];
    const chainColors = array(obj.chainColors);
    let totalPoints = 0;

    uniqueChains.forEach(function (chain, index) {
      const x = [], y = [], label = [], customdata = [];
      chains.forEach(function (name, row) {
        if (String(name) !== chain) return;
        const phi = phis[row];
        const psi = psis[row];
        // Missing angles are NA/null, not (0, 0).
        if (typeof phi !== "number" || typeof psi !== "number" ||
            !Number.isFinite(phi) || !Number.isFinite(psi)) return;
        x.push(phi);
        y.push(psi);
        const detail = [chain, residueIds[row], insertionCodes[row] || "",
                        amino[row]];
        customdata.push(detail);
        selectedCoordinates.set(selectionKey({
          chain: chain, resi: residueIds[row],
          insertion_code: insertionCodes[row] || ""
        }), { phi: phi, psi: psi });
        label.push(
          "<b>Chain " + escapeText(chain || "unassigned") + "</b><br>" +
          escapeText(amino[row]) + " " + escapeText(residueIds[row]) + "<br>" +
          "φ " + phi.toFixed(1) + "° · ψ " + psi.toFixed(1) + "°"
        );
      });
      if (!x.length) return;
      totalPoints += x.length;
      traces.push({
        type: x.length > 2500 ? "scattergl" : "scatter",
        mode: "markers",
        x: x, y: y,
        text: label, customdata: customdata,
        hovertemplate: "%{text}<extra></extra>",
        name: "Chain " + (chain || "unassigned"),
        // Large assemblies can contain hundreds of chains; avoid a legend
        // that covers most of the scientific plotting area.
        showlegend: uniqueChains.length <= 8,
        marker: {
          color: safeColor(chainColors[index],
                           defaultChainColors[index % defaultChainColors.length]),
          size: x.length > 2500 ? 5 : 7,
          opacity: 0.84,
          line: { color: "#ffffff", width: 0.7 }
        }
      });
    });

    // Keep a dedicated overlay trace for linked selection instead of
    // mutating the original chain markers or the scientific density layers.
    const active = selectionTraceData();
    traces.push({
      type: "scatter", mode: "markers",
      x: active.x, y: active.y,
      showlegend: false, hoverinfo: "skip",
      marker: {
        size: 17, color: "rgba(255,160,48,0.65)",
        line: { color: "#142b35", width: 2 }
      }
    });
    selectionTrace = traces.length - 1;

    const narrow = plot.clientWidth < 540;
    plot.style.height = plotHeight() + "px";
    const axis = {
      range: [-180, 180],
      tickvals: [-180, -90, 0, 90, 180],
      tickfont: { color: "#526875", size: narrow ? 10 : 11 },
      gridcolor: "#e7edf0",
      gridwidth: 1,
      zeroline: true,
      zerolinecolor: "#bdcbd0",
      zerolinewidth: 1,
      linecolor: "#b9c9ce",
      ticks: "outside",
      ticklen: 4,
      tickcolor: "#a8b9be",
      showline: true,
      automargin: true
    };
    const layout = {
      // Explicit dimensions prevent Plotly's built-in resize handler from
      // replacing the requested height with CSS min-height on narrow screens.
      autosize: false,
      width: plot.clientWidth,
      height: plotHeight(),
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#fafcfc",
      font: {
        family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif',
        color: "#263f4b", size: 12
      },
      margin: {
        l: narrow ? 49 : 62,
        r: 14, t: uniqueChains.length ? 72 : 24, b: 52
      },
      xaxis: Object.assign({}, axis, {
        // Keep the scientific -180°..180° domain instead of allowing Plotly
        // to extend the x-range when it enforces equal pixel scales.
        constrain: "domain",
        title: { text: "Phi (φ), degrees", standoff: 9,
                 font: { size: narrow ? 11 : 12 } }
      }),
      yaxis: Object.assign({}, axis, {
        title: { text: "Psi (ψ), degrees", standoff: 10,
                 font: { size: narrow ? 11 : 12 } },
        scaleanchor: "x", scaleratio: 1, constrain: "domain"
      }),
      legend: {
        orientation: "h", x: 0, y: 1.04, yanchor: "bottom",
        font: { size: 11, color: "#3c5461" },
        bgcolor: "rgba(255,255,255,0)", itemsizing: "constant"
      },
      hovermode: "closest",
      hoverlabel: {
        bgcolor: "#19303b", bordercolor: "#19303b",
        font: { color: "#ffffff", size: 12 }
      },
      uirevision: "ramplotr-geometry"
    };
    const config = {
      // The dedicated ResizeObserver sets both dimensions. Plotly's built-in
      // responsive listener races with it and restores the old desktop SVG.
      responsive: false,
      displaylogo: false,
      toImageButtonOptions: {
        format: "png",
        filename: String(obj.name || "RamplotR").replace(/[^a-z0-9_-]/gi, "_"),
        scale: 2
      }
    };

    if (empty) empty.hidden = true;
    if (currentStructure) {
      currentStructure.textContent =
        (obj.name ? obj.name + " · " : "") +
        totalPoints.toLocaleString() + " plotted residues";
    }
    Promise.resolve(window.Plotly.react(plot, traces, layout, config))
      .then(function () {
        if (typeof plot.on === "function" && !plot.__ramPickBound) {
          plot.on("plotly_click", onPlotClick);
          plot.__ramPickBound = true;
        }
        // Keep the full introductory form for first-time visitors. Once a
        // structure has rendered, the compact toolbar makes room for science.
        const app = document.querySelector && document.querySelector(".ram-app");
        if (app && !app.classList.contains("ram-has-data")) {
          app.classList.add("ram-has-data");
          const details = document.getElementById("ram-prediction-upload");
          if (details) details.open = false;
          scheduleResize();
        }
      })
      .catch(function () {
        plot.style.display = "none";
        if (empty) {
          empty.hidden = false;
          const message = empty.querySelector("p");
          if (message) message.textContent =
            "The plot could not be rendered. Choose a different structure or reload.";
        }
      });
  }

  let comparisonRows = null;
  let comparisonSelectedId = null;
  let comparisonOverlayIndices = [2, 3];
  function compareSelectedPoint(side) {
    if (!comparisonRows || comparisonSelectedId == null)
      return { x: [], y: [] };
    const ids = array(comparisonRows.rowIds);
    const index = ids.findIndex(id => Number(id) === Number(comparisonSelectedId));
    if (index < 0) return { x: [], y: [] };
    const phi = array(side === "a" ? comparisonRows.phiA : comparisonRows.phiB)[index];
    const psi = array(side === "a" ? comparisonRows.psiA : comparisonRows.psiB)[index];
    return Number.isFinite(phi) && Number.isFinite(psi) ?
      { x: [phi], y: [psi] } : { x: [], y: [] };
  }
  function emphasizeComparison() {
    if (!comparePlot || !window.Plotly || !comparePlot.data) return;
    const a = compareSelectedPoint("a"), b = compareSelectedPoint("b");
    const ia = comparisonOverlayIndices[0], ib = comparisonOverlayIndices[1];
    if (!Number.isInteger(ia) || !Number.isInteger(ib) ||
        comparePlot.data.length <= Math.max(ia, ib)) return;
    window.Plotly.restyle(comparePlot, {
      x: [a.x], y: [a.y],
      customdata: [a.x.length ? [[Number(comparisonSelectedId), "a"]] : []]
    }, [ia]);
    window.Plotly.restyle(comparePlot, {
      x: [b.x], y: [b.y],
      customdata: [b.x.length ? [[Number(comparisonSelectedId), "b"]] : []]
    }, [ib]);
  }

  function onComparePlotClick(event) {
    const item = event && event.points && event.points[0];
    const data = item && (item.customdata ||
      (item.data && item.data.customdata &&
       item.data.customdata[item.pointNumber]));
    const index = Array.isArray(data) ? Number(data[0]) : NaN;
    if (!Number.isInteger(index) || index < 1) return;
    comparisonSelectedId = index;
    emphasizeComparison();
    if (window.Shiny && window.Shiny.setInputValue)
      window.Shiny.setInputValue("ramComparePlotPick", index,
                                 { priority: "event" });
  }

  function drawComparison(obj) {
    if (!comparePlot || !obj) return;
    comparisonRows = obj;
    if (!window.Plotly) {
      deferredComparison = obj;
      if (!comparisonLoadPending) {
        comparisonLoadPending = true;
        window.ramLoadPlotly().then(function () {
          comparisonLoadPending = false;
          const latest = deferredComparison;
          deferredComparison = null;
          if (latest) drawComparison(latest);
        }).catch(function (error) {
          comparisonLoadPending = false;
          console.error("RamplotR comparison plot:", error);
        });
      }
      return;
    }
    const ids = array(obj.rowIds);
    const aPhi = array(obj.phiA), aPsi = array(obj.psiA);
    const bPhi = array(obj.phiB), bPsi = array(obj.psiB);
    const meta = {
      a: {chain:array(obj.chainA),pos:array(obj.posA),
          ins:array(obj.insA),amino:array(obj.aminoA),
          plddt:array(obj.plddtA),confidence:array(obj.confidenceA)},
      b: {chain:array(obj.chainB),pos:array(obj.posB),
          ins:array(obj.insB),amino:array(obj.aminoB),
          plddt:array(obj.plddtB),confidence:array(obj.confidenceB)}
    };
    const deltaPlddt = array(obj.deltaPlddt);
    const makeTrace = function (side, phi, psi, name, color, symbol) {
      const x = [], y = [], text = [], customdata = [];
      const m = meta[side];
      for (let i = 0; i < ids.length; i++) {
        if (!Number.isFinite(phi[i]) || !Number.isFinite(psi[i])) continue;
        const label = String(m.amino[i] || "") + " " +
          String(m.chain[i] || "") + ":" + String(m.pos[i] || "") +
          String(m.ins[i] || "");
        const partner = meta[side === "a" ? "b" : "a"];
        const other = partner.pos[i] == null ? "alignment gap" :
          String(partner.amino[i] || "") + " " + String(partner.chain[i] || "") +
          ":" + String(partner.pos[i]) + String(partner.ins[i] || "");
        const dphi = array(obj.deltaPhi)[i], dpsi = array(obj.deltaPsi)[i];
        const fmt = v => Number.isFinite(v) ? v.toFixed(1) + "°" : "N/A";
        x.push(phi[i]); y.push(psi[i]);
        const plddt = Number(m.plddt[i]);
        const conf = String(m.confidence[i] || "");
        const dplddt = Number(deltaPlddt[i]);
        const confidenceText = Number.isFinite(plddt)
          ? "<br>pLDDT " + plddt.toFixed(1) +
            (conf ? " · " + escapeText(conf) : "") : "";
        const deltaConfidenceText = Number.isFinite(dplddt)
          ? "<br>ΔpLDDT " + (dplddt >= 0 ? "+" : "") + dplddt.toFixed(1) : "";
        text.push(escapeText(label) + " ↔ " + escapeText(other) +
          "<br>φ " + fmt(phi[i]) + " · ψ " + fmt(psi[i]) +
          "<br>Δφ " + fmt(dphi) + " · Δψ " + fmt(dpsi) +
          confidenceText + deltaConfidenceText);
        customdata.push([Number(ids[i]), side]);
      }
      return {type:"scattergl",mode:"markers",name:name,x:x,y:y,
        text:text,customdata:customdata,
        hovertemplate:"%{text}<extra>" + escapeText(name) + "</extra>",
        marker:{color:color,size:8,opacity:.82,symbol:symbol}};
    };
    const frame = obj && obj.matrix;
    const limits = array(obj && obj.limits);
    const shades = array(obj && obj.backgroundColors);
    const traces = [];
    if (frame && frame.x && frame.y && frame.z && limits.length >= 3) {
      const background = [
        safeColor(shades[0], "#FFF8ED"),
        safeColor(shades[1], "#D4ECE7"),
        safeColor(shades[2], "#7DB9B5"),
        safeColor(shades[3], "#126E74")
      ];
      // Reuse the same scientific density map, contour thresholds and palette
      // as the individual Ramachandran view. Comparison points sit above it.
      traces.push(
        contourTrace(frame, "<", limits[2], background[1]),
        contourTrace(frame, "<", limits[1], background[2]),
        contourTrace(frame, "<", limits[0], background[3]),
        contourTrace(frame, ">", limits[2], background[0])
      );
    }
    traces.push(
      makeTrace("a",aPhi,aPsi,String(obj.nameA || "Primary"),"#CE6A4D","circle"),
      makeTrace("b",bPhi,bPsi,String(obj.nameB || "Comparison"),"#317E9A","diamond")
    );
    comparisonOverlayIndices = [traces.length, traces.length + 1];
    traces.push(
      {type:"scattergl",mode:"markers",showlegend:false,hoverinfo:"skip",
       x:[],y:[],marker:{color:"#ffc04a",size:17,opacity:.95,
                           line:{color:"#623c17",width:2}}},
      {type:"scattergl",mode:"markers",showlegend:false,hoverinfo:"skip",
       x:[],y:[],marker:{color:"#83e6f5",size:17,opacity:.95,
                           symbol:"diamond",line:{color:"#174958",width:2}}}
    );
    // Overlay markers retain the same aligned-pair identity when clicked.
    const index = ids.findIndex(id => Number(id) === Number(comparisonSelectedId));
    if (comparisonSelectedId != null && index >= 0) {
      for (const [side, traceIndex] of [
        ["a",comparisonOverlayIndices[0]],["b",comparisonOverlayIndices[1]]
      ]) {
        const trace = traces[traceIndex];
        const p = compareSelectedPoint(side);
        trace.x = p.x; trace.y = p.y;
        trace.customdata = p.x.length ? [[Number(ids[index]),side]] : [];
      }
    }
    const axis = {range:[-180,180],tickvals:[-180,-90,0,90,180],
      gridcolor:"#e3eeeb",zerolinecolor:"#a0bab9",constrain:"domain"};
    comparePlot.style.minHeight = "410px";
    Promise.resolve(window.Plotly.react(comparePlot,traces,{
      autosize:true,paper_bgcolor:"#ffffff",plot_bgcolor:"#fbfdfc",
      margin:{l:59,r:14,t:46,b:58},
      xaxis:Object.assign({},axis,{title:"Phi (°)"}),
      yaxis:Object.assign({},axis,{title:"Psi (°)",scaleanchor:"x",scaleratio:1}),
      legend:{orientation:"h",y:1.13,x:0},
      annotations: obj.backgroundName ? [{
        x:1,y:-0.16,xref:"paper",yref:"paper",xanchor:"right",
        showarrow:false,font:{size:10,color:"#6b7f84"},
        text:"Background: " + escapeText(String(obj.backgroundName))
      }] : [],
      height:Math.min(650,Math.max(410,comparePlot.clientWidth+45))
    },{responsive:true,displaylogo:false})).then(function () {
      if (typeof comparePlot.on === "function" && !comparePlot.__ramCompareBound) {
        comparePlot.on("plotly_click",onComparePlotClick);
        comparePlot.__ramCompareBound = true;
      }
    }).catch(function (error) {
      console.error("RamplotR comparison plot:",error);
    });
  }

  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler("process", drawPlot);
    window.Shiny.addCustomMessageHandler("ram-comparison", drawComparison);
    window.Shiny.addCustomMessageHandler("ram-comparison-selected", function (choice) {
      comparisonSelectedId = choice && !choice.clear ? Number(choice.rowId) : null;
      emphasizeComparison();
    });
    window.Shiny.addCustomMessageHandler("ram-selection", function (choice) {
      selectedResidue = choice && !choice.clear ? choice : null;
      const inspector = document.querySelector(".ram-global-inspector");
      if (inspector) inspector.classList.toggle("is-empty", !selectedResidue);
      const reviewButton = document.getElementById("nextReview");
      if (reviewButton) reviewButton.textContent =
        selectedResidue ? "Next issue" : "Review issues";
      refreshSelection();
      markSequenceSelection();
      syncNglSelection(true);
    });
    // Shiny requires every custom message handler to declare one argument.
    window.Shiny.addCustomMessageHandler("ram-bind-ngl", function (message) {
      bindNglPick();
      syncNglSelection(true);
      // NGLVieweR reuses the camera between renderValue calls. An old
      // residue-focused camera can clip a newly loaded multi-chain protein
      // even after its representations finish loading. Only reframe after
      // an actual structure load, not after every representation change.
      if (message && message.resetView === true && !selectedResidue &&
          typeof window.getNGLStage === "function") {
        const stage = window.getNGLStage("NGL");
        if (stage && typeof stage.autoView === "function") {
          if (typeof stage.handleResize === "function") stage.handleResize();
          stage.autoView(450);
        }
      }
    });
  }

  // Resize Plotly if the sidebar or viewport changes size. Do not recreate
  // the scientific traces or reset the current plot selection on resize.
  function resizePlot() {
    if (!plot || !window.Plotly ||
        !plot.classList.contains("js-plotly-plot")) return;
    const width = plot.clientWidth;
    if (!width) return;
    const height = plotHeight();
    // Give the host element a definite height. Plotly's responsive resize
    // otherwise remeasures CSS min-height (270px on mobile) and silently
    // undoes the intended square-plot dimensions.
    plot.style.height = height + "px";
    // Relayout both dimensions ourselves. Calling Plots.resize afterward
    // overrides the explicit height with the CSS minimum (270 px on mobile)
    // and collapses the scientific plotting area.
    if (Math.abs((plot._fullLayout && plot._fullLayout.height || 0) - height) >= 2 ||
        Math.abs((plot._fullLayout && plot._fullLayout.width || 0) - width) >= 2) {
      window.Plotly.relayout(plot, {
        autosize: false, width: width, height: height
      });
    }
  }

  // Both a window resize and a change in the plotting card can affect plot
  // geometry. Resize after layout has settled to avoid racing browser reflow
  // or the NGL widget's own resize handler.
  let pendingResize = null;
  function scheduleResize() {
    if (typeof window.requestAnimationFrame === "function") {
      window.requestAnimationFrame(function () {
        window.requestAnimationFrame(resizePlot);
      });
    } else {
      resizePlot();
    }
    if (typeof window.setTimeout === "function") {
      if (pendingResize !== null && typeof window.clearTimeout === "function")
        window.clearTimeout(pendingResize);
      pendingResize = window.setTimeout(function () {
        pendingResize = null;
        resizePlot();
      }, 160);
    }
  }

  if (plot && window.ResizeObserver) {
    let previousWidth = 0;
    const observer = new ResizeObserver(function () {
      const width = plot.clientWidth;
      if (!width || !plot.classList.contains("js-plotly-plot")) return;
      if (Math.abs(width - previousWidth) >= 2 ||
          !plot._fullLayout ||
          Math.abs(plot._fullLayout.width - width) >= 2) {
        previousWidth = width;
        scheduleResize();
      }
    });
    observer.observe(plot);
    if (plot.parentElement) observer.observe(plot.parentElement);
  }
  if (typeof window.addEventListener === "function")
    window.addEventListener("resize", scheduleResize);

  // When the viewport changes while the plot tab is hidden, the window
  // resize listener correctly skips its zero-width plot. Bootstrap can then
  // activate the tab without another resize event. Watch the pane itself so
  // the square plot is always measured again after it becomes visible.
  if (plot && typeof plot.closest === "function" &&
      typeof window.MutationObserver === "function") {
    const pane = plot.closest(".tab-pane");
    if (pane) {
      const tabObserver = new window.MutationObserver(function () {
        if (pane.classList.contains("active")) scheduleResize();
      });
      tabObserver.observe(pane, {
        attributes: true, attributeFilter: ["class"]
      });
    }
  }

  // Bootstrap hides the DT at initialisation. Recalculate *the same table's*
  // columns on visibility changes: do not enable DataTables scrollX, which
  // duplicates the header and creates cross-version alignment problems.
  function adjustVisibleTables() {
    if (comparePlot && window.Plotly &&
        comparePlot.classList.contains("js-plotly-plot") &&
        comparePlot.clientWidth > 0) {
      window.Plotly.Plots.resize(comparePlot);
    }
    if (window.jQuery && window.jQuery.fn &&
        window.jQuery.fn.dataTable) {
      const api = window.jQuery.fn.dataTable.tables({
        visible: true, api: true
      });
      if (api && typeof api.columns === "function") api.columns.adjust();
    }
  }
  // Show or hide analysis settings without duplicating the plots, breaking
  // Shiny inputs, or forcing users to scroll through the sidebar.
  const settingsToggle = document.getElementById("ram-toggle-settings");
  const ramApp = document.querySelector && document.querySelector(".ram-app");
  if (settingsToggle && ramApp) {
    settingsToggle.addEventListener("click", function () {
      const collapsed = ramApp.classList.toggle("ram-focus-mode");
      settingsToggle.textContent = collapsed ? "Show settings" : "Hide settings";
      settingsToggle.setAttribute("aria-expanded", String(!collapsed));
      settingsToggle.title = collapsed
        ? "Show structure input and analysis settings"
        : "Expand the plots by hiding analysis settings";
      scheduleResize();
      adjustVisibleTables();
    });
  }

  document.addEventListener("shown.bs.tab", function () {
    scheduleResize(); adjustVisibleTables(); markSequenceSelection();
  });
  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab", function () {
      scheduleResize();
      markSequenceSelection();
      // Wait for the Bootstrap pane to finish its layout before measuring DT.
      if (window.requestAnimationFrame)
        window.requestAnimationFrame(adjustVisibleTables);
      else adjustVisibleTables();
    });
  }

})();
