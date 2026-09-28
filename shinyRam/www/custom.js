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
  let selectedResidue = null;
  let selectedCoordinates = new Map();
  let selectionTrace = -1;
  let boundStage = null;
  let boundPickHandler = null;
  let focusedKey = "";
  let focusedStage = null;

  // A residue identity is made of chain, sequence position and insertion
  // code. Do not focus an arbitrary selector received from a browser event.
  function nglSelection(item) {
    if (!item || item.chain == null || item.resi == null) return null;
    const chain = String(item.chain);
    const ins = String(item.insertion_code || "");
    const resi = Number(item.resi);
    if (!Number.isInteger(resi) || !/^[A-Za-z0-9_-]*$/.test(chain) ||
        !/^[A-Za-z0-9]*$/.test(ins)) return null;
    return String(resi) + (ins ? "^" + ins : "") +
      (chain ? ":" + chain : "");
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

    const key = sele ? selectionKey(selectedResidue) : "";
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

  function markSequenceSelection() {
    if (typeof document.querySelectorAll !== "function") return;
    const selected = selectionKey(selectedResidue);
    document.querySelectorAll(".ram-seq-res").forEach(function (button) {
      const key = selectionKey({
        chain: button.dataset.chain, resi: Number(button.dataset.resi),
        insertion_code: button.dataset.insertion || ""
      });
      button.setAttribute("aria-pressed", String(!!selected && selected === key));
    });
  }

  function changeSource() {
    const selected = document.querySelector('input[name="inputSource"]:checked');
    const upload = selected && selected.value === "upload";
    const pdbWrap = document.getElementById("ram-pdb-wrap");
    const uploadWrap = document.getElementById("ram-upload-wrap");
    if (pdbWrap) pdbWrap.classList.toggle("is-hidden", upload);
    if (uploadWrap) uploadWrap.classList.toggle("is-hidden", !upload);
  }

  function initSourceControl() {
    changeSource();
    document.addEventListener("change", function (event) {
      if (event.target && event.target.name === "inputSource") changeSource();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initSourceControl);
  } else {
    initSourceControl();
  }

  function plotHeight() {
    // The height follows the available panel width. The figure itself retains
    // equal scaling on phi and psi so the geometry is never distorted.
    return Math.max(plot.clientWidth < 540 ? 285 : 385,
                    Math.min(690, Math.round(plot.clientWidth +
                      (plot.clientWidth < 540 ? 35 : 75))));
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
      if (empty) {
        empty.hidden = false;
        const message = empty.querySelector("p");
        if (message) message.textContent =
          "The plotting library could not load. Check your connection and reload.";
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

  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler("process", drawPlot);
    window.Shiny.addCustomMessageHandler("ram-selection", function (choice) {
      selectedResidue = choice && !choice.clear ? choice : null;
      refreshSelection();
      markSequenceSelection();
      syncNglSelection(true);
    });
    // Shiny requires every custom message handler to declare one argument.
    window.Shiny.addCustomMessageHandler("ram-bind-ngl", function (_message) {
      bindNglPick();
      syncNglSelection(true);
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
    if (window.jQuery && window.jQuery.fn &&
        window.jQuery.fn.dataTable) {
      const api = window.jQuery.fn.dataTable.tables({
        visible: true, api: true
      });
      if (api && typeof api.columns === "function") api.columns.adjust();
    }
  }
  document.addEventListener("shown.bs.tab", function () {
    scheduleResize(); adjustVisibleTables();
  });
  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab", function () {
      scheduleResize();
      // Wait for the Bootstrap pane to finish its layout before measuring DT.
      if (window.requestAnimationFrame)
        window.requestAnimationFrame(adjustVisibleTables);
      else adjustVisibleTables();
    });
  }

})();
