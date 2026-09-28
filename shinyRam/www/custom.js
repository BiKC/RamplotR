/* RamplotR plot and small interface interactions.
 * The server owns scientific calculations; this file only draws their results.
 */
(function () {
  "use strict";

  const defaultChainColors = [
    "#1c7777", "#9366aa", "#cf823a", "#5a8eb8",
    "#65915f", "#c25b76", "#a77e45", "#737d94"
  ];
  const plot = document.getElementById("plotly");
  const empty = document.getElementById("plot-empty");
  const currentStructure = document.getElementById("ram-current-structure");
  const selectedLabel = document.getElementById("ram-selected-residue");
  let lastPoints = [];
  let selectedResidue = null;
  let overlayIndex = -1;
  let clickBound = false;

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

    const background = [
      safeColor(shades[0], "#F1EEF6"),
      safeColor(shades[1], "#BDC9E1"),
      safeColor(shades[2], "#74A9CF"),
      safeColor(shades[3], "#0570B0")
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
    const phis = array(df.phi);
    const psis = array(df.psi);
    const insertionCodes = array(df.insertion_code);
    const points = [];
    lastPoints = points;
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
        const residue = {
          chain: String(chain), resi: Number(residueIds[row]),
          insertion_code: String(insertionCodes[row] || ""),
          resn: String(amino[row] || ""), phi: phi, psi: psi
        };
        points.push(residue);
        customdata.push(residue);
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

    // Keep one lightweight overlay for linked selections. Updating it does
    // not trigger scientific reclassification or a full contour redraw.
    overlayIndex = traces.length;
    traces.push({
      type: "scatter", mode: "markers",
      x: [], y: [], showlegend: false, hoverinfo: "skip",
      marker: { size: 15, symbol: "circle-open", color: "#ec4c2c",
                line: { width: 3, color: "#ec4c2c" } }
    });

    const narrow = plot.clientWidth < 540;
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
      autosize: true,
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
      responsive: true,
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
        if (!clickBound && typeof plot.on === "function") {
          plot.on("plotly_click", function (event) {
            const point = event && event.points && event.points[0];
            if (!point || !point.customdata || !window.Shiny) return;
            window.Shiny.setInputValue("ramplotr_point", point.customdata,
                                       { priority: "event" });
          });
          clickBound = true;
        }
        applySelection(selectedResidue);
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

  function sameResidue(a, b) {
    return a && b && String(a.chain) === String(b.chain) &&
      Number(a.resi) === Number(b.resi) &&
      String(a.insertion_code || "") === String(b.insertion_code || "");
  }

  function applySelection(residue) {
    selectedResidue = residue || null;
    if (selectedLabel) {
      selectedLabel.textContent = residue
        ? "Selected: " + (residue.resn || "Residue") + " " + residue.resi +
          (residue.insertion_code || "") + " · chain " + residue.chain
        : "Click a point or residue-table row to highlight it in 3D";
    }
    if (!plot || overlayIndex < 0 || !plot.data || !window.Plotly) return;
    const point = lastPoints.find(x => sameResidue(x, residue));
    window.Plotly.restyle(plot, {
      x: [point ? [point.phi] : []],
      y: [point ? [point.psi] : []]
    }, [overlayIndex]);
  }

  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler("process", drawPlot);
    window.Shiny.addCustomMessageHandler("ramplotr-select", applySelection);
  }

  // Resize Plotly if the sidebar or viewport changes size. Do not recreate
  // the scientific traces or reset the current plot selection on resize.
  function resizePlot() {
    if (!plot || !window.Plotly ||
        !plot.classList.contains("js-plotly-plot")) return;
    const width = plot.clientWidth;
    if (!width) return;
    const height = plotHeight();
    // Plotly can retain both the original desktop SVG width and height.
    // Update both dimensions to keep the full -180°..180° square visible.
    if (Math.abs((plot._fullLayout && plot._fullLayout.height || 0) - height) >= 2 ||
        Math.abs((plot._fullLayout && plot._fullLayout.width || 0) - width) >= 2) {
      window.Plotly.relayout(plot, { width: width, height: height }).then(function () {
        window.Plotly.Plots.resize(plot);
      });
    } else {
      window.Plotly.Plots.resize(plot);
    }
  }

  if (plot && window.ResizeObserver) {
    let previousWidth = 0;
    const observer = new ResizeObserver(function () {
      const width = plot.clientWidth;
      if (!width || Math.abs(width - previousWidth) < 5 ||
          !plot.classList.contains("js-plotly-plot")) return;
      previousWidth = width;
      resizePlot();
    });
    observer.observe(plot.parentElement || plot);
  }

  // Bootstrap 3 dispatches tab events through jQuery, not DOM EventTarget.
  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab", function () {
      resizePlot();
    });
  }
})();
