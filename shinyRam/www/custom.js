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
    return Math.max(385, Math.min(690, Math.round(plot.clientWidth + 75)));
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
    const uniqueChains = [...new Set(chains.map(String))];
    const chainColors = array(obj.chainColors);
    let totalPoints = 0;

    uniqueChains.forEach(function (chain, index) {
      const x = [], y = [], label = [];
      chains.forEach(function (name, row) {
        if (String(name) !== chain) return;
        const phi = phis[row];
        const psi = psis[row];
        // Missing angles are NA/null, not (0, 0).
        if (typeof phi !== "number" || typeof psi !== "number" ||
            !Number.isFinite(phi) || !Number.isFinite(psi)) return;
        x.push(phi);
        y.push(psi);
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
        text: label,
        hovertemplate: "%{text}<extra></extra>",
        name: "Chain " + (chain || "unassigned"),
        showlegend: true,
        marker: {
          color: safeColor(chainColors[index],
                           defaultChainColors[index % defaultChainColors.length]),
          size: x.length > 2500 ? 5 : 7,
          opacity: 0.84,
          line: { color: "#ffffff", width: 0.7 }
        }
      });
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

    plot.style.display = "block";
    if (empty) empty.hidden = true;
    if (currentStructure) {
      currentStructure.textContent =
        (obj.name ? obj.name + " · " : "") +
        totalPoints.toLocaleString() + " plotted residues";
    }
    Promise.resolve(window.Plotly.react(plot, traces, layout, config))
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
  }

  // Resize Plotly if the sidebar or viewport changes size. Do not recreate
  // the scientific traces or reset the current plot selection on resize.
  if (plot && window.ResizeObserver) {
    let previousWidth = 0;
    const observer = new ResizeObserver(function () {
      const width = plot.clientWidth;
      if (!width || Math.abs(width - previousWidth) < 5 ||
          !plot.classList.contains("js-plotly-plot")) return;
      previousWidth = width;
      if (window.Plotly) window.Plotly.Plots.resize(plot);
    });
    observer.observe(plot.parentElement || plot);
  }

  document.addEventListener("shown.bs.tab", function () {
    if (plot && window.Plotly && plot.classList.contains("js-plotly-plot")) {
      window.Plotly.Plots.resize(plot);
    }
    window.dispatchEvent(new Event("resize"));
  });
})();
