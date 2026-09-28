/* Optional predicted-model confidence. Independent of scientific classification.
 * AlphaFold/ESMFold tracks appear only when the user loads prediction data.
 */
(function () {
  "use strict";
  let payload = null;
  let selected = null;
  let boundPlot = null;

  function key(value) {
    return value && [String(value.chain || ""), String(value.resi),
                     String(value.insertion_code || "")].join("|");
  }

  function selectedIndex() {
    if (!payload || !Array.isArray(payload.residues)) return -1;
    return payload.residues.findIndex(function (residue) {
      return key(residue) === key(selected);
    });
  }

  function render() {
    const node = document.getElementById("ram-pae-plot");
    const panel = document.getElementById("ram-confidence-panel");
    if (!node || !panel || !panel.open || !payload ||
        !Array.isArray(payload.z) || !window.Plotly) return;

    const n = payload.residues.length;
    const ticks = [];
    const labels = [];
    const step = Math.max(1, Math.ceil(n / 12));
    for (let i = 0; i < n; i += step) {
      ticks.push(i + 1);
      labels.push(String(payload.labels[i] || (i + 1)));
    }
    const chosen = selectedIndex();
    const traces = [{
      type: "heatmap",
      z: payload.z,
      x: Array.from({length: n}, (_, i) => i + 1),
      y: Array.from({length: n}, (_, i) => i + 1),
      colorscale: [[0, "#f5faf8"], [0.25, "#b6d8d3"],
                   [0.55, "#4e9f9e"], [1, "#173b53"]],
      zmin: 0, zmax: Math.max(10, Math.min(100,
        Math.ceil(Math.max.apply(null, payload.z.map(row =>
          Math.max.apply(null, row)))))),
      colorbar: {title: {text: "Å"}, thickness: 11,
                 tickfont: {size: 10}},
      hovertemplate:
        "Aligned position %{y}<br>Target position %{x}" +
        "<br>PAE %{z:.1f} Å<extra></extra>"
    }, {
      type: "scatter",
      mode: "markers",
      x: chosen < 0 ? [] : [chosen + 1],
      y: chosen < 0 ? [] : [chosen + 1],
      marker: {size: 12, symbol: "x", color: "#e67639", line: {width: 2}},
      name: "Selected residue", showlegend: false, hoverinfo: "skip"
    }];
    const width = Math.min(590, node.parentElement.clientWidth || 450);
    const layout = {
      width: width, height: Math.min(630, Math.max(330, width + 35)),
      margin: {l: 66, r: 65, t: 15, b: 70},
      paper_bgcolor: "#ffffff", plot_bgcolor: "#f8fbfa",
      font: {family: "system-ui, sans-serif", size: 11, color: "#34545c"},
      xaxis: {title: "Target residue", tickvals: ticks,
              ticktext: labels, tickangle: -40, constrain: "domain"},
      yaxis: {title: "Aligned residue", tickvals: ticks, ticktext: labels,
              autorange: "reversed", scaleanchor: "x", scaleratio: 1},
      uirevision: String(payload.total_tokens),
      showlegend: false
    };
    Promise.resolve(window.Plotly.react(node, traces, layout, {
      displaylogo: false, responsive: true
    })).then(function () {
      if (boundPlot !== node && typeof node.on === "function") {
        node.on("plotly_click", function (event) {
          if (!event || !event.points || !event.points.length || !payload) return;
          const index = Math.round(Number(event.points[0].x)) - 1;
          if (index < 0 || index >= payload.residues.length ||
              !window.Shiny || !window.Shiny.setInputValue) return;
          window.Shiny.setInputValue("ramPaePick", payload.residues[index],
                                     {priority: "event"});
        });
        boundPlot = node;
      }
    }).catch(function (error) {
      const note = document.getElementById("ram-pae-note");
      if (note) note.textContent = "PAE could not be rendered: " + error.message;
    });
    const note = document.getElementById("ram-pae-note");
    if (note) note.textContent = payload.downsampled
      ? "Visualisation downsampled to " + payload.displayed_tokens +
        " of " + payload.total_tokens + " mapped protein tokens; the original matrix is not modified."
      : "PAE is directional: rows are the aligned residues, columns are target residues.";
  }

  function deferRender() {
    if (typeof window.requestAnimationFrame === "function")
      window.requestAnimationFrame(render);
    else render();
  }

  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler("ram-confidence", function (data) {
      payload = data && !data.clear ? data : null;
      boundPlot = null;
      deferRender();
    });
    window.Shiny.addCustomMessageHandler("ram-confidence-selected", function (data) {
      selected = data && !data.clear ? data : null;
      const node = document.getElementById("ram-pae-plot");
      if (node && payload && window.Plotly && node.data) {
        const index = selectedIndex();
        window.Plotly.restyle(node, {
          x: [index < 0 ? [] : [index + 1]],
          y: [index < 0 ? [] : [index + 1]]
        }, [1]);
      }
    });
  }
  document.addEventListener("shiny:value", function (event) {
    if (event.target && event.target.id === "predictionPanel") deferRender();
  });
  document.addEventListener("toggle", function (event) {
    if (event.target && event.target.id === "ram-confidence-panel") deferRender();
  }, true);
  window.addEventListener("resize", function () {
    const node = document.getElementById("ram-pae-plot");
    if (node && payload && node.data) deferRender();
  });
})();
