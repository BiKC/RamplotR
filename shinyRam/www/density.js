/* Local cryo-EM density overlay: the user's File is parsed by NGL in the
 * browser, not uploaded to an external map-fitting service. No quantitative
 * fit assessment is inferred from simply displaying an isosurface.
 */
(function () {
  "use strict";
  let component = null, representation = null, owningStage = null;
  let generation = 0;
  const bytesLimit = 64 * 1024 * 1024;
  const find = id => document.getElementById(id);
  const status = text => {
    const node = find("ram-density-status");
    if (node) node.textContent = text;
  };
  const currentStage = () =>
    typeof window.getNGLStage === "function" ?
      window.getNGLStage("NGL") : null;
  const level = () => {
    const input = find("ram-density-level");
    const n = Number(input ? input.value : 2);
    return Number.isFinite(n) && n >= 0.5 && n <= 5 ? n : 2;
  };
  function clearMap(announce = true) {
    generation++;
    if (component && owningStage &&
        typeof owningStage.removeComponent === "function") {
      try { owningStage.removeComponent(component); } catch (_) {}
    }
    component = representation = owningStage = null;
    if (announce) status("Map removed.");
  }
  async function loadMap() {
    const input = find("ram-density-file");
    const file = input && input.files && input.files[0];
    if (!file) return status("Choose a local CCP4/MRC map first.");
    if (!/\.(map|mrc|ccp4)$/i.test(file.name))
      return status("Only .map, .mrc and .ccp4 files are supported.");
    if (!file.size || file.size > bytesLimit)
      return status("Choose a nonempty map of 64 MB or less for interactive rendering.");
    const stage = currentStage();
    if (!stage || typeof stage.loadFile !== "function")
      return status("Load a molecular structure before displaying its map.");
    clearMap(false);
    const token = generation;
    status("Reading local density map in NGL…");
    let next = null;
    try {
      // NGL's CCP4 parser supports binary CCP4/MRC density maps. Let the
      // browser read the File object directly; no server URL is generated.
      next = await stage.loadFile(file, {ext:"ccp4",defaultRepresentation:false});
      if (token !== generation || stage !== currentStage()) {
        if (next && typeof stage.removeComponent==="function")
          stage.removeComponent(next);
        return;
      }
      if (!next || typeof next.addRepresentation !== "function")
        throw new Error("NGL did not return a compatible volume component.");
      const chosen = level();
      const added = next.addRepresentation("surface", {
        isolevelType:"sigma",isolevel:chosen,opacity:0.30,
        color:"#188e9a",useWorker:true
      });
      component = next;
      owningStage = stage;
      representation = added && typeof added.setParameters==="function" ?
        added : next.reprList && next.reprList[0];
      status("Local map displayed at " + chosen.toFixed(2) +
             "σ (visual overlay, not a validation score).");
    } catch (error) {
      if (next && typeof stage.removeComponent==="function") {
        try { stage.removeComponent(next); } catch (_) {}
      }
      if (token === generation)
        status("Unable to display this map: " + (error && error.message || error));
    }
  }
  document.addEventListener("click", function (event) {
    if (!event.target || !event.target.closest) return;
    if (event.target.closest("#ram-density-load")) {
      event.preventDefault();
      loadMap();
    } else if (event.target.closest("#ram-density-clear")) {
      event.preventDefault();
      clearMap();
    }
  });
  document.addEventListener("change", function (event) {
    if (!event.target || event.target.id !== "ram-density-level") return;
    const value = level();
    const display = find("ram-density-value");
    if (display) display.textContent = value.toFixed(2) + "σ";
    if (!component || !representation ||
        typeof representation.setParameters !== "function") return;
    try {
      representation.setParameters({isolevelType:"sigma",isolevel:value});
      status("Local map displayed at " + value.toFixed(2) +
             "σ (visual overlay, not a validation score).");
    } catch (error) {
      status("Could not update the map threshold.");
    }
  });
  if (window.Shiny && typeof window.Shiny.addCustomMessageHandler==="function")
    window.Shiny.addCustomMessageHandler("ram-clear-density",function (_message) {
      clearMap(false);
      status("No map loaded.");
    });
})();
