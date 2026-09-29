/* Linked 3D comparison. The NGL widget owns loading and superposition;
 * RamplotR only frames the two selected chains and highlights aligned pairs.
 */
(function (window, document) {
  "use strict";
  let config = null;
  let selected = null;
  let stage = null;
  let components = null;
  let boundStage = null;
  let pickHandler = null;
  let needsFit = true;
  let retry = 0;
  let fitScheduled = false;
  let framedComponents = null;
  let highlightedComponents = null;
  let highlightKey = "";

  function ready() {
    if (typeof window.getNGLStage !== "function" ||
        typeof window.getNGLStructure !== "function") return false;
    const currentStage = window.getNGLStage("NGLCompare");
    const loaded = window.getNGLStructure("NGLCompare");
    // NGLVieweR replaces stage components asynchronously when the two inputs
    // or the selected chains change. Do not operate on its previous render.
    if (!currentStage || !Array.isArray(loaded) || loaded.length < 2 ||
        !loaded[0] || !loaded[1] || !loaded[0].structure ||
        !loaded[1].structure) return false;
    if (Array.isArray(currentStage.compList) &&
        (!currentStage.compList.includes(loaded[0]) ||
         !currentStage.compList.includes(loaded[1]))) return false;
    stage = currentStage;
    components = loaded;
    return true;
  }

  function validChain(value) {
    return typeof value === "string" && /^[A-Za-z0-9_-]{1,16}$/.test(value);
  }

  function chainSelection(side) {
    if (!config) return null;
    const chain = side === "a" ? config.chainA : config.chainB;
    if (!validChain(chain)) return null;
    const model = side === "a" ? config.modelA : config.modelB;
    const multiple = side === "a" ? config.multipleA : config.multipleB;
    const index = Number(model);
    return ":" + chain + " and protein" +
      (multiple && Number.isInteger(index) && index >= 1
        ? " and /" + (index - 1) : "");
  }

  function residueSelection(item) {
    if (!item || !validChain(item.chain)) return null;
    const number = Number(item.resi);
    const insertion = String(item.insertion_code || "");
    if (!Number.isInteger(number) || !/^[A-Za-z0-9]*$/.test(insertion))
      return null;
    const model = Number(item.modelIndex);
    return String(number) + (insertion ? "^" + insertion : "") +
      ":" + item.chain + " and protein" +
      (item.multipleModels && Number.isInteger(model) && model >= 1
        ? " and /" + (model - 1) : "");
  }

  // getBoundingBox honours NGL's atom selection, unlike autoView() over the
  // entire downloaded PDB/mmCIF, which can include many hidden chains.
  function selectedBox(selections) {
    if (!stage || !components || !window.NGL || !window.NGL.Selection)
      return null;
    let result = null;
    for (let i = 0; i < 2; i++) {
      if (!selections[i]) continue;
      const structure = components[i].structure;
      if (!structure || typeof structure.getBoundingBox !== "function") continue;
      const box = structure.getBoundingBox(new window.NGL.Selection(selections[i]));
      if (!box || (typeof box.isEmpty === "function" && box.isEmpty()))
        continue;
      if (!result) result = box.clone();
      else result.union(box);
    }
    return result;
  }

  function frame(selections, duration, padding) {
    if (!ready() || !config) return false;
    const holder = document.getElementById("NGLCompare");
    if (!holder || holder.clientWidth < 10 || holder.clientHeight < 10)
      return false; // Bootstrap tab or checkbox has not become visible yet.
    if (typeof stage.handleResize === "function") stage.handleResize();
    const box = selectedBox(selections);
    if (!box || !stage.animationControls ||
        typeof stage.animationControls.zoomMove !== "function" ||
        typeof stage.getZoomForBox !== "function") return false;
    const center = box.getCenter(stage.getCenter());
    const zoom = stage.getZoomForBox(box);
    if (!Number.isFinite(zoom)) return false;
    stage.animationControls.zoomMove(center, zoom * padding, duration);
    framedComponents = components;
    needsFit = false;
    return true;
  }

  function fitBoth() {
    const a = chainSelection("a");
    const b = chainSelection("b");
    return frame([a, b], 460, 1.14);
  }

  function focusPair() {
    if (!selected) return fitBoth();
    const a = residueSelection(selected.a);
    const b = residueSelection(selected.b);
    return frame([a, b], 560, 1.55) || fitBoth();
  }

  function setHighlight(name, selection) {
    if (!stage || typeof stage.getRepresentationsByName !== "function")
      return;
    const list = stage.getRepresentationsByName(name);
    if (list && typeof list.setSelection === "function")
      list.setSelection(selection || "none");
  }

  function paintSelection(zoom) {
    if (!ready()) return false;
    const a = residueSelection(selected && selected.a);
    const b = residueSelection(selected && selected.b);
    const key = (a || "none") + "|" + (b || "none");
    // setSelection may rebuild WebGL representations and emit NGL_rendering.
    // Do not set identical highlights after every readiness notification.
    if (components !== highlightedComponents || key !== highlightKey) {
      setHighlight("ram-compare-highlight-a", a);
      setHighlight("ram-compare-highlight-b", b);
      highlightedComponents = components;
      highlightKey = key;
    }
    if (zoom) return focusPair();
    return true;
  }

  function bindClicks() {
    if (!ready() || !stage.signals || !stage.signals.clicked) return;
    if (boundStage === stage) return;
    if (boundStage && pickHandler && boundStage.signals &&
        boundStage.signals.clicked) {
      boundStage.signals.clicked.remove(pickHandler);
    }
    boundStage = stage;
    pickHandler = function (proxy) {
      const atom = proxy && (proxy.atom || proxy.closestBondAtom);
      if (!atom || atom.resno == null || !components) return;
      const side = proxy.component === components[0] ? "a" :
        proxy.component === components[1] ? "b" : null;
      if (!side || !config) return;
      const chain = String(atom.chainname || atom.chainid || "");
      const expected = side === "a" ? config.chainA : config.chainB;
      const position = Number(atom.resno);
      if (chain !== expected || !Number.isInteger(position)) return;
      if (window.Shiny && window.Shiny.setInputValue) {
        window.Shiny.setInputValue("ramCompareNglPick", {
          side: side, chain: chain, resi: position,
          insertion_code: String(atom.inscode || "")
        }, {priority: "event"});
      }
    };
    stage.signals.clicked.add(pickHandler);
  }

  function refresh() {
    fitScheduled = false;
    if (ready()) {
      bindClicks();
      if (needsFit) {
        paintSelection(true);
      } else {
        paintSelection(false);
      }
      retry = 0;
    } else if (retry++ < 12 && typeof window.setTimeout === "function") {
      window.setTimeout(scheduleFit, 170);
    }
  }

  function scheduleFit() {
    if (fitScheduled) return;
    fitScheduled = true;
    if (typeof window.requestAnimationFrame === "function") {
      window.requestAnimationFrame(function () {
        window.requestAnimationFrame(refresh);
      });
    } else refresh();
  }

  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler("ram-compare-config", function (value) {
      config = value || null;
      needsFit = true;
      retry = 0;
      scheduleFit();
    });
    window.Shiny.addCustomMessageHandler("ram-compare-ready", function (_value) {
      const loaded = typeof window.getNGLStructure === "function" ?
        window.getNGLStructure("NGLCompare") : null;
      // Rendering and highlight changes can both emit readiness. Reframe
      // only for new structure components or an explicitly pending focus.
      if (!loaded || loaded !== framedComponents) needsFit = true;
      retry = 0;
      scheduleFit();
    });
    window.Shiny.addCustomMessageHandler("ram-compare-pair", function (value) {
      selected = value && !value.clear ? value : null;
      needsFit = true;
      scheduleFit();
    });
  }

  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest &&
      event.target.closest("#compareResetView");
    if (!button) return;
    // Preserve the selected pair/highlights while fitting both chains.
    fitBoth();
  });

  // The widget may have finished loading while the Compare tab was hidden.
  // Fit it as soon as Bootstrap reveals its canvas.
  function onCompareShown() {
    needsFit = true;
    scheduleFit();
  }
  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab", function (event) {
      const link = event && event.target;
      if (link && link.getAttribute &&
          link.getAttribute("data-value") === "compare") onCompareShown();
    });
  }
  if (window.ResizeObserver) {
    const holder = document.querySelector(".ram-compare-ngl");
    if (holder) {
      let lastWidth = 0;
      new window.ResizeObserver(function () {
        const width = holder.clientWidth;
        if (width && !lastWidth) onCompareShown();
        lastWidth = width;
        if (width && stage && typeof stage.handleResize === "function")
          stage.handleResize();
      }).observe(holder);
    }
  }

  // Also support the first frame before an explicit readiness message arrives.
  document.addEventListener("shiny:value", function (event) {
    if (event.target && event.target.id === "NGLCompare") {
      needsFit = true;
      scheduleFit();
    }
  });
})(window, document);
