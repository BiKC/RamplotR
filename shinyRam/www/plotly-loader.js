/* The first visit should show the input form without downloading Plotly.
 * Shiny calls this when a plot is actually needed. In-flight requests share
 * one promise; an unsuccessful download can be retried.
 */
(function (window, document) {
  "use strict";
  let loading = null;
  const url = "https://cdn.plot.ly/plotly-2.14.0.min.js";

  window.ramLoadPlotly = function () {
    if (window.Plotly) return Promise.resolve(window.Plotly);
    if (loading) return loading;

    loading = new Promise(function (resolve, reject) {
      const script = document.createElement("script");
      script.src = url;
      script.async = true;
      script.onload = function () {
        if (window.Plotly) resolve(window.Plotly);
        else reject(new Error("Plotly loaded without its global API"));
      };
      script.onerror = function () {
        if (script.parentNode) script.parentNode.removeChild(script);
        reject(new Error("Could not download Plotly. Check your internet connection."));
      };
      document.head.appendChild(script);
    }).catch(function (error) {
      loading = null;
      throw error;
    });
    return loading;
  };

  // Start the download while Shiny parses the selected structure. The form
  // remains fast to display and the plotting code keeps its normal fallback
  // for keyboard-triggered or programmatic submissions.
  document.addEventListener("click", function (event) {
    const control = event.target && typeof event.target.closest === "function"
      ? event.target.closest("button") : null;
    if (control && (control.id === "submit" ||
                    control.id === "compareSubmit")) {
      window.ramLoadPlotly().catch(function () {
        // The plot itself displays the actionable error if loading fails.
      });
    }
  }, true);
})(window, document);
