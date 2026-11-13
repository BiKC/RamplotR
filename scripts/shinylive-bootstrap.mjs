// App-local bootstrap: wait for worker control before importing Shinylive.
// The upstream loader immediately reloads when ready resolves before control,
// which can repeat indefinitely in Gecko. Never reload automatically here.
export function waitForController(serviceWorkers, { timeoutMs = 15000, signal,
  expectedScriptURL } = {}) {
  const controlled = () => serviceWorkers.controller && (!expectedScriptURL ||
    serviceWorkers.controller.scriptURL === expectedScriptURL);
  if (controlled()) return Promise.resolve(serviceWorkers.controller);
  return new Promise((resolve, reject) => {
    let timer;
    const cleanup = () => {
      clearTimeout(timer);
      serviceWorkers.removeEventListener('controllerchange', changed);
      if (signal) signal.removeEventListener('abort', aborted);
    };
    const changed = () => {
      if (controlled()) { cleanup(); resolve(serviceWorkers.controller); }
    };
    const aborted = () => { cleanup(); reject(new Error('Startup cancelled.')); };
    serviceWorkers.addEventListener('controllerchange', changed);
    if (signal) signal.addEventListener('abort', aborted, { once: true });
    timer = setTimeout(() => {
      cleanup(); reject(new Error('The browser did not give this page control of its service worker.'));
    }, timeoutMs);
    if (signal && signal.aborted) aborted(); else changed();
  });
}

export async function start() {
  const root = document.getElementById('root');
  const status = document.createElement('p');
  status.setAttribute('role', 'status');
  status.textContent = 'Starting RamplotR…';
  root.replaceChildren(status);
  const abort = new AbortController();
  window.addEventListener('pagehide', () => abort.abort(), { once: true });
  try {
    if (!navigator.serviceWorker) throw new Error('Service workers are disabled or unavailable. Use HTTPS or localhost.');
    const url = new URL('../shinylive-sw.js', import.meta.url);
    await navigator.serviceWorker.register(url.href, { type: 'module' });
    await waitForController(navigator.serviceWorker, { signal: abort.signal, expectedScriptURL: url.href });
    if (abort.signal.aborted) return;
    const { runExportedApp } = await import('../shinylive/shinylive.js');
    if (abort.signal.aborted) return;
    await runExportedApp({ id: 'root', appEngine: 'r', relPath: '../' });
  } catch (error) {
    if (abort.signal.aborted) return;
    status.textContent = 'RamplotR could not start: ' + error.message;
    const retry = document.createElement('button');
    retry.type = 'button'; retry.textContent = 'Reload once';
    retry.addEventListener('click', () => window.location.reload());
    root.replaceChildren(status, retry);
  }
}

if (typeof window !== 'undefined' && typeof document !== 'undefined') start();
