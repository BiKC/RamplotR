import assert from 'node:assert/strict';
import { waitForController } from '../scripts/shinylive-bootstrap.mjs';
class Workers extends EventTarget {
  controller = null;
  listeners = 0;
  addEventListener(type,callback,...args) { super.addEventListener(type,callback,...args); this.listeners++; }
  removeEventListener(type,callback,...args) { super.removeEventListener(type,callback,...args); this.listeners--; }
}
const workers = new Workers();
const ready = waitForController(workers,{timeoutMs:1000,expectedScriptURL:'worker.js'});
workers.controller={scriptURL:'other.js'};workers.dispatchEvent(new Event('controllerchange'));
workers.controller={scriptURL:'worker.js'};workers.dispatchEvent(new Event('controllerchange'));
assert.equal(await ready,workers.controller);assert.equal(workers.listeners,0);
assert.equal(await waitForController(workers),workers.controller);
const missing = new Workers();
await assert.rejects(waitForController(missing,{timeoutMs:5}),/did not give this page control/);
assert.equal(missing.listeners,0);
const abort = new AbortController();
const cancelled = waitForController(missing,{timeoutMs:1000,signal:abort.signal});
abort.abort();await assert.rejects(cancelled,/cancelled/);assert.equal(missing.listeners,0);
console.log('Delayed worker control, timeout and cancellation clean up without automatic reloads.');
