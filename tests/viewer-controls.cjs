const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
(async function () {
const messages = {}, events = {}, visibility = {}, motion = {}, frames = [];
const status = { setAttribute() {} };
const progress = { setAttribute() {} };
let createdBusy = false;
const marker = { style:{}, isConnected:true, setAttribute() {} };
const checkboxes = { dna: {checked:false}, rna: {checked:false}, ligands: {checked:false} };
const layerBuilds = [];
const addedRepresentationHandlers = [];
const component = { structure: {}, signals: { representationAdded: { add(callback) {addedRepresentationHandlers.push(callback);} } }, addRepresentation(type,params) {
  layerBuilds.push({type,params});
  return { visible:true, setVisibility(value) { this.visible=value; }, setSelection(value) { this.selection=value; } };
} };
const sections = ['A', 'B'].map(chain => ({ getAttribute() { return chain; } }));
const plot = { data: [], style: {}, clientWidth: 600,
  _fullLayout:{ xaxis:{l2p:x=>x+180,_length:360,_offset:50}, yaxis:{l2p:y=>180-y,_length:360,_offset:40} },
  appendChild() {},
  classList: { contains() { return true; } }, on() {} };
let renders = 0, restyles = 0;
const stage = {
  getRepresentationsByName(name) {
    return { setVisibility(value) { visibility[name] = value; } };
  },
  tasks: { onZeroOnce(callback) { callback(); } },
  setRock(value) { motion.rock = value; }, setSpin(value) { motion.spin = value; }
};
const document = {
  getElementById(id) { return id === 'plotly' ? plot : checkboxes[id] || null; },
  querySelector() { return null; }, querySelectorAll(selector) { return selector === "[data-ram-chain]" ? sections : []; },
  createElement(tag) { if(tag === "progress") return progress; if(!createdBusy) { createdBusy=true;return status; } return marker; }, body: { appendChild() {} },
  addEventListener(name, callback) { (events[name] ||= []).push(callback); }
};
const jqueryEvents = {};
const jq = { on(name, callback) { jqueryEvents[name] = callback; return jq; } };
const window = {
  getNGLStructure() { return [component]; },
  getNGLStage() { return stage; }, jQuery() { return jq; },
  requestAnimationFrame(callback) { frames.push(callback); },
  Shiny: { addCustomMessageHandler(name, callback) { messages[name] = callback; } },
  Plotly: {
    react(element, traces) { renders++; element.data = traces; return Promise.resolve(); },
    restyle(element, update, indices) {
      restyles++;
      indices.forEach((index, i) => { element.data[index].visible = update.visible[i]; });
      return Promise.resolve();
    },
    relayout() {}
  }
};
const settle = async () => { for (let i = 0; i < 8; i++) await Promise.resolve(); };
const paint = async () => { while (frames.length) frames.shift()(); await settle(); };
vm.runInNewContext(fs.readFileSync('shinyRam/www/custom.js', 'utf8'), { window, document, console });
// Pick before the first server handshake and before representations exist.
const earlyPicker = { id:'chainselection', options:[{value:'A'},{value:'B'}],
  selectedOptions:[{value:'B'}], closest() {return {getAttribute() {return 'startup';}};} };
jqueryEvents['changed.bs.select']({target:earlyPicker});await paint();
assert.equal(visibility['ram-chain-A'],false);
const stageLookup=window.getNGLStage;
window.getNGLStage=()=>{throw new Error('Widget not bound yet');};
assert.doesNotThrow(()=>messages['ram-chain-visibility']({structure:'startup',chains:['A','B'],selected:['A','B']}));
window.getNGLStage=stageLookup;
await settle();assert.equal(visibility['ram-chain-A'],false,'first handshake preserves early selection');
let lateVisible=true;
addedRepresentationHandlers.forEach(callback=>callback({name:'ram-chain-A',setVisibility(value){lateVisible=value;}}));
assert.equal(lateVisible,false,'late initial representations use the current selection');
messages['ram-chain-visibility']({structure:'startup',chains:['A','B'],selected:['B']});await settle();
assert.equal(status.hidden,true);
messages['ram-chain-visibility']({ chains: ['A', 'B'], selected: ['A', 'B'] });
messages.process({ name: 'test', matrix: { x: [-180,180], y: [-180,180], z: [[0,1],[1,0]] },
  limits: [.1,.2,.3], df: process.argv[2] ? JSON.parse(fs.readFileSync(process.argv[2], 'utf8')) :
    { chain: ['A','B'], resn: ['ALA','GLY'], resi: [1,2], phi: [1,2], psi: [3,4] } });
await settle();
assert.equal(renders, 1);
const contourObjects = plot.data.slice(0,4);
const trace = plot.data.find(t=>t.customdata?.length);
const detail = trace.customdata[0]; const beforeSelection=restyles;
messages['ram-selection']({chain:detail[0],resi:detail[1],insertion_code:detail[2]});
assert.equal(restyles,beforeSelection,'residue marker never triggers a Plotly restyle');
assert.equal(marker.hidden,false);
assert.equal(marker.style.left,(trace.x[0]+230)+'px');
messages['ram-selection']({clear:true}); assert.equal(marker.hidden,true);
// The actual Bootstrap picker event must show feedback before expensive work.
jqueryEvents['changed.bs.select']({ target: { id: 'chainselection', selectedOptions: [{ value: 'B' }] } });
assert.equal(status.hidden, false);
assert.equal(status.textContent, 'Updating chains…');
assert.equal(progress.hidden, false);
assert.equal(visibility['ram-chain-A'], true, 'NGL work waits for feedback to paint');
await paint();
assert.equal(visibility['ram-chain-A'], false);
assert.equal(visibility['ram-chain-B'], true);
assert.equal(sections[0].hidden, true);
assert.equal(plot.data.find(t => t.meta?.ramChain === 'A').visible, false);
assert.equal(status.hidden, false, 'feedback remains until server acknowledges');
jqueryEvents['shiny:idle'](); assert.equal(status.hidden, false);
messages['ram-chain-visibility']({ chains: ['A', 'B'], selected: 'B' });
await settle(); assert.equal(status.hidden, true);
// Rapid selection changes: delayed acknowledgements must not undo the latest choice.
jqueryEvents['shiny:inputchanged']({ name: 'chainselection:shiny.picker', value: ['A'] });
jqueryEvents['shiny:inputchanged']({ name: 'chainselection', value: [] });
messages['ram-chain-visibility']({ chains: ['A', 'B'], selected: ['A'] });
await paint();
assert.equal(visibility['ram-chain-A'], false);
assert.equal(visibility['ram-chain-B'], false);
assert.equal(status.hidden, false);
messages['ram-chain-visibility']({ chains: ['A', 'B'], selected: [] });
await settle(); assert.equal(status.hidden, true);
// Returning to the acknowledged value can be deduplicated by Shiny: no reply arrives.
jqueryEvents['shiny:inputchanged']({ name: 'chainselection', value: ['A'] });
jqueryEvents['shiny:inputchanged']({ name: 'chainselection', value: [] });
await paint(); assert.equal(status.hidden, true);
messages['ram-chain-visibility']({ chains: ['A','B'], selected: ['A'] });
await settle(); assert.equal(visibility['ram-chain-A'], false);
assert.equal(renders, 1, 'chain changes never reconstruct the plot');
assert.ok(restyles > 0);
contourObjects.forEach((trace,i) => assert.equal(plot.data[i], trace));
for (const callback of events.change) callback({ target: { id: 'rocking', checked: false } });
assert.equal(motion.rock, false);
for (const callback of events.change) callback({ target: { id: 'spinning', checked: true } });
assert.deepEqual(motion, { rock: false, spin: true });
document.querySelector = selector => selector === '#regions table' ? {} : null;
for(const callback of events.click) callback({ target: { closest(selector) {
  return selector === 'a[data-value="residues"]' ? {} : null;
} } });
assert.equal(status.textContent,'Loading residue list…'); assert.equal(progress.hidden,false);
jqueryEvents['draw.dt init.dt']({target:{closest() {return {};}}});
assert.equal(status.hidden,true);
document.querySelector=()=>null;
messages['ram-bind-ngl']({ large:true, modelIndex:1, multipleModels:true, layers:{} });
for (const checked of [true,false,true]) {
  checkboxes.dna.checked=checked;
  for(const callback of events.change) callback({target:{id:'dna',checked}});
  assert.equal(status.hidden,false);
  await paint(); assert.equal(status.hidden,true);
}
assert.equal(layerBuilds.length,1,'DNA buffers survive toggling off/on');
assert.equal(layerBuilds[0].params.quality,'low');
assert.equal(layerBuilds[0].params.sele,'dna and /0');
for(const callback of events.change) callback({target:{name:'nglRepresentation',value:'ribbon',checked:true}});
assert.equal(status.hidden,false);
messages['ram-style-complete']({style:'ribbon'});await paint();assert.equal(status.hidden,true);
messages['ram-layer-model']({modelIndex:2,multipleModels:true});await paint();
assert.equal(layerBuilds.length,1,'model change reuses the representation');
jqueryEvents['shiny:busy'](); assert.equal(status.hidden, false);
jqueryEvents['shiny:idle'](); assert.equal(status.hidden, true);
console.log('Startup/settled chain selection, late representations, feedback, plot reuse, cached layers and animations passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
