// Node-only linked comparison checks; no R, network or real WebGL required.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

const handlers = {};
const events = {};
const moves = [];
const highlights = [];
const sent = [];
let clicked = null;
let holder = {clientWidth: 630, clientHeight: 410};
const newBox = (lo,hi) => ({
  lo, hi,
  isEmpty() {return this.lo > this.hi;},
  clone() {return newBox(this.lo,this.hi);},
  union(other) {this.lo=Math.min(this.lo,other.lo);
                this.hi=Math.max(this.hi,other.hi);return this;},
  getCenter(target) {target.x=(this.lo+this.hi)/2;return target;}
});
function component(offset) {
  return {
    structure: {
      getBoundingBox(selection) {
        assert.ok(selection.value.includes(" and protein"));
        const pair = /104:A|119:B/.test(selection.value);
        return newBox(offset + (pair?4:0), offset + (pair?6:12));
      }
    }
  };
}
const components = [component(0),component(3)];
const stage = {
  compList: components,
  signals: {clicked: {add(fn) {clicked=fn;},remove(fn) {
    if (clicked===fn) clicked=null;
  }}},
  handleResize() {},
  getCenter() {return {x:0};},
  getZoomForBox(box) {return -10-(box.hi-box.lo);},
  animationControls: {zoomMove(center,zoom,duration) {
    moves.push({center:center.x,zoom,duration});
  }},
  getRepresentationsByName(name) {
    return {setSelection(sele) {highlights.push({name,sele});}};
  }
};
const window = {
  Shiny: {
    addCustomMessageHandler(name,callback) {
      assert.equal(callback.length,1,name+" should accept an argument");
      handlers[name]=callback;
    },
    setInputValue(name,value) {sent.push({name,value});}
  },
  getNGLStage() {return stage;},
  getNGLStructure() {return components;},
  NGL: {Selection: class Selection {constructor(value){this.value=value;}}},
  requestAnimationFrame(callback) {callback();}
};
const document = {
  addEventListener(name,callback) {events[name]=callback;},
  getElementById(id) {return id==="NGLCompare"?holder:null;},
  querySelector() {return null;}
};
vm.runInNewContext(fs.readFileSync("shinyRam/www/compare.js","utf8"),
  {window,document,Number,Array});
handlers["ram-compare-config"]({
  chainA:"A",chainB:"B",modelA:1,modelB:1,
  multipleA:false,multipleB:false
});
assert.ok(moves.length, "Viewer should fit the two visible chains");
assert.equal(moves.at(-1).center,7.5,
  "Combined viewport must use both selected chains, not the full structures");
handlers["ram-compare-pair"]({
  a:{chain:"A",resi:104,insertion_code:"",modelIndex:1,multipleModels:false},
  b:{chain:"B",resi:119,insertion_code:"",modelIndex:1,multipleModels:false}
});
assert.equal(highlights.at(-2).sele,"104:A and protein");
assert.equal(highlights.at(-1).sele,"119:B and protein");
assert.equal(moves.at(-1).center,6.5,
  "Pair selection should center over corresponding highlighted residues");
assert.equal(typeof clicked,"function","NGL atom picking must be linked");
clicked({component:components[1],atom:{
  chainname:"B",resno:119,inscode:"A"
}});
assert.equal(sent.at(-1).name,"ramCompareNglPick");
assert.equal(sent.at(-1).value.side,"b");
assert.equal(sent.at(-1).value.insertion_code,"A");
const before=moves.length;
holder.clientWidth=0;
handlers["ram-compare-pair"]({a:{
  chain:"A",resi:104,insertion_code:"",modelIndex:1,multipleModels:false
},b:null});
assert.equal(moves.length,before,"Do not frame an invisible viewer");
holder.clientWidth=630;
handlers["ram-compare-ready"]({});
assert.ok(moves.length>before,
  "Once the Compare canvas becomes visible, focus the pending selection");
events.click({target:{closest(){return {id:"compareResetView"};}}});
assert.equal(moves.at(-1).center,7.5,
  "Fit both chains button must reset the comparison camera");
assert.ok(highlights.at(-1).sele==="none",
  "Alignment gaps should clear a missing partner's highlight");
console.log("Linked comparison framing and paired NGL picks passed.");
