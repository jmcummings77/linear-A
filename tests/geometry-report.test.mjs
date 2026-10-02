import assert from 'node:assert/strict';
import test from 'node:test';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const source=readFileSync(new URL('../benchmarks/geometry-report.mjs',import.meta.url),'utf8');
const template=readFileSync(new URL('../benchmarks/report.html',import.meta.url),'utf8');
const styles=['open','solid','chevron','tapered'];
const colors=['transition','cyan','blue','coral','violet','white'];
const motions=['static','flow','pulse'];
const widths=['1','1.5','2.25'];

class Element {
  constructor(tag){this.tagName=tag.toUpperCase();this.children=[];this.parentElement=null;this.dataset={};this.style={};this.listeners=new Map();this.attributes={};this._text='';this._value='';this.selectedIndex=-1;this.disabled=false;this.checked=false;this.hidden=false;this.className='';this.classList={toggle:(name,on)=>{const classes=new Set(this.className.split(' ').filter(Boolean));if(on)classes.add(name);else classes.delete(name);this.className=[...classes].join(' ');}};}
  set textContent(value){this.children=[];this._text=String(value);}
  get textContent(){return this._text+this.children.map(child=>child.textContent).join('');}
  set value(value){if(this.tagName==='SELECT')this.selectedIndex=this.options.findIndex(option=>option.value===String(value));else this._value=String(value);}
  get value(){return this.tagName==='SELECT'?(this.options[this.selectedIndex]?.value??''):this._value;}
  get options(){return this.children.filter(child=>child.tagName==='OPTION');}
  append(...children){for(const child of children){child.parentElement=this;this.children.push(child);if(this.tagName==='SELECT'&&this.selectedIndex<0)this.selectedIndex=0;}}
  replaceChildren(...children){this.children=[];this._text='';this.append(...children);}
  setAttribute(name,value){this.attributes[name]=String(value);}
  addEventListener(type,handler){if(!this.listeners.has(type))this.listeners.set(type,[]);this.listeners.get(type).push(handler);}
  dispatch(type,extra={}){const event={type,target:this,preventDefault(){},...extra};this['on'+type]?.(event);for(const callback of this.listeners.get(type)||[])callback(event);}
  querySelectorAll(tag){return this.children.filter(child=>child.tagName===tag.toUpperCase());}
  getBoundingClientRect(){return {width:720,height:430};}
  setPointerCapture(){}
  showModal(){this.open=true;}
  close(){this.open=false;}
}

/** Record drawing operations and reject invalid coordinates like a canvas would. */
class CanvasContext {
  constructor(){this.strokeStyle='#000';this.fillStyle='#000';this.lineWidth=1;this.globalAlpha=1;this.stack=[];this.clearRect();}
  finite(...values){assert.ok(values.every(Number.isFinite),`nonfinite canvas coordinates: ${values}`);}
  setTransform(...values){this.finite(...values);}
  clearRect(){this.strokes=[];this.fills=[];this.labels=[];this.path=[];}
  beginPath(){this.path=[];}
  moveTo(...values){this.finite(...values);this.path.push(['move',...values]);}
  lineTo(...values){this.finite(...values);this.path.push(['line',...values]);}
  closePath(){this.path.push(['close']);}
  arc(x,y,radius,start,end,anticlockwise=false){this.finite(x,y,radius,start,end);assert.ok(radius>=0);this.path.push(['arc',x,y,radius,start,end,anticlockwise]);}
  stroke(){this.finite(this.lineWidth,this.globalAlpha);assert.ok(this.lineWidth>0);this.strokes.push({path:structuredClone(this.path),color:this.strokeStyle,width:this.lineWidth,alpha:this.globalAlpha});}
  fill(){this.finite(this.globalAlpha);this.fills.push({path:structuredClone(this.path),color:this.fillStyle,alpha:this.globalAlpha});}
  fillText(text,x,y){this.finite(x,y);this.labels.push({text,x,y});}
  save(){this.stack.push({strokeStyle:this.strokeStyle,fillStyle:this.fillStyle,lineWidth:this.lineWidth,globalAlpha:this.globalAlpha,lineCap:this.lineCap,lineJoin:this.lineJoin});}
  restore(){assert.ok(this.stack.length,'unbalanced canvas restore');Object.assign(this,this.stack.pop());}
  setLineDash(values){this.finite(...values);}
  drawing(){return JSON.stringify({strokes:this.strokes,fills:this.fills});}
}

function sceneFixture(){
  const identity=[1,0,0,0,1,0,0,0,1],result=identity.map(value=>value*2);
  // The zero vector and a vector along the initial camera axis exercise arrow
  // heads whose projected length is zero, alongside ordinary nonzero vectors.
  const cameraAxis=[-Math.sin(-.68)*Math.cos(.43),Math.sin(.43),Math.cos(-.68)*Math.cos(.43)];
  const points=[[0,0,0],[.5,0,0],[0,.5,0],cameraAxis,[1e-8,0,0]];
  const steps=[2,4,6].map((sum,index)=>({row:0,col:0,k:index,left:1,right:2,product:2,sum,matrix:[sum,0,0,0,0,0,0,0,0],source_line:100+index,call_path:['m_multiply']}));
  return {operation:'multiply',m:identity,n:result,result,points,input_vectors:points,output_vectors:points.map(point=>point.map(value=>value*2)),rotation_frames:[],steps,trace_available:true,shared_checks:{passed:119,total:119},determinant_m:1,determinant_result:8,eigen_m:null};
}

function report({available=true,reducedMotion=false,scene=sceneFixture()}={}){
  const nodes=new Map(),listeners=new Map(),mediaListeners=[],frames=new Map(),timers=new Map(),workers=[];
  const ctx=new CanvasContext();let sequence=0,now=0;
  const attributes=text=>Object.fromEntries([...text.matchAll(/([\w-]+)(?:="([^"]*)")?/g)].map(match=>[match[1],match[2]??'']));
  for(const match of template.matchAll(/<([a-z][\w-]*)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){
    const [,tag,raw,id]=match,node=new Element(tag),attrs=attributes(raw);nodes.set(id,node);
    node.parentElement=new Element('label');node.disabled='disabled'in attrs;node.checked='checked'in attrs;node.hidden='hidden'in attrs;
    if('value'in attrs)node.value=attrs.value;
    if(tag==='select'){
      const inner=template.slice(match.index+match[0].length).split('</select>')[0];
      for(const optionMatch of inner.matchAll(/<option\b([^>]*)>([^<]*)<\/option>/g)){
        const option=new Element('option'),optionAttrs=attributes(optionMatch[1]);option.value=optionAttrs.value??optionMatch[2];option.textContent=optionMatch[2];option.disabled='disabled'in optionAttrs;node.append(option);if('selected'in optionAttrs)node.selectedIndex=node.options.length-1;
      }
    }
  }
  const get=id=>{assert.ok(nodes.has(id),`unknown DOM id ${id}`);return nodes.get(id);};
  for(const value of ['0','0','1']){const input=new Element('input');input.value=value;get('geometry-axis').append(input);}
  get('geometry-canvas').getContext=()=>ctx;
  const bundle={available,geometry_worker_source:'export {};',reason:'No geometry module in this report.'};
  get('live-data').textContent=JSON.stringify(bundle);get('data').textContent=JSON.stringify({profiles:[]});
  const document={hidden:false,getElementById:get,createElement:tag=>new Element(tag),createTextNode:text=>{const element=new Element('#text');element.textContent=text;return element;},addEventListener(type,callback){if(!listeners.has(type))listeners.set(type,[]);listeners.get(type).push(callback);}};
  const media={matches:reducedMotion,addEventListener(type,callback){assert.equal(type,'change');mediaListeners.push(callback);},addListener(callback){mediaListeners.push(callback);}};
  const matchMedia=query=>{assert.match(query,/prefers-reduced-motion/);return media;};
  const window={matchMedia,dispatchEvent(){}};
  class Worker {
    constructor(){workers.push(this);this.messages=[];this.terminated=false;}
    postMessage(message){this.messages.push(message);}
    terminate(){this.terminated=true;}
    reply(result){this.onmessage({data:{type:'geometry',result}});}
  }
  vm.runInNewContext(source,{document,window,matchMedia,Worker,WebAssembly:{},Blob:class {},URL:{createObjectURL:()=>`blob:geometry-${++sequence}`,revokeObjectURL(){}},Event:class {},ResizeObserver:class {observe(){}},devicePixelRatio:1,performance:{now:()=>0},requestAnimationFrame:callback=>{const id=++sequence;frames.set(id,callback);return id;},cancelAnimationFrame:id=>frames.delete(id),setTimeout:(callback,delay=0)=>{const id=++sequence;timers.set(id,{callback,due:now+delay});return id;},clearTimeout:id=>timers.delete(id)});
  if(available){assert.equal(workers.length,1,'initial computation should start once');workers[0].reply(scene);}
  const snapshot=()=>JSON.stringify({matrix:get('geometry-matrix').children.map(cell=>cell.textContent),m:get('geometry-m').children.map(input=>input.value),n:get('geometry-n').children.map(input=>input.value),progress:get('geometry-progress').value,step:get('geometry-step').value,mode:get('geometry-mode').value,term:get('geometry-term').textContent,volume:get('geometry-volume').textContent});
  return {get,ctx,workers,frames,scene,snapshot,timers,
    advance(milliseconds){const end=now+milliseconds;while(true){const next=[...timers].filter(([,timer])=>timer.due<=end).sort((a,b)=>a[1].due-b[1].due)[0];if(!next)break;now=next[1].due;timers.delete(next[0]);next[1].callback();}now=end;},
    choose(id,value){get(id).value=value;get(id).dispatch('change');},
    input(id,value){get(id).value=value;get(id).dispatch('input');},
    check(id,checked){get(id).checked=checked;get(id).dispatch('change');},
    click(id){get(id).dispatch('click');},
    frame(time){const pending=[...frames];frames.clear();for(const [,callback]of pending)callback(time);},
    hidden(value){document.hidden=value;for(const callback of listeners.get('visibilitychange')||[])callback();},
    reduce(value){media.matches=value;for(const callback of mediaListeners)callback({matches:value});},
  };
}

test('geometry controls start with the grid visible and static arrows without an idle animation loop',()=>{
  const page=report();assert.equal(page.get('geometry-show-grid').checked,true);
  for(const [id,values]of [['geometry-arrow-style',styles],['geometry-arrow-color',colors],['geometry-arrow-motion',motions],['geometry-arrow-width',widths]])assert.deepEqual(page.get(id).options.map(option=>option.value),values);
  assert.equal(page.get('geometry-arrow-motion').value,'static');assert.equal(page.frames.size,0);
  assert.ok(page.get('geometry-arrow-note').textContent.length>0);
});

test('rapid matrix and sample edits recompute once after 300 ms of quiet using the latest values',()=>{
  const page=report(),input=page.get('geometry-m').children[0];
  assert.doesNotMatch(template,/id="geometry-compute"/);
  input.value='2';input.dispatch('input');page.advance(200);
  input.value='3';input.dispatch('input');page.advance(200);
  page.input('geometry-grid-x','7');page.input('geometry-grid-y','3');page.input('geometry-grid-z','1');
  page.advance(299);assert.equal(page.workers.length,1);
  page.advance(1);assert.equal(page.workers.length,2);
  const config=page.workers[1].messages[0].config;
  assert.equal(config.m[0],3);assert.deepEqual(Array.from(config.grid),[7,3,1]);
  const scene=sceneFixture();scene.points=Array.from({length:21},()=>[0,0,0]);scene.input_vectors=scene.points;scene.output_vectors=scene.points;
  page.workers[1].reply(scene);
  assert.equal(page.get('geometry-vector-count').textContent,'21');
  assert.match(page.get('geometry-status').textContent,/21 field vectors/);
  page.advance(30000);assert.equal(page.workers.length,2);assert.equal(page.timers.size,0);
});

test('new inputs cancel obsolete workers and late replies cannot replace the new field',()=>{
  const page=report();page.input('geometry-grid-x','3');page.advance(300);
  const stale=page.workers[1];page.input('geometry-grid-x','4');assert.equal(stale.terminated,true);
  stale.reply(sceneFixture());assert.equal(page.get('geometry-play').disabled,true);
  page.advance(300);const current=page.workers[2];current.reply(sceneFixture());
  const state=page.snapshot(),status=page.get('geometry-status').textContent;
  stale.reply({...sceneFixture(),result:Array(9).fill(99)});
  stale.onmessage({data:{type:'error',message:'obsolete error'}});
  stale.onerror({preventDefault(){}});
  assert.equal(page.snapshot(),state);assert.equal(page.get('geometry-status').textContent,status);
  assert.equal(page.get('geometry-canvas').attributes['aria-busy'],'false');
});

test('every mathematical input schedules a field update while display controls do not',()=>{
  const page=report();
  const changes=[()=>page.choose('geometry-preset','reflection'),()=>page.click('geometry-randomize'),
    ()=>page.choose('geometry-operation','scale'),()=>page.input('geometry-scalar','2'),
    ()=>page.choose('geometry-operation','rotate'),()=>page.input('geometry-angle','45'),
    ()=>{const axis=page.get('geometry-axis').children[0];axis.value='1';axis.dispatch('input');},
    ()=>{const n=page.get('geometry-n').children[0];n.value='2';n.dispatch('input');}];
  for(const change of changes){const count=page.workers.length;change();page.advance(299);assert.equal(page.workers.length,count);page.advance(1);assert.equal(page.workers.length,count+1);page.workers.at(-1).reply(sceneFixture());}
  const count=page.workers.length;
  page.check('geometry-show-box',false);page.choose('geometry-arrow-color','blue');page.check('geometry-loop',true);
  page.advance(1000);assert.equal(page.workers.length,count);
});

test('invalid inputs explain the problem and automatically recover when corrected',()=>{
  const page=report(),input=page.get('geometry-m').children[0];
  input.value='';input.dispatch('input');page.advance(300);
  assert.equal(page.workers.length,1);assert.match(page.get('geometry-status').textContent,/Every matrix cell needs a number/);
  input.value='2';input.dispatch('input');page.input('geometry-grid-x','16');page.advance(300);
  assert.equal(page.workers.length,1);assert.match(page.get('geometry-status').textContent,/sample counts from 1 to 15/);
  page.input('geometry-grid-x','1');page.advance(300);assert.equal(page.workers.length,2);
  page.workers[1].reply(sceneFixture());assert.equal(page.get('geometry-play').disabled,false);
  const unavailable=report({available:false});unavailable.click('geometry-randomize');unavailable.advance(1000);
  assert.equal(unavailable.workers.length,0);assert.equal(unavailable.timers.size,0);
  assert.match(unavailable.get('geometry-status').textContent,/No geometry module/);
});

test('looping wraps transformation and calculation playback, then stops when disabled',()=>{
  for(const mode of ['morph','steps']){
    const page=report();page.choose('geometry-mode',mode);page.choose('geometry-duration','10');
    assert.equal(page.get('geometry-loop').checked,false);page.check('geometry-loop',true);
    page.click('geometry-play');page.frame(0);
    for(let time=100;time<=10100;time+=100)page.frame(time);
    assert.equal(page.frames.size,1);assert.ok(Number(page.get('geometry-progress').value)<50);
    if(mode==='steps')assert.equal(page.get('geometry-step').value,'0');
    page.click('geometry-play');assert.equal(page.frames.size,0);
    page.click('geometry-play');page.frame(10200);page.check('geometry-loop',false);
    for(let time=10300;time<=21000;time+=100)page.frame(time);
    assert.equal(page.get('geometry-progress').value,'1000');assert.equal(page.frames.size,0);
    page.check('geometry-loop',true);page.click('geometry-reset');
    assert.equal(page.get('geometry-loop').checked,true);assert.equal(page.frames.size,0);
  }
});

test('settings dialog contains overlays and appearance and preserves the calculation',()=>{
  const page=report(),state=page.snapshot(),dialog=page.get('geometry-settings');
  assert.ok(!dialog.open);page.click('geometry-settings-open');assert.equal(dialog.open,true);
  page.check('geometry-show-box',false);page.choose('geometry-arrow-style','solid');
  page.click('geometry-settings-close');assert.equal(dialog.open,false);assert.equal(page.snapshot(),state);
  page.click('geometry-settings-open');assert.equal(page.get('geometry-show-box').checked,false);
  assert.equal(page.get('geometry-arrow-style').value,'solid');
  const content=template.match(/<dialog\b[^>]*id="geometry-settings"[^>]*>([\s\S]*?)<\/dialog>/)[1];
  assert.match(content,/id="geometry-show-axes"/);assert.match(content,/id="geometry-arrow-style"/);
  assert.match(template,/<footer[\s\S]*© 2026 J\.M\. Cummings\.[\s\S]*<\/footer>/);
});

test('grid toggle changes only its eighteen strokes and preserves computation and scrub position',()=>{
  const page=report();page.input('geometry-progress','375');
  const before=page.snapshot(),withGrid=structuredClone(page.ctx.strokes),fills=structuredClone(page.ctx.fills);
  page.check('geometry-show-grid',false);
  assert.equal(page.ctx.strokes.length,withGrid.length-18);
  assert.deepEqual(page.ctx.fills,fills);assert.equal(page.snapshot(),before);assert.equal(page.workers.length,1);
  page.check('geometry-show-grid',true);assert.deepEqual(page.ctx.strokes,withGrid);assert.equal(page.snapshot(),before);
});

test('arrow appearance preserves matrix inputs, transformation progress, and captured calculation step',()=>{
  const page=report();page.input('geometry-progress','420');let before=page.snapshot();
  for(const [id,value]of [['geometry-arrow-style','solid'],['geometry-arrow-color','violet'],['geometry-arrow-width','2.25']]){page.choose(id,value);assert.equal(page.snapshot(),before);}
  page.choose('geometry-mode','steps');page.input('geometry-step','2');before=page.snapshot();
  for(const [id,value]of [['geometry-arrow-style','tapered'],['geometry-arrow-color','coral'],['geometry-arrow-width','1'],['geometry-arrow-motion','flow']]){page.choose(id,value);assert.equal(page.snapshot(),before);}
  page.frame(0);page.frame(250);assert.equal(page.snapshot(),before);assert.equal(page.workers.length,1);
  page.choose('geometry-arrow-motion','static');assert.equal(page.frames.size,0);
});

test('each arrow shape, palette, and width draws finite zero-length and projected-short vectors',()=>{
  const page=report();for(const id of ['geometry-show-grid','geometry-show-axes','geometry-show-box','geometry-show-ellipse'])page.check(id,false);
  page.input('geometry-progress','500');const original=JSON.stringify(page.scene),shapeDrawings=new Set(),paletteDrawings=new Set(),swatches=new Set(),widthDrawings=new Set();
  for(const style of styles){page.choose('geometry-arrow-color','cyan');page.choose('geometry-arrow-width','1.5');page.choose('geometry-arrow-style',style);shapeDrawings.add(page.ctx.drawing());for(const color of colors){page.choose('geometry-arrow-color',color);for(const width of widths){page.choose('geometry-arrow-width',width);assert.ok(page.ctx.strokes.length+page.ctx.fills.length>0);assert.equal(page.ctx.stack.length,0);}}}
  assert.equal(shapeDrawings.size,styles.length,'styles must visibly change arrow paths');
  for(const color of colors){page.choose('geometry-arrow-color',color);paletteDrawings.add(page.ctx.drawing());swatches.add(page.get('geometry-field-swatch').style.background);}
  assert.equal(paletteDrawings.size,colors.length);assert.equal(swatches.size,colors.length);
  for(const width of widths){page.choose('geometry-arrow-width',width);widthDrawings.add(page.ctx.drawing());}
  assert.equal(widthDrawings.size,widths.length);assert.equal(JSON.stringify(page.scene),original);assert.equal(page.frames.size,0);
});

test('flow and pulse animate the arrows while a paused transformation stays fixed',()=>{
  for(const motion of ['flow','pulse']){
    const page=report();page.input('geometry-progress','350');const before=page.snapshot();page.choose('geometry-arrow-motion',motion);
    assert.equal(page.frames.size,1);page.frame(0);const initial=page.ctx.drawing();page.frame(250);
    assert.notEqual(page.ctx.drawing(),initial);assert.equal(page.snapshot(),before);assert.equal(page.frames.size,1);
    page.choose('geometry-arrow-color','blue');page.choose('geometry-arrow-width','2.25');assert.equal(page.frames.size,1,'appearance updates must not fork animation loops');
    page.choose('geometry-arrow-motion','static');assert.equal(page.frames.size,0);assert.equal(page.snapshot(),before);
  }
});

test('arrow effects stop when hidden or invalidated, and unavailable geometry stays idle',()=>{
  const page=report();page.choose('geometry-arrow-motion','flow');assert.equal(page.frames.size,1);
  page.hidden(true);assert.equal(page.frames.size,0);page.hidden(false);page.choose('geometry-arrow-motion','pulse');assert.equal(page.frames.size,1);
  const input=page.get('geometry-m').children[0];input.value='2';input.dispatch('input');assert.equal(page.frames.size,0);assert.equal(page.get('geometry-play').disabled,true);
  const unavailable=report({available:false});for(const motion of motions){unavailable.choose('geometry-arrow-motion',motion);assert.equal(unavailable.frames.size,0);}
  assert.equal(unavailable.workers.length,0);assert.match(unavailable.get('geometry-status').textContent,/no geometry|unavailable/i);
});

test('reduced motion suppresses arrow effects but keeps explicit transformation playback usable',()=>{
  const page=report({reducedMotion:true});page.choose('geometry-arrow-motion','flow');assert.equal(page.frames.size,0);
  page.choose('geometry-arrow-motion','pulse');assert.equal(page.frames.size,0);const before=page.snapshot();
  page.click('geometry-play');assert.equal(page.frames.size,1);page.frame(0);page.frame(500);assert.notEqual(page.snapshot(),before);
  page.click('geometry-play');assert.equal(page.frames.size,0);
  page.reduce(false);assert.equal(page.frames.size,1);page.reduce(true);assert.equal(page.frames.size,0);
});

function axisLabels(page){
  const labels=page.ctx.labels.filter(label=>['x','y','z'].includes(label.text));
  assert.equal(labels.length,3,'coordinate labels must remain visible');
  return structuredClone(labels);
}

function assertSameCamera(actual,expected){
  assert.deepEqual(actual.map(label=>label.text),expected.map(label=>label.text));
  for(let index=0;index<actual.length;index++)for(const coordinate of ['x','y'])assert.ok(Math.abs(actual[index][coordinate]-expected[index][coordinate])<1e-7,`${actual[index].text} ${coordinate} changed from ${expected[index][coordinate]} to ${actual[index][coordinate]}`);
}

test('wheel and keyboard zoom crop past the old limit and share limits and reset behavior',()=>{
  const page=report(),canvas=page.get('geometry-canvas'),initial=axisLabels(page),state=page.snapshot();
  const magnified=factor=>initial.map(label=>({...label,x:360+(label.x-360)*factor,y:215+(label.y-215)*factor}));
  canvas.dispatch('wheel',{deltaY:-1e5});
  assertSameCamera(axisLabels(page),magnified(25));
  assert.ok(axisLabels(page).some(({x,y})=>x<0||x>720||y<0||y>430),'close-ups must allow cropping');
  canvas.dispatch('keydown',{key:'+'});assertSameCamera(axisLabels(page),magnified(25));
  assert.equal(page.snapshot(),state);
  page.click('geometry-reset');assertSameCamera(axisLabels(page),initial);
  for(let i=0;i<100;i++)canvas.dispatch('keydown',{key:'+'});
  assertSameCamera(axisLabels(page),magnified(25));
  canvas.dispatch('wheel',{deltaY:1e5});assertSameCamera(axisLabels(page),magnified(.45));
  canvas.dispatch('keydown',{key:'-'});assertSameCamera(axisLabels(page),magnified(.45));
  canvas.dispatch('keydown',{key:'+'});assertSameCamera(axisLabels(page),magnified(.45*1.1));
  page.click('geometry-reset');assertSameCamera(axisLabels(page),initial);
  assert.equal(page.workers.length,1);
});

test('hidden transformed overlays do not shrink the field in morph, rotation, or calculation views',()=>{
  const baseline=report();for(const id of ['geometry-show-axes','geometry-show-box','geometry-show-ellipse','geometry-show-grid'])baseline.check(id,false);
  baseline.choose('geometry-arrow-style','solid');
  const anchors=page=>page.ctx.strokes.filter(stroke=>stroke.color.startsWith('rgb(')).map(stroke=>stroke.path[0]);
  const expected=anchors(baseline);
  for(const mode of ['morph','rotation','steps']){
    const scene=sceneFixture();scene.m=scene.m.map(value=>value*10);scene.result=scene.result.map(value=>value*10);
    scene.input_vectors=scene.input_vectors.map(vector=>vector.map(value=>value*10));
    scene.output_vectors=scene.output_vectors.map(vector=>vector.map(value=>value*10));
    scene.steps=scene.steps.map(frame=>({...frame,matrix:frame.matrix.map(value=>value*10)}));
    if(mode==='rotation'){scene.rotation_extent=40;scene.rotation_frames=[scene.m,scene.result];}
    if(mode==='steps')scene.steps[0].matrix=scene.m;
    const page=report({scene});page.choose('geometry-arrow-style','solid');
    if(mode==='steps'){page.choose('geometry-mode','steps');page.input('geometry-step','1');}
    for(const id of ['geometry-show-axes','geometry-show-box','geometry-show-ellipse','geometry-show-grid'])page.check(id,false);
    assert.deepEqual(anchors(page),expected,`${mode}: hidden matrices must not change field projection`);
    const state=page.snapshot();
    for(const id of ['geometry-show-box','geometry-show-axes']){
      page.check(id,true);assert.notDeepEqual(anchors(page),expected,`${mode}: visible geometry should affect fitting`);
      page.check(id,false);assert.deepEqual(anchors(page),expected);assert.equal(page.snapshot(),state);
    }
  }
});

test('close-up orbit stays finite for every arrow style with overlays visible or hidden',()=>{
  for(const overlays of [true,false]){
    const page=report(),canvas=page.get('geometry-canvas');
    for(const id of ['geometry-show-axes','geometry-show-box','geometry-show-ellipse','geometry-show-grid'])page.check(id,overlays);
    page.input('geometry-progress','420');const state=page.snapshot();
    canvas.dispatch('wheel',{deltaY:-1e5});page.check('geometry-auto-orbit',true);
    for(const style of styles){
      page.choose('geometry-arrow-style',style);page.choose('geometry-arrow-motion','flow');
      for(let i=0;i<=300;i++)page.frame(i*100);
      canvas.dispatch('pointerdown',{clientX:100,clientY:100,pointerId:1});
      canvas.dispatch('pointermove',{clientX:400,clientY:300});canvas.dispatch('pointerup');
      assert.equal(page.ctx.stack.length,0);assert.equal(page.snapshot(),state);
    }
  }
});

test('auto-orbit starts off and completes a full camera turn in thirty active seconds without changing the matrix',()=>{
  const page=report();assert.equal(page.get('geometry-auto-orbit').checked,false);assert.equal(page.frames.size,0);
  assert.ok(page.get('geometry-orbit-note').textContent.length>0);
  page.input('geometry-progress','420');const state=page.snapshot(),initial=axisLabels(page),drawing=page.ctx.drawing();
  const quadrants=new Set(),recordQuadrant=()=>{const labels=axisLabels(page);quadrants.add(['x','z'].map(axis=>Math.sign(labels.find(label=>label.text===axis).x-360)).join(','));};
  recordQuadrant();page.check('geometry-auto-orbit',true);assert.equal(page.frames.size,1);page.frame(0);
  for(let time=100;time<=30000;time+=100){page.frame(time);assert.equal(page.frames.size,1);if(time%7500===0&&time<30000){recordQuadrant();assert.notEqual(page.ctx.drawing(),drawing);}}
  assert.equal(quadrants.size,4,'one turn must visit all four camera quadrants');
  assertSameCamera(axisLabels(page),initial);assert.equal(page.snapshot(),state);assert.equal(page.workers.length,1);
  page.check('geometry-auto-orbit',false);assert.equal(page.frames.size,0);
});

test('orbit preserves captured calculation steps and shares one frame loop with arrows and transformation playback',()=>{
  const page=report();page.choose('geometry-mode','steps');page.input('geometry-step','2');
  const state=page.snapshot(),camera=axisLabels(page);page.check('geometry-auto-orbit',true);page.frame(0);page.frame(100);
  assert.notDeepEqual(axisLabels(page),camera);assert.equal(page.snapshot(),state);
  page.choose('geometry-arrow-motion','flow');page.click('geometry-play');assert.equal(page.frames.size,1);
  page.frame(200);page.frame(300);assert.equal(page.frames.size,1);
  page.click('geometry-play');assert.equal(page.frames.size,1,'camera and arrows continue after playback pauses');
  page.check('geometry-auto-orbit',false);assert.equal(page.frames.size,1,'arrow effects retain the shared frame loop');
  page.choose('geometry-arrow-motion','static');assert.equal(page.frames.size,0);
});

test('orbit selection survives reset and recomputation, with no idle loop while the scene is unavailable',()=>{
  const page=report();page.check('geometry-auto-orbit',true);page.frame(0);page.frame(100);
  page.click('geometry-reset');assert.equal(page.get('geometry-auto-orbit').checked,true);assert.equal(page.frames.size,1);
  page.input('geometry-grid-x','3');page.advance(300);assert.equal(page.frames.size,0);assert.equal(page.get('geometry-auto-orbit').checked,true);
  assert.equal(page.workers.length,2);page.workers[1].reply(sceneFixture());assert.equal(page.frames.size,1);
  const input=page.get('geometry-m').children[0];input.value='2';input.dispatch('input');assert.equal(page.frames.size,0);assert.equal(page.get('geometry-auto-orbit').checked,true);
  const unavailable=report({available:false});unavailable.check('geometry-auto-orbit',true);assert.equal(unavailable.frames.size,0);assert.equal(unavailable.workers.length,0);
});

test('hidden and reduced-motion states suppress orbit without accumulating skipped camera time',()=>{
  const page=report();page.check('geometry-auto-orbit',true);page.frame(0);page.frame(100);const camera=axisLabels(page);
  page.hidden(true);assert.equal(page.frames.size,0);page.hidden(false);assert.equal(page.frames.size,1);page.frame(100000);assertSameCamera(axisLabels(page),camera);
  page.reduce(true);assert.equal(page.frames.size,0);page.click('geometry-play');assert.equal(page.frames.size,1);
  const state=page.snapshot();page.frame(100100);page.frame(100200);assert.notEqual(page.snapshot(),state);assertSameCamera(axisLabels(page),camera);
  page.click('geometry-play');assert.equal(page.frames.size,0);page.reduce(false);assert.equal(page.frames.size,1);page.frame(200000);assertSameCamera(axisLabels(page),camera);
  page.frame(200100);assert.notDeepEqual(axisLabels(page),camera);
  const initiallyReduced=report({reducedMotion:true});initiallyReduced.check('geometry-auto-orbit',true);assert.equal(initiallyReduced.frames.size,0);
});

test('manual drag pauses orbit and pointer release resumes from the dragged camera without a jump',()=>{
  for(const endEvent of ['pointerup','pointercancel','lostpointercapture']){
    const page=report(),canvas=page.get('geometry-canvas');page.check('geometry-auto-orbit',true);page.frame(0);page.frame(100);
    const state=page.snapshot(),before=axisLabels(page);canvas.dispatch('pointerdown',{clientX:100,clientY:100,pointerId:1});assert.equal(page.frames.size,0,endEvent);
    canvas.dispatch('pointermove',{clientX:140,clientY:115,pointerId:1});const dragged=axisLabels(page);assert.notDeepEqual(dragged,before);assert.equal(page.frames.size,0);
    canvas.dispatch(endEvent,{pointerId:1});assert.equal(page.frames.size,1,endEvent);page.frame(100000);assertSameCamera(axisLabels(page),dragged);
    page.frame(100100);assert.notDeepEqual(axisLabels(page),dragged);assert.equal(page.snapshot(),state);
  }
});

test('unchecking and rechecking orbit keeps the current view and resets its animation clock',()=>{
  const page=report();page.check('geometry-auto-orbit',true);page.frame(0);page.frame(100);const camera=axisLabels(page);
  page.check('geometry-auto-orbit',false);assert.equal(page.frames.size,0);assertSameCamera(axisLabels(page),camera);
  page.check('geometry-auto-orbit',true);assert.equal(page.frames.size,1);page.frame(100000);assertSameCamera(axisLabels(page),camera);
  page.frame(100100);assert.notDeepEqual(axisLabels(page),camera);
});
