import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/amd-report.mjs',import.meta.url),'utf8');
function setup({width=700,results}={}){
 const elements=new Map(),timers=new Map(),workers=[],frames=new Map(),events={};let id=0;
 const draws=[],transforms=[];const ctx=new Proxy({},{get:(o,k)=>o[k]??(k==='fillText'?((...args)=>draws.push({args,font:o.font})):k==='setTransform'?((...args)=>transforms.push(args)):(()=>{}))});
 const element=(...args)=>({tagName:args.at(-1),clientWidth:width,clientHeight:width,value:'',checked:false,disabled:false,textContent:'',children:[],width:700,height:700,handlers:{},addEventListener(k,f){this.handlers[k]=f;},append(...v){this.children.push(...v);},replaceChildren(...v){this.children=v;},setAttribute(k,v){this[k]=v;},getContext:()=>ctx});
 const get=k=>{if(!elements.has(k))elements.set(k,element());return elements.get(k);};
 const recorded={implementations:[{id:'c',name:'C',status:'passed'}],created_at:'2026-10-03',revision:'abcdef',machine:{os:'test',architecture:'arm64'},methodology:{},results:[{implementation:'c',problem:'grid',size:12,operation:'total',ordering:'natural',status:'passed',median_ns:100,mad_ns:1},{implementation:'c',problem:'grid',size:12,operation:'total',ordering:'rcm',status:'failed',error:'zero pivot'}]};
 if(results)recorded.results=results;
 get('data').textContent=JSON.stringify(recorded);get('live').textContent=JSON.stringify({available:true,capabilities:['amd']});
 for(const [k,v] of Object.entries({language:'c',problem:'grid',size:'12',metric:'total','live-size':'12','live-problem':'scrambled',contrast:'2',speed:'4',progress:'0','graph-order':'amd'}))get(k).value=v;
 class Worker{constructor(){workers.push(this);this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:element,createElementNS:element,body:{},documentElement:{}},window:{devicePixelRatio:2,addEventListener:(k,f)=>events[k]=f},Worker,Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},MutationObserver:class{observe(){}},getComputedStyle:()=>({color:'#123'}),performance:{now:()=>0},requestAnimationFrame:f=>{frames.set(++id,f);return id;},cancelAnimationFrame:i=>frames.delete(i),setTimeout:(f,delay)=>{timers.set(++id,{f,delay});return id;},clearTimeout:i=>timers.delete(i)};
 vm.runInNewContext(source,sandbox);
 const flush=delay=>{for(const [i,t] of [...timers])if(t.delay===delay){timers.delete(i);t.f();}};
 const reply=(id)=>{const solver={ok:true,order:[0,1],trace:{edges:[[0,1,-1]],frames:[{pivot:0,neighbors:[1],fill:[]},{pivot:1,neighbors:[],fill:[]}]},points:[[0,0,-1,2],[1,0,0,1],[1,1,-1,2]],nnz:3,fillCount:1,residual:1e-14,reconstructionError:1e-14,symbolicMilliseconds:0,factorMilliseconds:1,solveMilliseconds:1,totalMilliseconds:2};workers[0].onmessage({data:{type:'result',id,result:{n:2,points:[[0,0,1],[1,1,2]],p:[1,0],inverse:[1,0],width:0,newWidth:0,natural:solver,rcm:solver,amd:solver}}});};
 return {get,workers,timers,frames,events,flush,reply,draws,transforms};
}
test('AMD controls debounce, ignore stale replies and retain recorded failures',()=>{
 const t=setup();assert.ok(t.get('saved-values').children.some(row=>row.children.some(cell=>String(cell.textContent).includes('zero pivot'))));
 t.workers[0].onmessage({data:{type:'ready'}});t.get('live-size').value='8';t.get('live-size').handlers.change();t.get('live-size').value='16';t.get('live-size').handlers.change();
 assert.equal([...t.timers.values()].filter(x=>x.delay===300).length,1);t.flush(300);assert.equal(t.workers[0].sent[1].config.size,16);
 t.reply(1);assert.equal(t.get('play').disabled,true);t.reply(3);assert.equal(t.get('play').disabled,false);assert.match(t.get('live-state').textContent,/All three factorizations finished/);
});
test('animation loops, scrubbing stops playback, teardown clears timers and workers',()=>{
 const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);
 t.get('loop').checked=true;t.get('play').handlers.click();let [id,f]=[...t.frames][0];t.frames.delete(id);f(11000);assert.equal(t.get('progress').value,'1000');assert.equal(t.frames.size,1);
 t.get('progress').value='300';t.get('progress').handlers.input();assert.equal(t.frames.size,0);assert.equal(t.get('play').textContent,'Animate elimination');
 t.events.pagehide();assert.ok(t.workers[0].terminated);assert.equal(t.timers.size,0);
});
test('worker timeout is visible and terminates work',()=>{const t=setup();t.flush(30000);assert.ok(t.workers[0].terminated);assert.match(t.get('live-state').textContent,/initialization timed out/);});
test('HTML IDs are unique, including executable data and live section',()=>{
 const html=readFileSync(new URL('../benchmarks/amd-report.html',import.meta.url),'utf8'),ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);assert.equal(new Set(ids).size,ids.length);
});

test('pivot buttons step the actual elimination history',()=>{const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);t.get('next').handlers.click();assert.equal(t.get('progress').value,'500');assert.match(t.get('frame-code').textContent,/eliminate original vertex 0/);t.get('previous').handlers.click();assert.equal(t.get('progress').value,'0');});

test('graph ordering and past-edge controls redraw without recomputing',()=>{const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);const count=t.workers[0].sent.length;t.get('graph-order').value='rcm';t.get('graph-order').handlers.change();assert.match(t.get('graph-note').textContent,/RCM/);t.get('past-edges').checked=true;t.get('past-edges').handlers.change();assert.equal(t.workers[0].sent.length,count);});

test('natural ordering is explicitly unnecessary rather than a failed timing',()=>{const t=setup();t.get('metric').value='ordering';t.get('metric').handlers.change();assert.ok(t.get('saved-chart').children.some(e=>e.textContent==='No ordering required'));});

function successfulRows(){
 return ['natural','rcm','amd'].flatMap((ordering,i)=>['total','solve','setup','symbolic','factor','ordering'].map(operation=>({implementation:'c',problem:'grid',size:12,operation,ordering,status:'passed',median_ns:100/(i+1),mad_ns:1,solver_iterations:8-i,bandwidth:10-i,fill_count:4-i,factor_nnz:20-i,logical_factor_bytes:480-i*24})));
}
test('chart explanation follows the selected metric while measured total stays distinct',()=>{
 const t=setup({results:successfulRows()}),total=t.get('total-tradeoff').textContent;
 assert.match(total,/Measured total tradeoff/);assert.match(total,/below 1/);
 for(const [metric,label] of [['fill_count','Fill entries'],['logical_factor_bytes','Logical factor bytes'],['factor','Numerical factorization time']]){t.get('metric').value=metric;t.get('metric').handlers.change();assert.ok(t.get('saved-summary').textContent.startsWith(label));assert.ok(t.get('saved-chart')['aria-label'].startsWith(label));assert.equal(t.get('total-tradeoff').textContent,total);}
});
test('narrow saved charts use readable CSS-pixel labels above proportional bars',()=>{
 const t=setup({width:302,results:successfulRows()}),svg=t.get('saved-chart'),texts=svg.children.filter(e=>e.tagName==='text'),bars=svg.children.filter(e=>e.tagName==='rect');
 assert.equal(Number(svg.viewBox.split(' ')[2]),302);assert.ok(Number(svg.height)<230);
 assert.ok(texts.every(e=>e['font-size']==='13'&&Number(e.x)<=302));
 assert.ok(Number(bars[0].y)>Number(texts[0].y));assert.equal(Number(bars[0].width),302);assert.equal(Number(bars[1].width),151);
 svg.clientWidth=800;t.events.resize();assert.equal(Number(svg.viewBox.split(' ')[2]),800);assert.ok(svg.children.filter(e=>e.tagName==='text').every(e=>e['font-size']==='13'));
});
test('missing metrics and a zero total remain finite and explicit',()=>{
 const rows=successfulRows().map(r=>({...r,median_ns:0,fill_count:undefined,bandwidth:undefined}));
 const t=setup({width:302,results:rows});assert.match(t.get('total-tradeoff').textContent,/ratio is not defined/);
 t.get('metric').value='fill_count';t.get('metric').handlers.change();assert.match(t.get('saved-summary').textContent,/unavailable/);assert.equal(t.get('saved-chart').children.filter(e=>e.tagName==='rect').length,0);
});
test('live plots use device-scaled backing pixels and readable CSS font sizes',()=>{
 const t=setup({width:302});t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);
 assert.equal(t.get('pattern').width,604);assert.equal(t.get('pattern').height,604);assert.ok(t.transforms.some(args=>args[0]===2&&args[3]===2));
 assert.ok(t.draws.length>0);assert.ok(t.draws.every(draw=>parseFloat(draw.font)>=12&&parseFloat(draw.font)<=14));
 t.get('pattern').clientWidth=420;t.events.resize();assert.equal(t.get('pattern').width,840);
});
test('report identity, sibling routes and section links support direct navigation',()=>{
 const html=readFileSync(new URL('../benchmarks/amd-report.html',import.meta.url),'utf8'),ids=new Set([...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]));
 assert.ok(html.includes('<h1>Approximate minimum degree (AMD)</h1>'));assert.match(html,/class="report-topbar"/);assert.match(html,/class="report-nav"/);
 for(const [,anchor] of html.matchAll(/href="#([^"]+)"/g))assert.ok(ids.has(anchor),`Missing target ${anchor}`);
 for(const sibling of ['ordering','cholesky','amd','ic0'].filter(name=>name!=='amd'))assert.ok(html.includes(`href="../${sibling}/"`));
 assert.ok(html.indexOf('id="play"')<html.indexOf('<canvas'));assert.doesNotMatch(html,/id="(?:frame-code|bandwidth-note)" aria-live="polite"/);
});

test('mobile pattern selection changes the visible order without resetting playback or solving again',()=>{
 const t=setup({width:302});t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);
 const requests=t.workers[0].sent.length;t.get('progress').value='500';t.get('pattern-view').value='rcm';t.get('pattern-view').handlers.change();
 assert.equal(t.get('pattern-comparison')['data-view'],'rcm');assert.equal(t.get('progress').value,'500');assert.equal(t.workers[0].sent.length,requests);
});
