import assert from 'node:assert/strict';
import {test} from 'node:test';
import {readFileSync} from 'node:fs';
import {createMatrixAPI} from '../ports/wasm/matrix.mjs';
import {exploreMultigrid} from '../benchmarks/multigrid-live.mjs';
const api=await createMatrixAPI();
test('live multigrid captures both legs and true level residuals',()=>{
 for(const width of [3,7,15,31])for(const pattern of ['mixed','smooth','rough','broad']){
  const r=exploreMultigrid(api,{width,pattern});assert.equal(r.frames.length,4*(r.levels-1)+2);
  const first=r.frames[0],last=r.frames.at(-1);assert.deepEqual(first.error,r.truth);assert.equal(last.phase,'post_smooth');assert.equal(last.level,0);
  for(const f of r.frames){assert.equal(f.x.length,f.width*f.width);assert.equal(f.b.length,f.x.length);assert.equal(f.error===null,f.level!==0);}
  for(const run of Object.values(r.runs)){assert.equal(run.residuals.length,run.iterations+1);assert.ok(run.residuals.at(-1)<=run.residuals[0]*1.00001e-9);}
  if(pattern==='broad'&&width>=7)assert.ok(r.runs.multigrid.iterations<r.runs.cg.iterations);
 }
});
test('live grid validation and report controls',()=>{
 assert.throws(()=>exploreMultigrid(api,{width:9}));assert.throws(()=>exploreMultigrid(api,{pattern:'invalid'}));
 const html=readFileSync('benchmarks/multigrid-report.html','utf8'),script=readFileSync('benchmarks/multigrid-report.mjs','utf8');
 for(const id of ['live-size','pattern','progress','loop','play','previous','next','comparison','report-theme'])assert.ok(html.includes(`id="${id}"`));
 assert.ok(script.includes('},300)'));assert.ok(script.includes('m.id!==request'));assert.ok(html.includes('Copyright © 2026 J.M. Cummings.'));
});
