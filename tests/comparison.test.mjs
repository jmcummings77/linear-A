import assert from 'node:assert/strict';
import {test} from 'node:test';
import {normalize,compare,exportComparison} from '../benchmarks/comparison/model.mjs';
const raw=(times=[10,20,30])=>({schema_version:1,revision:'a'.repeat(40),source_sha256:'b'.repeat(64),dirty:false,seed:17,suite:'quick',machine:{os:'Test',architecture:'arm64',logical_cpus:8},methodology:{numeric_type:'double',matrix_layout:'row major',timing:'kernel',warmup:'5 calls',sampling:'serial',sample_unit:'fresh_process_per_batch',workload_version:'v1'},implementations:[{id:'c',name:'C',status:'passed',toolchain:'clang 1',build_commands:[['clang','-O3']]}],results:[{implementation:'c',operation:'multiply',size:16,status:'passed',median_ns:9999,samples:times.map(v=>({elapsed_ns:v*10,iterations:10,ns_per_op:v}))}]});
test('ratio orientation, recomputed median/MAD and observed range envelope',()=>{
 const a=normalize(raw()),b=normalize(raw([20,40,60])),r=compare(a,b).rows[0];
 assert.equal(a.rows[0].median,20);assert.equal(a.rows[0].mad,10);assert.equal(r.ratio,2);assert.equal(r.change,100);assert.deepEqual(r.envelope,[20/30,60/10]);
 assert.equal(compare(b,a).rows[0].ratio,.5);
});
test('missing and failed rows are never represented as zero timings',()=>{
 const a=raw(),b=raw();b.results[0].size=32;const r=compare(normalize(a),normalize(b));assert.equal(r.rows.length,2);assert.ok(r.rows.every(row=>row.ratio===null));
 b.results[0].size=16;b.results[0].status='failed';assert.equal(compare(normalize(a),normalize(b)).rows[0].ratio,null);
});
test('seed and timing-contract mismatches withhold ratios; toolchains warn',()=>{
 for(const edit of [d=>d.seed++,d=>d.methodology.timing='different',d=>d.methodology.workload_version='v2',d=>delete d.methodology.warmup]){
  const b=raw();edit(b);assert.equal(compare(normalize(raw()),normalize(b)).rows[0].ratio,null);
 }
 const b=raw();b.implementations[0].toolchain='clang 2';const r=compare(normalize(raw()),normalize(b));assert.equal(r.rows[0].ratio,1);assert.ok(r.rows[0].warnings.includes('Toolchain changed.'));
});
test('legacy provenance is visibly uncertain and verification failure blocks ratios',()=>{
 const b=raw();delete b.methodology.sample_unit;delete b.methodology.workload_version;b.dirty=true;
 const r=compare(normalize(raw()),normalize(b));assert.ok(r.warnings.some(w=>w.includes('Legacy')));assert.ok(r.warnings.some(w=>w.includes('independence')));assert.ok(r.warnings.some(w=>w.includes('dirty')));
 b.implementations[0].status='failed';assert.equal(compare(normalize(raw()),normalize(b)).rows[0].ratio,null);
});
test('malformed, duplicate, nonfinite or inconsistent measurements are rejected',()=>{
 for(const edit of [d=>d.results.push(d.results[0]),d=>d.results[0].samples[0].ns_per_op=NaN,d=>d.results[0].samples[0].iterations=0,d=>d.results[0].samples[0].elapsed_ns=1,d=>d.results[0].samples=[],d=>d.results[0].size=true,d=>d.results[0].status='/Users/private',d=>d.schema_version=2]){
  const d=raw();edit(d);assert.throws(()=>normalize(d));
 }
});
test('exports omit arbitrary metadata, command paths and environment values',()=>{
 const d=raw();d.env={TOKEN:'private'};d.implementations[0].toolchain='/Users/name/private/token';d.results[0].error='secret';
 const a=normalize(d),output=JSON.stringify(exportComparison(a,a,compare(a,a)));assert.ok(!output.includes('/Users/'));assert.ok(!output.includes('private'));assert.ok(!output.includes('secret'));assert.deepEqual(JSON.parse(output).rows[0].baseline.samples,a.rows[0].samples);
});
test('numeric overflow cannot silently become an exported null measurement',()=>{
 const a=normalize(raw([1e-300])),b=normalize(raw([1e300]));const r=compare(a,b).rows[0];assert.equal(r.ratio,null);assert.match(r.state,/finite range/);
});

test('even medians avoid overflow and tiny inconsistent timings are rejected',()=>{
 const d=raw([1e307,1e307]);d.results[0].samples.forEach(s=>{s.elapsed_ns=1e308;s.iterations=10;});
 assert.equal(normalize(d).rows[0].median,1e307);
 const bad=raw([1e-20]);bad.results[0].samples[0].ns_per_op=1e-11;assert.throws(()=>normalize(bad));
});
