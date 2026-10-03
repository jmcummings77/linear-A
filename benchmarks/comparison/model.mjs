/** Descriptive comparisons of saved batches. No significance tests or performance gates. */
const finite=v=>typeof v==='number'&&Number.isFinite(v);
const text=v=>typeof v==='string'?v.slice(0,2000):'';
export const median=values=>{const a=[...values].sort((x,y)=>x-y),m=Math.floor(a.length/2);return a.length%2?a[m]:a[m-1]+(a[m]-a[m-1])/2;};
export function normalize(raw,label='Snapshot'){
 if(!raw||raw.schema_version!==1||!Array.isArray(raw.results)||raw.results.length>20000||!Array.isArray(raw.implementations))throw new Error('Expected a schema-version-1 benchmark report.');
 const ids=new Set(),implementations=raw.implementations.map(i=>{
  if(!i||typeof i.id!=='string'||!/^[a-z0-9_-]{1,40}$/.test(i.id)||ids.has(i.id))throw new Error('Invalid or duplicate implementation.');ids.add(i.id);
  return {id:i.id,name:text(i.name)||i.id,status:text(i.status),toolchain:text(i.toolchain),build_commands:Array.isArray(i.build_commands)?i.build_commands.filter(Array.isArray).map(c=>c.map(text)):[]};
 });
 const keys=new Set();const rows=raw.results.map(r=>{
  if(!r||!ids.has(r.implementation)||typeof r.operation!=='string'||!/^[a-z0-9_-]{1,80}$/.test(r.operation)||!Number.isSafeInteger(r.size)||r.size<0)throw new Error('Invalid workload identity.');
  const key=JSON.stringify([r.implementation,r.operation,r.size]);if(keys.has(key))throw new Error('Duplicate workload identity.');keys.add(key);
  if(!['passed','failed','skipped','unavailable'].includes(r.status))throw new Error('Invalid row status.');
  const row={implementation:r.implementation,operation:r.operation,size:r.size,status:text(r.status),samples:[]};
  if(row.status!=='passed')return row;
  if(!Array.isArray(r.samples)||r.samples.length<1||r.samples.length>10000)throw new Error('A passed row needs recorded samples.');
  row.samples=r.samples.map(s=>{
   if(!s||!finite(s.elapsed_ns)||s.elapsed_ns<=0||!Number.isSafeInteger(s.iterations)||s.iterations<=0||!finite(s.ns_per_op)||s.ns_per_op<=0)throw new Error('Invalid sample timing.');
   const derived=s.elapsed_ns/s.iterations;if(Math.abs(derived-s.ns_per_op)>Math.max(Number.MIN_VALUE,derived*1e-9))throw new Error('Sample timing disagrees with elapsed time / iterations.');
   return {elapsed_ns:s.elapsed_ns,iterations:s.iterations,ns_per_op:s.ns_per_op};
  });
  const values=row.samples.map(s=>s.ns_per_op);row.median=median(values);row.mad=median(values.map(v=>Math.abs(v-row.median)));row.min=Math.min(...values);row.max=Math.max(...values);return row;
 });
 const machine={};for(const k of ['os','release','architecture','logical_cpus'])if(raw.machine?.[k]!==undefined)machine[k]=typeof raw.machine[k]==='number'?raw.machine[k]:text(raw.machine[k]);
 const methodology={};for(const k of ['numeric_type','matrix_layout','timing','warmup','sampling','sample_unit','workload_version'])methodology[k]=text(raw.methodology?.[k]);
 return {label:text(label),revision:/^[a-f0-9]{40,64}$/.test(raw.revision)?raw.revision:null,source_sha256:/^[a-f0-9]{64}$/.test(raw.source_sha256)?raw.source_sha256:null,dirty:typeof raw.dirty==='boolean'?raw.dirty:null,created_at:text(raw.created_at),suite:text(raw.suite),seed:Number.isSafeInteger(raw.seed)?raw.seed:null,machine,methodology,implementations,rows};
}
const same=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
export function compare(base,candidate){
 const warnings=[],blocked=[];
 if(base.seed===null||candidate.seed===null)blocked.push('Workload seed is unknown.');else if(base.seed!==candidate.seed)blocked.push('Workload seeds differ.');
 for(const field of ['numeric_type','matrix_layout','timing','warmup']){
  if(!base.methodology[field]||!candidate.methodology[field])blocked.push(`${field} metadata is missing.`);
  else if(base.methodology[field]!==candidate.methodology[field])blocked.push(`${field} differs.`);
 }
 const bv=base.methodology.workload_version,cv=candidate.methodology.workload_version;
 if(bv&&cv&&bv!==cv)blocked.push('Workload versions differ.');else if(!bv||!cv)warnings.push('Legacy workload version is unspecified; inspect source/protocol changes before attributing a difference to code.');
 if(!same(base.machine,candidate.machine))warnings.push('Recorded hardware or operating-system details differ.');
 warnings.push('Basic hardware metadata cannot establish identical machines or controlled thermal/scheduling conditions.');
 if(base.suite!==candidate.suite)warnings.push('Suite/calibration settings differ.');
 if(base.methodology.sampling!==candidate.methodology.sampling)warnings.push('Sampling descriptions differ.');
 if(!base.methodology.sample_unit||!candidate.methodology.sample_unit)warnings.push('Process independence is unspecified in these saved reports; rows count recorded timing batches.');
 if(base.dirty!==false||candidate.dirty!==false)warnings.push('At least one source state is dirty or unspecified; a commit alone does not identify all measured code.');
 if(!base.revision||!candidate.revision)warnings.push('Recorded Git revision is missing.');
 else if(base.revision===candidate.revision)warnings.push('Both snapshots record the same Git revision; this is not evidence of a commit-to-commit change.');
 const key=r=>JSON.stringify([r.implementation,r.operation,r.size]),a=new Map(base.rows.map(r=>[key(r),r])),b=new Map(candidate.rows.map(r=>[key(r),r]));
 const rows=[...new Set([...a.keys(),...b.keys()])].map(k=>{
  const before=a.get(k),after=b.get(k),identity=before||after,row={implementation:identity.implementation,operation:identity.operation,size:identity.size,before,after,warnings:[],ratio:null,change:null,envelope:null};
  if(!before||!after){row.state=!before?'missing baseline':'missing candidate';return row;}
  if(before.status!=='passed'||after.status!=='passed'){row.state=`baseline ${before.status}; candidate ${after.status}`;return row;}
  if(blocked.length){row.state='incompatible metadata';return row;}
  const ai=base.implementations.find(i=>i.id===row.implementation),bi=candidate.implementations.find(i=>i.id===row.implementation);
  if(ai.status!=='passed'||bi.status!=='passed'){row.state='implementation verification not passed';return row;}
  if(!ai.toolchain||!bi.toolchain)row.warnings.push('Toolchain version missing.');else if(ai.toolchain!==bi.toolchain)row.warnings.push('Toolchain changed.');
  if(!same(ai.build_commands,bi.build_commands))row.warnings.push('Build commands differ.');
  if(Math.min(before.samples.length,after.samples.length)<5)row.warnings.push('Fewer than five recorded batches on at least one side.');
  row.ratio=after.median/before.median;row.change=100*(row.ratio-1);
  row.envelope=[after.min/before.max,after.max/before.min];
  if(![row.ratio,row.change,...row.envelope].every(finite)){row.ratio=null;row.change=null;row.envelope=null;row.state='ratio outside finite range';return row;}
  row.state='descriptive comparison';return row;
 });
 return {warnings,blocked,rows};
}
// Exports contain numerical results and validated hashes only, never arbitrary input metadata.
export function exportComparison(base,candidate,result){
 return {schema_version:1,kind:'descriptive_snapshot_comparison',baseline:{revision:base.revision,source_sha256:base.source_sha256,dirty:base.dirty},candidate:{revision:candidate.revision,source_sha256:candidate.source_sha256,dirty:candidate.dirty},
  warnings:result.warnings,blocked:result.blocked,rows:result.rows.map(r=>({implementation:r.implementation,operation:r.operation,size:r.size,state:r.state,warnings:r.warnings,ratio:r.ratio,change_percent:r.change,observed_range_envelope:r.envelope,
   baseline:r.before?{status:r.before.status,samples:r.before.samples}:null,candidate:r.after?{status:r.after.status,samples:r.after.samples}:null}))};
}
