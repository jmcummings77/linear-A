#!/usr/bin/env python3
"""IC(0) setup, application and repeated CG solves under three symmetric orderings."""
import argparse,datetime,json,math,os,platform,re,time
from pathlib import Path
import run,sparse
from sparse_reference import diffusion,multiply
from ordering_reference import rcm,permute
from amd_reference import amd
from ic0_reference import factor,apply,fixtures
from cholesky_reference import analyze
from publication import PublicSanitizer
from report import json_for_html
from report_design import apply_report_design
HERE=Path(__file__).resolve().parent

def cases(size,problem='grid'):
 a=diffusion(size,2);n=a['rows'];truth=[1+(i%7)/10 for i in range(n)];b=multiply(a,truth)
 for ordering in ('natural','rcm','amd'):
  p=list(range(n)) if ordering=='natural' else rcm(a) if ordering=='rcm' else amd(a)
  q=permute(a,p);rhs=[b[i] for i in p];x=[truth[i] for i in p];rp,ci,v=factor(q,exact=False)
  stats=dict(factor_nnz=len(v),fill_count=0,logical_factor_bytes=16*len(v)+8*(n+1))
  options=dict(rtol=1e-10,atol=0,limit=2000,capture=0,jacobi=0)
  for operation,op,mode,expected in [('ic0_setup','ic0_factor',0,rp+ci+v),('ic0_apply','ic0_apply',0,apply(q,rhs,exact=False)),('cg','cg',0,x),('jacobi','cg',1,x),('ic0_solve','cg',2,x),('ic0_total','cg',3,x),('cholesky_solve','chol_solve',0,x),('cholesky_total','chol_total',0,x)]:
   rowstats=stats.copy()
   if operation.startswith('cholesky'):
    _,c,steps=analyze(q);rowstats.update(factor_nnz=len(c),fill_count=sum(k>=0 for k in steps),logical_factor_bytes=16*len(c)+8*(n+1))
   elif operation in ('cg','jacobi'):rowstats.update(factor_nnz=0,logical_factor_bytes=0)
   yield operation,ordering,dict(name='IC0 '+operation,op=op,a=q,b=rhs,expected=expected,invalid=False,reason='converged',options={**options,'jacobi':mode}),rowstats

def render(data,destination,*,live_override=None, published=False):
 data = (PublicSanitizer.for_published() if published else PublicSanitizer()).report(data)
 live=sparse.live_bundle(data,{'ic0_source':(HERE/'ic0-live.mjs').read_text()}) if live_override is None else live_override
 replacements={'DATA':json_for_html(data),'LIVE':json_for_html(live),'SCRIPT':(HERE/'ic0-report.mjs').read_text()}
 html=re.sub(r'@@(DATA|LIVE|SCRIPT)@@',lambda m:replacements[m[1]],(HERE/'ic0-report.html').read_text())
 Path(destination).write_text(apply_report_design(html))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--no-build',action='store_true');p.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));p.add_argument('--sizes',nargs='+',type=int,default=[8]);p.add_argument('--samples',type=int,default=3);p.add_argument('--output',type=Path,default=HERE/'reports/ic0');p.add_argument('--render-only',type=Path);args=p.parse_args()
 if args.render_only:render(json.loads(args.render_only.read_text()),args.render_only.with_name('index.html'));return 0
 if not 1<=args.samples<=100 or any(n<2 or n>16 for n in args.sizes):p.error('sizes must be 2–16, samples 1–100')
 implementations=[];start=time.perf_counter()
 for name in args.languages:
  impl=dict(id=name,name=run.LABELS[name])
  try:impl=run.build_one(name,no_build=args.no_build);sparse.verify(impl,fixtures())
  except Exception as e:impl.update(status='failed',error=str(e))
  implementations.append(impl);print(name,impl['status'],flush=True)
 available=[i for i in implementations if i['status']=='passed']
 # Reuse the serial sampler, restricted to one explicitly described problem family.
 results=benchmark(available,args.sizes,args.samples)
 data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='sparse-ic0-v1',seed=2026,total_seconds=time.perf_counter()-start,machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],methodology=dict(numeric_type='IEEE 754 binary64',matrix_layout='SPD five-point diffusion with conductivity 1–100. Same system under natural, RCM and AMD permutations. Size is grid width; unknowns = size squared.',timing='Serial in-process timings. IC0 setup includes factor construction, lower export and checksum. IC0 apply reuses a factor. Prepared CG and Cholesky solves reuse factors across repeated calls with the same RHS. IC0 total constructs a fresh factor and solves; Cholesky total analyzes, factors and solves. Totals are directly measured on already-permuted systems; ordering and permutation costs are excluded from every row. Plain/Jacobi CG include their own validation and diagonal preparation. Parsing, process startup, compilation and initial CSR creation excluded. Explicit disposal included; managed GC follows runtime policy.',warmup='Three warmups per process; calibration targets 20 ms, capped at 10,000 calls.',sampling='Three serial samples by default; seeded randomized language order. Median and MAD. Each result is validated before timing, including true residuals for CG.',storage='Logical factor bytes: 8-byte values and indices plus 8-byte row offsets. Excludes source matrix, transient symbolic plan, solver workspace, permutation and runtime overhead. Plain/Jacobi factor storage is reported as zero; Jacobi still uses a diagonal vector.',profiles='No sampled profiles in this suite. Live frames show accepted algorithm iterates, not stack frames or wall-clock time.',limitations='Educational unshifted IC(0), with no pivoting or fill beyond the original symmetric stored pattern. Nonpositive or nonfinite pivots fail explicitly, even for some SPD inputs. CG requires SPD input; positive IC0 pivots alone do not certify it. Ordering affects both approximation quality and convergence. One small workload and three samples do not predict general performance. Live timings are separate single observations.'))
 data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');render(data,args.output/'index.html')
 return int(len(available)!=len(implementations) or any(r['status']!='passed' for r in results))

def benchmark(implementations,sizes,samples,case_factory=None):
 import random,statistics
 from sparse_reference import protocol
 rows=[];rng=random.Random(2026)
 for size in sizes:
  for operation,ordering,case,stats in (case_factory or cases)(size):
   current=[]
   for impl in implementations:
    row=dict(implementation=impl['id'],operation=operation,ordering=ordering,problem='grid',size=size,unknowns=size*size,nnz=len(case['a']['values']),**stats,status='passed',samples=[])
    try:
     args,stdin=protocol(case);r=run.execute(impl['runner']+args,stdin=stdin,env=impl['env'],timeout=120)
     if r.returncode:raise ValueError(r.stderr)
     actual=json.loads(r.stdout);sparse.check_result(actual,case)
     if case['op']=='cg':row['solver_iterations']=int(actual['values'][1])
     calibration=sparse.measure(impl,case,1);row['iterations']=max(1,min(10000,math.ceil(20e6/calibration['elapsed_ns'])))
    except Exception as e:row.update(status='failed',error=str(e) or 'incorrect result')
    current.append((impl,row))
   for _ in range(samples):
    rng.shuffle(current)
    for impl,row in current:
     if row['status']=='passed':
      try:row['samples'].append(sparse.measure(impl,case,row['iterations']))
      except Exception as e:row.update(status='failed',error=str(e))
   for _,row in current:
    if row['status']=='passed':
     values=[s['ns_per_op'] for s in row['samples']];median=statistics.median(values);row.update(median_ns=median,min_ns=min(values),max_ns=max(values),mad_ns=statistics.median(abs(v-median) for v in values))
    rows.append(row);print(row['implementation'],size,ordering,operation,row['status'],row.get('error',''),flush=True)
 return rows
if __name__=='__main__':raise SystemExit(main())
