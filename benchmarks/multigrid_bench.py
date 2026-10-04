#!/usr/bin/env python3
"""Structured multigrid, IC(0), and CG scaling with reusable and fresh setup."""
import argparse,datetime,json,math,os,platform,re,time
from pathlib import Path
import run,sparse
from multigrid_reference import poisson,truth,fixtures,reference
from sparse_reference import multiply
from publication import PublicSanitizer
from report import json_for_html
from report_design import apply_report_design
from ic0_bench import benchmark as sample_cases
HERE=Path(__file__).resolve().parent

def stencil_ic0(width):
 """Independent no-triangle stencil recurrence, not a production factor call."""
 n=width*width;rp=[0];ci=[];v=[];diag=[]
 for i in range(n):
  terms=[]
  for j in ([i-width] if i>=width else [])+([i-1] if i%width else []):
   value=-1/diag[j];ci.append(j);v.append(value);terms.append(value*value)
  diag.append(math.sqrt(4-sum(terms)));ci.append(i);v.append(diag[-1]);rp.append(len(v))
 return rp,ci,v

def cases(size):
 a=poisson(size);n=size*size;x=truth(size);b=multiply(a,x);rp,ci,v=stencil_ic0(size)
 opts=dict(rtol=1e-10,atol=0,limit=2000,capture=0,jacobi=0)
 operations=[('mg_setup','mg_setup',0,[n,(size+1).bit_length()-1]),('ic0_setup','ic0_factor',0,rp+ci+v),('cg','cg',0,x),('jacobi','cg',1,x),('ic0_solve','cg',2,x),('ic0_total','cg',3,x),('mg_solve','cg',4,x),('mg_total','cg',5,x)]
 if size==7:operations.append(('mg_apply','mg_apply',0,list(map(float,reference(size,b)))))
 for operation,op,mode,expected in operations:
  storage=16*len(v)+8*(n+1) if operation.startswith('ic0') else 8 if operation.startswith('mg') else 0
  # Live arrays differ by language: report a portable mathematical workspace model.
  workspace=8*(5*n if op=='cg' else 0)
  if operation.startswith('mg'):
   w=size;cycle_values=0
   while w:cycle_values+=w*w+2*(w//2)**2;w//=2
   workspace+=8*(n+cycle_values)
  yield operation,'natural',dict(name='Multigrid '+operation,a=a,b=b,op=op,expected=expected,invalid=False,reason='converged',options={**opts,'jacobi':mode}),dict(factor_nnz=len(v) if operation.startswith('ic0') else 0,logical_factor_bytes=storage,logical_workspace_bytes=workspace)

def benchmark(implementations,sizes,samples):
 return sample_cases(implementations,sizes,samples,case_factory=cases)

def render(data,destination,*,live_override=None, published=False):
 data = (PublicSanitizer.for_published() if published else PublicSanitizer()).report(data)
 live=sparse.live_bundle(data,{'multigrid_source':(HERE/'multigrid-live.mjs').read_text()}) if live_override is None else live_override
 replacements={'DATA':json_for_html(data),'LIVE':json_for_html(live),'SCRIPT':(HERE/'multigrid-report.mjs').read_text()}
 html=re.sub(r'@@(DATA|LIVE|SCRIPT)@@',lambda m:replacements[m[1]],(HERE/'multigrid-report.html').read_text())
 Path(destination).write_text(apply_report_design(html))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--no-build',action='store_true');p.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));p.add_argument('--sizes',nargs='+',type=int,default=[7,15,31]);p.add_argument('--samples',type=int,default=3);p.add_argument('--output',type=Path,default=HERE/'reports/multigrid');p.add_argument('--render-only',type=Path);args=p.parse_args()
 if args.render_only:render(json.loads(args.render_only.read_text()),args.render_only.with_name('index.html'));return 0
 if not 1<=args.samples<=100 or any(n not in (1,3,7,15,31) for n in args.sizes):p.error('sizes must be nested widths through 31; samples 1–100')
 implementations=[];start=time.perf_counter()
 for name in args.languages:
  impl=dict(id=name,name=run.LABELS[name])
  try:impl=run.build_one(name,no_build=args.no_build);sparse.verify(impl,fixtures())
  except Exception as e:impl.update(status='failed',error=str(e))
  implementations.append(impl);print(name,impl['status'],flush=True)
 available=[i for i in implementations if i['status']=='passed'];results=benchmark(available,args.sizes,args.samples)
 data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='geometric-multigrid-v1',seed=2026,total_seconds=time.perf_counter()-start,machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],methodology=dict(numeric_type='IEEE 754 binary64',matrix_layout='Unit constant-coefficient five-point Dirichlet Laplacian on nested square interior grids. Width 2^k−1; zero boundary; unscaled diagonal 4 and neighbor entries −1. Deterministic broad-spectrum solution, b=A x. Natural ordering only.',timing='Serial in-process measurements. Prepared solves reuse the preconditioner. Total rows construct a fresh preconditioner and solve, measured directly. IC0 setup includes lower export and checksum; MG setup validates a width and returns metadata for its implicit hierarchy. Matrix assembly, process startup and parsing excluded. Apply allocates its workspaces on every call. No BLAS or external solver.',warmup='Three warmups per process. Calibration targets 20 ms, capped at 10,000 calls.',sampling='Three serial samples by default, with seeded randomized language order; median and MAD. Correctness checked before timing.',storage='Logical retained bytes: IC0 values and indices at 8 bytes plus row offsets; MG one 8-byte width with implicit levels and transfers. Workspace is a mathematical array model: CG five n-vectors; MG output plus each active level residual and coarse RHS/correction. Excludes source CSR, runtime objects, allocator overhead, language-specific temporaries and capture frames. Setup rows show the apply workspace model for comparison, not setup peak allocation.',profiles='No sampled profiles. Live V-cycle snapshots come from actual C/WASM function callbacks. Playback is algorithm progress, not measured CPU time.',limitations='Structured isotropic constant-coefficient diffusion only. Two weighted-Jacobi sweeps before and after correction, weight 2/3. Bilinear P and restriction Pᵀ with rediscretized coarse five-point operators, not Galerkin operators. The unscaled stencil requires transpose restriction rather than normalized full weighting. Coarsest solve is exact. Fixed symmetric V-cycle is SPD and usable with CG; it is not an algebraic multigrid solver. Jacobi is scalar scaling for this operator, so it does not improve CG iteration counts. Small-grid timings and three samples do not establish general performance.'))
 data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');render(data,args.output/'index.html')
 return int(len(available)!=len(implementations) or any(r['status']!='passed' for r in results))
if __name__=='__main__':raise SystemExit(main())
