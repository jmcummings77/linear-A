#!/usr/bin/env python3
"""Measure natural and RCM ordering, including complete ILU/GMRES solves."""
import argparse,base64,datetime,json,math,os,platform,random,re,statistics,time
from pathlib import Path
import run,sparse
from publication import PublicSanitizer
from report import json_for_html
from report_design import apply_report_design
from ordering_reference import rcm,permute,bandwidth,fixtures,check_solve
from gmres_reference import transport,check_result as check_gmres
from sparse_reference import multiply,protocol
HERE=Path(__file__).resolve().parent

def matrix(size,problem):
    a=transport(size,2,.2,4,30)
    if problem=='scrambled':
        p=list(range(a['rows']));random.Random(2026).shuffle(p);a=permute(a,p)
    return a

def cases(size,problem):
    a=matrix(size,problem);n=a['rows'];p=rcm(a);q=permute(a,p)
    # Nonconstant solution catches accidentally omitting either permutation.
    x=[1+(i%7)/10 for i in range(n)];b=multiply(a,x)
    o=dict(restart=20,rtol=1e-10,atol=0,limit=3000,jacobi=2,capture=0)
    def case(op,a,b,expected,options=o):return dict(name=problem+' '+op,op=op,a=a,b=b,expected=expected,invalid=False,reason='converged',options=options)
    yield 'rcm','rcm',case('rcm',a,[0]*n,p,{**o,'jacobi':0}),bandwidth(q)
    yield 'permute','rcm',case('permute',a,p,q['offsets']+q['indices']+q['values'],{**o,'jacobi':0}),bandwidth(q)
    for order,mat,rhs,truth in [('natural',a,b,x),('rcm',q,[b[i] for i in p],[x[i] for i in p])]:
        yield 'setup',order,case('ilu_setup',mat,[0]*n,[n,len(mat['values'])],{**o,'jacobi':0}),bandwidth(mat)
        yield 'solve',order,case('gmres',mat,rhs,truth),bandwidth(mat)
    for order,op in [('natural','ilu_solve'),('rcm','rcm_solve')]:
        yield 'total',order,case(op,a,b,x,{**o,'jacobi':0}),bandwidth(a if order=='natural' else q)

def benchmark(implementations,sizes,samples):
    rows=[];rng=random.Random(2026)
    for size in sizes:
      for problem in ('grid','scrambled'):
       for operation,ordering,case,width in cases(size,problem):
        current=[];a=case['a'];n=a['rows'];nnz=len(a['values'])
        for impl in implementations:
            row=dict(implementation=impl['id'],operation=operation,ordering=ordering,problem=problem,size=size,unknowns=n,nnz=nnz,bandwidth=width,status='passed',samples=[])
            timed_case=case
            try:
                args,stdin=protocol(case);r=run.execute(impl['runner']+args,stdin=stdin,env=impl['env'],timeout=120)
                if r.returncode:raise ValueError(r.stderr)
                actual=json.loads(r.stdout);checked=sparse.check_result(actual,case)
                if case['op']=='gmres':row['solver_iterations']=check_gmres(actual,case)['iterations']
                elif case['op'].endswith('_solve'):
                    row['solver_iterations']=check_solve(actual,case)['iterations']
                    timed_case={**case,'expected':[row['solver_iterations']]+case['expected']}
                calibration=sparse.measure(impl,timed_case,1);row['iterations']=max(1,min(10000,math.ceil(20e6/calibration['elapsed_ns'])))
            except Exception as e:row.update(status='failed',error=str(e) or 'incorrect result')
            current.append((impl,row,timed_case))
        for _ in range(samples):
            rng.shuffle(current)
            for impl,row,timed_case in current:
                if row['status']=='passed':
                    try:row['samples'].append(sparse.measure(impl,timed_case,row['iterations']))
                    except Exception as e:row.update(status='failed',error=str(e))
        for _,row,_ in current:
            if row['status']=='passed':
                values=[s['ns_per_op'] for s in row['samples']];median=statistics.median(values);row.update(median_ns=median,min_ns=min(values),max_ns=max(values),mad_ns=statistics.median(abs(v-median) for v in values))
            rows.append(row);print(row['implementation'],problem,size,operation,ordering,row['status'],row.get('error',''),flush=True)
    return rows

def render(data,destination,*,live_override=None):
    data=PublicSanitizer().report(data);live=sparse.live_bundle(data,{'ordering_source':(HERE/'ordering-live.mjs').read_text()}) if live_override is None else live_override
    template=(HERE/'ordering-report.html').read_text()
    replacements={'DATA':json_for_html(data),'LIVE':json_for_html(live),'SCRIPT':(HERE/'ordering-report.mjs').read_text()}
    html=re.sub(r'@@(DATA|LIVE|SCRIPT)@@',lambda m:replacements[m[1]],template)
    Path(destination).write_text(apply_report_design(html))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--no-build',action='store_true');p.add_argument('--verify-only',action='store_true');p.add_argument('--require-all',action='store_true');p.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));p.add_argument('--sizes',nargs='+',type=int,default=[12]);p.add_argument('--samples',type=int,default=3);p.add_argument('--output',type=Path,default=HERE/'reports/ordering');p.add_argument('--render-only',type=Path);args=p.parse_args()
    if args.render_only:render(json.loads(args.render_only.read_text()),args.render_only.with_name('index.html'));return 0
    if not 1<=args.samples<=100 or any(n<2 or n>32 for n in args.sizes):p.error('sizes must be 2–32 and samples 1–100')
    start=time.perf_counter();implementations=[]
    for name in args.languages:
        impl=dict(id=name,name=run.LABELS[name])
        try:impl=run.build_one(name,no_build=args.no_build);sparse.verify(impl,fixtures())
        except run.ToolchainUnavailable as e:impl.update(status='unavailable',error=str(e))
        except Exception as e:impl.update(status='failed',error=str(e))
        implementations.append(impl);print(name,impl['status'],[c for c in impl.get('checks',[]) if not c['passed']],flush=True)
    available=[i for i in implementations if i['status']=='passed'];results=[] if args.verify_only else benchmark(available,args.sizes,args.samples)
    if not args.verify_only:
        data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='rcm-ordering-v1',seed=2026,total_seconds=time.perf_counter()-start,machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],methodology=dict(numeric_type='IEEE 754 binary64',matrix_layout='Five-point transport: conductivity 1–100, diffusivity 0.2, speed 4, angle 30°. Grid order and seeded scrambled order represent the same operator. Stored zeros count as graph edges; diagonal ignored; nonsymmetric pattern symmetrized by union.',timing='In-process monotonic clock. RCM: ordering plus checksum. Permute: CSR permutation plus packed-array checksum. Setup: ILU construction and shape checksum. Solve: GMRES with factor prepared outside timing. Total: measured ordering (RCM only), CSR and RHS permutation (RCM only), ILU construction, GMRES, restoring original variable order (RCM only), checksum and explicit cleanup. Managed GC follows runtime policy. Parsing, process startup, compilation and input CSR creation excluded.',warmup='Three warmup calls per process, 20 ms calibration target capped at 10,000 repetitions.',sampling='Three serial samples by default, seeded randomized port order; median and MAD. Each repetition starts from zero; restart 20, rtol 1e-10, limit 3000. Nonconstant known solution and independent true residual in original coordinates checked before timing.',profiles='No profiles collected for this suite.',limitations='RCM uses minimum-degree component seeds, degree/index neighbor ties, and reverses the whole traversal; no pseudo-peripheral search. It is a heuristic: neither smaller bandwidth nor faster convergence is guaranteed. ILU(0) has no pivoting or fill; failures are retained. Preconditioned solve times exclude preprocessing; total times include it and are independently measured. Logical CSR/ILU storage is unchanged by permutation; permutation arrays and graph work add temporary storage. Small workloads, JIT/GC effects and three samples limit generalization.'))
        data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');render(data,args.output/'index.html')
    return int(not available or any(i['status']=='failed' or (args.require_all and i['status']=='unavailable') for i in implementations) or any(r['status']!='passed' for r in results))
if __name__=='__main__':raise SystemExit(main())
