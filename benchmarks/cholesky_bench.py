#!/usr/bin/env python3
"""Measure sparse Cholesky analysis, factorization, repeated solves and fill."""
import argparse,datetime,json,math,os,platform,random,re,statistics,time
from pathlib import Path
import run,sparse
from publication import PublicSanitizer
from report import json_for_html
from report_design import apply_report_design
from ordering_reference import rcm,permute
from cholesky_reference import analyze,factor,fixtures
from sparse_reference import diffusion
from sparse_reference import multiply,protocol
HERE=Path(__file__).resolve().parent

def matrix(size,problem):
    a=diffusion(size,2)
    if problem=='scrambled':
        p=list(range(a['rows']));random.Random(2026).shuffle(p);a=permute(a,p)
    return a

def cases(size,problem):
    a=matrix(size,problem);n=a['rows'];p=rcm(a);q=permute(a,p)
    x=[1+(i%7)/10 for i in range(n)];b=multiply(a,x)
    for order,mat,rhs,truth in [('natural',a,b,x),('rcm',q,[b[i] for i in p],[x[i] for i in p])]:
        rp,ci,steps=analyze(mat);stats=dict(factor_nnz=len(ci),fill_count=sum(k>=0 for k in steps),logical_factor_bytes=16*len(ci)+8*(n+1))
        for operation,op,expected in [('symbolic','chol_symbolic',rp+ci+steps),('factor','chol_factor',sum(factor(mat),[])),('solve','chol_solve',truth)]:
            yield operation,order,dict(name=problem+' '+op,op=op,a=mat,b=rhs,expected=expected,invalid=False),stats
        op='chol_total' if order=='natural' else 'chol_rcm_total'
        yield 'total',order,dict(name=problem+' '+op,op=op,a=a,b=b,expected=x,invalid=False),stats

def benchmark(implementations,sizes,samples):
    rows=[];rng=random.Random(2026)
    for size in sizes:
      for problem in ('grid','scrambled'):
       for operation,ordering,case,stats in cases(size,problem):
        current=[];a=case['a'];n=a['rows'];nnz=len(a['values'])
        for impl in implementations:
            row=dict(implementation=impl['id'],operation=operation,ordering=ordering,problem=problem,size=size,unknowns=n,nnz=nnz,**stats,status='passed',samples=[])
            timed_case=case
            try:
                args,stdin=protocol(case);r=run.execute(impl['runner']+args,stdin=stdin,env=impl['env'],timeout=120)
                if r.returncode:raise ValueError(r.stderr)
                actual=json.loads(r.stdout);checked=sparse.check_result(actual,case)
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
    data=PublicSanitizer().report(data);live=sparse.live_bundle(data,{'cholesky_source':(HERE/'cholesky-live.mjs').read_text()}) if live_override is None else live_override
    template=(HERE/'cholesky-report.html').read_text()
    replacements={'DATA':json_for_html(data),'LIVE':json_for_html(live),'SCRIPT':(HERE/'cholesky-report.mjs').read_text()}
    html=re.sub(r'@@(DATA|LIVE|SCRIPT)@@',lambda m:replacements[m[1]],template)
    Path(destination).write_text(apply_report_design(html))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--no-build',action='store_true');p.add_argument('--verify-only',action='store_true');p.add_argument('--require-all',action='store_true');p.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));p.add_argument('--sizes',nargs='+',type=int,default=[12]);p.add_argument('--samples',type=int,default=3);p.add_argument('--output',type=Path,default=HERE/'reports/cholesky');p.add_argument('--render-only',type=Path);args=p.parse_args()
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
        data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='sparse-cholesky-v1',seed=2026,total_seconds=time.perf_counter()-start,machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],methodology=dict(numeric_type='IEEE 754 binary64',matrix_layout='Symmetric positive-definite five-point diffusion, conductivity 1–100. Grid order and seeded scrambled order describe the same system.',timing='In-process monotonic clock. Symbolic includes analysis and packed pattern checksum. Factor includes numerical factorization with an existing plan, lower-factor export and checksum. Solve reuses an existing factor, including solution allocation and checksum. Total directly measures analysis, factorization, solve and, for RCM, ordering, matrix/RHS permutation and restoring original coordinates. Explicit cleanup is included; managed GC follows runtime policy. Parsing, process startup, compilation and initial CSR creation excluded.',warmup='Three warmups per process. Calibration targets 20 ms, capped at 10,000 repetitions.',sampling='Three serial samples by default with seeded randomized port order. Median and MAD. Independent dense Cholesky, relative reconstruction error and original-coordinate solution residual checked before timing.',storage='Logical lower-factor storage: 8-byte values plus 8-byte indices and row offsets. Excludes symbolic plan, source-pattern snapshot, temporary graph, permutation and runtime overhead. Actual native index widths and object layouts vary.',profiles='No profiles collected for this suite. The live animation shows actual symbolic elimination steps, not sampled execution time.',limitations='Educational row-oriented Cholesky with sparse row intersections; no supernodes, pivoting, shifts or BLAS. Exact numeric symmetry required. Nonfinite arithmetic and nonpositive computed pivots rejected. Positive pivots do not certify exact-real definiteness under arbitrary roundoff. RCM is a bandwidth heuristic and need not minimize fill or runtime; AMD is not included. Small workloads and three samples limit generalization.'))
        data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');render(data,args.output/'index.html')
    return int(not available or any(i['status']=='failed' or (args.require_all and i['status']=='unavailable') for i in implementations) or any(r['status']!='passed' for r in results))
if __name__=='__main__':raise SystemExit(main())
