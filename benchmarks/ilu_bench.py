#!/usr/bin/env python3
"""ILU(0) setup, application and GMRES reuse comparison across eleven ports."""
import argparse,datetime,json,math,os,platform,random,statistics,time
from pathlib import Path
import run,sparse
from publication import PublicSanitizer
from gmres_bench import workspace_bytes
from gmres_reference import transport,check_result as check_gmres
from ilu_reference import fixtures,reference_apply
from sparse_reference import multiply,protocol
HERE=Path(__file__).resolve().parent

def cases(size):
    a=transport(size,2,.2,4,30);n=a['rows'];b=multiply(a,[1.]*n)
    options=dict(restart=20,rtol=1e-10,atol=0,limit=3000,jacobi=0,capture=0)
    base=dict(name='ILU benchmark',a=a,b=b,invalid=False,reason='converged')
    yield 'ilu_setup',{**base,'op':'ilu_setup','expected':[n,len(a['values'])],'options':options}
    yield 'ilu_apply',{**base,'op':'ilu_apply','expected':reference_apply(a,b,exact=False),'options':options}
    for mode,name in enumerate(('gmres_none','gmres_jacobi','gmres_ilu_reused','gmres_ilu_total')):
        yield name,{**base,'op':'gmres','expected':[1.]*n,'options':{**options,'jacobi':mode}}

def benchmark(implementations,sizes,samples):
    results=[];rng=random.Random(2026)
    for size in sizes:
        for operation,case in cases(size):
            a=case['a'];n=a['rows'];nnz=len(a['values']);rows=[]
            for impl in implementations:
                row=dict(implementation=impl['id'],operation=operation,size=size,unknowns=n,nnz=nnz,contrast=2,restart=20,logical_dense_bytes=8*n*n,logical_csr_bytes=16*nnz+8*(n+1),logical_preconditioner_bytes=16*nnz+16*n+8 if 'ilu' in operation else 0,logical_workspace_bytes=workspace_bytes(n,20,3000) if case['op']=='gmres' else 8*n if case['op']=='ilu_apply' else 0,status='passed',samples=[])
                try:
                    args,stdin=protocol(case);r=run.execute(impl['runner']+args,stdin=stdin,env=impl['env'],timeout=120)
                    if r.returncode:raise ValueError(r.stderr)
                    actual=json.loads(r.stdout);sparse.check_result(actual,case)
                    if case['op']=='gmres':
                        checked=check_gmres(actual,case);row.update(solver_iterations=checked['iterations'],restart_count=len(checked['restarts']))
                    calibration=sparse.measure(impl,case,1);row['iterations']=max(1,min(10000,math.ceil(20e6/calibration['elapsed_ns'])))
                except Exception as error:row.update(status='failed',error=str(error) or 'incorrect result')
                rows.append((impl,row))
            for _ in range(samples):
                rng.shuffle(rows)
                for impl,row in rows:
                    if row['status']!='passed':continue
                    try:row['samples'].append(sparse.measure(impl,case,row['iterations']))
                    except Exception as error:row.update(status='failed',error=str(error))
            for _,row in rows:
                if row['status']=='passed':
                    values=[s['ns_per_op'] for s in row['samples']];median=statistics.median(values);row.update(median_ns=median,min_ns=min(values),max_ns=max(values),mad_ns=statistics.median(abs(v-median) for v in values))
                results.append(row);print(row['implementation'],operation,size,row['status'],row.get('error',''),flush=True)
    return results

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--no-build',action='store_true');p.add_argument('--verify-only',action='store_true');p.add_argument('--require-all',action='store_true');p.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));p.add_argument('--sizes',nargs='+',type=int,default=[8,16]);p.add_argument('--samples',type=int,default=3);p.add_argument('--output',type=Path,default=HERE/'reports/ilu');p.add_argument('--render-only',type=Path);args=p.parse_args()
    if args.render_only:sparse.render(json.loads(args.render_only.read_text()),args.render_only.with_name('index.html'));return 0
    if not 1<=args.samples<=100 or any(n<2 or n>32 for n in args.sizes):p.error('samples must be 1–100 and sizes 2–32')
    started=time.perf_counter();implementations=[]
    for name in args.languages:
        i=dict(id=name,name=run.LABELS[name])
        try:i=run.build_one(name,no_build=args.no_build);sparse.verify(i,fixtures())
        except run.ToolchainUnavailable as error:i.update(status='unavailable',error=str(error))
        except Exception as error:i.update(status='failed',error=str(error))
        implementations.append(i);print(name,i['status'],[c for c in i.get('checks',[]) if not c['passed']],flush=True)
    available=[i for i in implementations if i['status']=='passed'];results=[] if args.verify_only else benchmark(available,args.sizes,args.samples)
    failed=not available or any(i['status']=='failed' or(args.require_all and i['status']=='unavailable') for i in implementations) or any(r['status']!='passed' for r in results)
    if not args.verify_only:
        data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='ilu-reuse-v1',seed=2026,total_seconds=time.perf_counter()-started,machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],methodology=dict(
            numeric_type='IEEE 754 binary64',matrix_layout='Canonical CSR, 1–100 horizontal conductivity, diffusivity 0.2, speed 4, angle 30 degrees clockwise from right; zero Dirichlet boundaries. Known x=ones defines b=A*x.',
            timing='In-process monotonic clock; excludes parsing, compilation, process startup and CSR construction. ILU setup includes factor construction, a two-value shape checksum, and explicit destruction where applicable; managed collection follows runtime policy. ILU apply includes triangular solves, validation, result allocation and checksum; factor setup is outside timing. Reused GMRES also excludes factor setup. One-shot total includes factor setup and the solve. No-preconditioner and Jacobi runs include all their solver preparation.',
            warmup='Three excluded warmup calls in each process, including Julia checksum-loop compilation. Separate calibration targets 20 ms, capped at 10,000 repetitions.',sampling='Three samples by default; serial seeded randomized port order, median and MAD. Timed repetitions use the same RHS, with no captured frames. Reuse across different RHS is tested separately.',profiles='No profiles collected for this suite.',
            limitations='GMRES restart 20, rtol 1e-10, atol 0, limit 3000; true residual recomputed every accepted update. ILU(0) has no fill, pivoting, permutation, shifts or fallback. It can break down on nonsingular matrices. Factor memory is a logical model with float64 values and 64-bit indices: 16*nnz+16*n+8 bytes; excludes objects, allocator overhead and transient copies. Solver workspace uses the GMRES logical model; no RSS measurements. Setup and solve samples are independent, so their medians need not sum to the measured one-shot median. Amortized charts are estimates from separate measurements, not measured multi-RHS timings. Small systems and JIT/GC effects limit generalization.'))
        data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');sparse.render(data,args.output/'index.html')
    return int(failed)
if __name__=='__main__':raise SystemExit(main())
