#!/usr/bin/env python3
"""Compare restarted GMRES on nonsymmetric upwind transport systems."""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import time
import run
import sparse
from sparse_reference import multiply,protocol
from gmres_reference import fixtures,check_result,transport
from publication import PublicSanitizer

HERE=Path(__file__).resolve().parent

def workspace_bytes(n,restart,limit):
    m=min(n,restart,limit)
    return 8*((m+1)*n+(m+1)*m+8*n+4*(m+1))


def benchmark(implementations,sizes,restart_lengths,samples):
    results=[];rng=random.Random(2026)
    for size in sizes:
        a=transport(size,1,.2,4,30);n=a['rows'];b=multiply(a,[1.]*n)
        for restart in restart_lengths:
            for jacobi in [False,True]:
                operation='gmres_restart_%d%s'%(restart,'_jacobi' if jacobi else '')
                case=dict(name='GMRES benchmark',a=a,b=b,op='gmres',expected=[1.]*n,reason='converged',options=dict(restart=restart,rtol=1e-10,atol=0,limit=1500,jacobi=int(jacobi),capture=0))
                rows=[]
                for impl in implementations:
                    row=dict(implementation=impl['id'],operation=operation,size=size,unknowns=n,nnz=len(a['values']),restart=restart,logical_workspace_bytes=workspace_bytes(n,restart,1500),logical_dense_bytes=8*n*n,logical_csr_bytes=16*len(a['values'])+8*(n+1),contrast=1,status='passed',samples=[])
                    try:
                        args,stdin=protocol(case);r=run.execute(impl['runner']+args,stdin=stdin,env=impl['env'],timeout=120)
                        if r.returncode:raise ValueError(r.stderr)
                        checked=check_result(json.loads(r.stdout),case);row.update(solver_iterations=checked['iterations'],restart_count=len(checked['restarts']))
                        calibration=sparse.measure(impl,case,1);row['iterations']=max(1,min(10000,math.ceil(20e6/calibration['elapsed_ns'])))
                    except Exception as error:row.update(status='failed',error=str(error) or 'incorrect GMRES result')
                    rows.append((impl,row))
                for _ in range(samples):
                    rng.shuffle(rows)
                    for impl,row in rows:
                        if row['status']!='passed':continue
                        try:row['samples'].append(sparse.measure(impl,case,row['iterations']))
                        except Exception as error:row.update(status='failed',error=str(error))
                for _,row in rows:
                    if row['status']=='passed':
                        values=[s['ns_per_op'] for s in row['samples']];median=statistics.median(values)
                        row.update(median_ns=median,min_ns=min(values),max_ns=max(values),mad_ns=statistics.median(abs(v-median) for v in values))
                    results.append(row);print(row['implementation'],operation,size,row['status'],row.get('error',''),flush=True)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-build',action='store_true');parser.add_argument('--verify-only',action='store_true');parser.add_argument('--require-all',action='store_true')
    parser.add_argument('--languages',nargs='+',choices=run.NAMES,default=list(run.NAMES));parser.add_argument('--sizes',nargs='+',type=int,default=[8,16]);parser.add_argument('--restarts',nargs='+',type=int,default=[5,20]);parser.add_argument('--samples',type=int,default=3)
    parser.add_argument('--output',type=Path,default=HERE/'reports/gmres');parser.add_argument('--render-only',type=Path)
    args=parser.parse_args()
    if args.render_only:sparse.render(json.loads(args.render_only.read_text()),args.render_only.with_name('index.html'));return 0
    if not 1<=args.samples<=100 or any(s<2 or s>32 for s in args.sizes) or any(r<1 or r>100 for r in args.restarts):parser.error('samples 1–100, grid widths 2–32, restart lengths 1–100')
    started=time.perf_counter();implementations=[]
    for name in args.languages:
        i=dict(id=name,name=run.LABELS[name])
        try:i=run.build_one(name,no_build=args.no_build);sparse.verify(i,fixtures())
        except run.ToolchainUnavailable as error:i.update(status='unavailable',error=str(error))
        except Exception as error:i.update(status='failed',error=str(error))
        implementations.append(i);print(name,i['status'],[c for c in i.get('checks',[]) if not c['passed']],flush=True)
    available=[i for i in implementations if i['status']=='passed'];results=[] if args.verify_only else benchmark(available,args.sizes,args.restarts,args.samples)
    failed=not available or any(i['status']=='failed' or (args.require_all and i['status']=='unavailable') for i in implementations) or any(r['status']!='passed' for r in results)
    if not args.verify_only:
        data=dict(schema_version=1,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=run.require(['git','rev-parse','HEAD']).strip(),dirty=bool(run.require(['git','status','--porcelain']).strip()),source_sha256=run.fingerprint(),suite='gmres-transport-v1',seed=2026,total_seconds=time.perf_counter()-started,
                  machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),implementations=implementations,results=results,profiles=[],
                  methodology=dict(numeric_type='IEEE 754 binary64',matrix_layout='Canonical CSR, five-point diffusion plus first-order upwind transport. Conductivity 1–10, diffusivity 0.2, speed 4, angle 30 degrees clockwise from right, zero Dirichlet boundaries. Known x=ones defines b=A*x.',
                      timing='In-process monotonic timer. Include solver validation, allocations, copies, true-residual checks and checksum; exclude parsing, CSR construction, compilation, serialization and process startup.',
                      warmup='Three warmup calls per fresh process; Julia warms the identical checksum loop. One excluded calibration batch targets 20 ms, capped at 10,000 calls.',sampling='Serial samples, seeded randomized port order each round; median and MAD. No frame capture during timed runs.',profiles='No profiles collected for this suite.',
                      limitations='GMRES uses two-pass modified Gram–Schmidt, Givens rotations, right Jacobi when selected, x0=0, rtol=1e-10, atol=0, total limit=1500. True residual is recomputed every accepted step. Workspace is a logical float64 model: 8*((m+1)*n+(m+1)*m+8*n+4*(m+1)) bytes, m=min(restart,n,limit). Excludes matrix storage, results, histories, snapshots, runtime objects, allocator overhead and transient copies; not measured RSS. Restart can cause stagnation; nonsingular input does not guarantee convergence. Small workloads, JIT/GC and machine scheduling limit generalization; these are implementation comparisons, not language rankings.'))
        data=PublicSanitizer().report(data);args.output.mkdir(parents=True,exist_ok=True);(args.output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');sparse.render(data,args.output/'index.html')
    return int(failed)

if __name__=='__main__':raise SystemExit(main())
