#!/usr/bin/env python3
"""Reproducible multiplication-locality study; Python standard library only."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = HERE / '.build'
sys.path.insert(0, str(ROOT / 'benchmarks'))
from publication import PublicSanitizer
from report_design import apply_report_design
from profile_publication import sanitize_macos_sample
from profiles import parse_macos_sample

VARIANTS = {
    'ijk': ('Scalar i-j-k', 'One dot product per output; strided reads from B.'),
    'ikj': ('Scalar i-k-j', 'Stream across a row of B and C for each A entry.'),
    'blocked': ('32 × 32 tiles', 'Bound the active submatrices; scalar inner loops.'),
    'neon': ('Tiles + ARM64 NEON', 'Two adjacent output columns per explicit SIMD operation.'),
    'blas': ('Accelerate BLAS', 'Optional vendor implementation; internal threading and allocation are opaque.'),
}
SANITIZER = PublicSanitizer(root=ROOT)


def public(value, sanitizer=None):
    sanitizer = SANITIZER if sanitizer is None else sanitizer
    if isinstance(value, dict): return {key: public(item, sanitizer) for key,item in value.items()}
    if isinstance(value, (list,tuple)): return [public(item, sanitizer) for item in value]
    return sanitizer.text(value) if isinstance(value,str) else value


def execute(command, stdin=None, timeout=120):
    # Applied only to child processes; the exported schema never captures environment.
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS='1')
    result = subprocess.run(list(map(str,command)), input=stdin, text=True, capture_output=True, env=env, timeout=timeout, cwd=ROOT)
    if result.returncode:
        raise RuntimeError(SANITIZER.text((result.stderr or result.stdout)[-3000:]))
    return result.stdout


def build(blas=False, sanitize=False):
    compiler = shutil.which('clang')
    if not compiler: raise RuntimeError('Clang is required')
    if blas and platform.system() != 'Darwin': raise ValueError('Accelerate requires macOS; omit --blas')
    BUILD.mkdir(parents=True,exist_ok=True)
    target = BUILD / ('sanitized' if sanitize else 'runner')
    flags = ['-std=c11','-O1' if sanitize else '-O3','-g','-ffp-contract=off','-fno-vectorize','-fno-slp-vectorize',
             '-Wall','-Wextra','-Werror','-fno-omit-frame-pointer',
             '-ffile-prefix-map='+str(ROOT)+'=.','-fdebug-prefix-map='+str(ROOT)+'=.']
    if sanitize: flags += ['-fsanitize=address,undefined']
    if blas: flags += ['-DHAVE_ACCELERATE','-DACCELERATE_NEW_LAPACK','-framework','Accelerate']
    command = [compiler,*flags,HERE/'src/runner.c',HERE/'src/kernels.c','-o',target]
    execute(command)
    variants=['ijk','ikj','blocked']
    if platform.machine().lower() in ('arm64','aarch64'): variants.append('neon')
    if blas: variants.append('blas')
    return target, variants, {'compiler':execute([compiler,'--version']).splitlines()[0], 'command':SANITIZER.command(command)}


def reference(a,b,m,k,n):
    av,bv=list(map(Fraction,a)),list(map(Fraction,b))
    return [sum((av[i*k+p]*bv[p*n+j] for p in range(k)),Fraction()) for i in range(m) for j in range(n)]


def verify(runner, variants):
    rng=random.Random(41027)
    shapes=[(0,0,0),(0,3,5),(3,0,5),(3,5,0),(1,1,1),(2,3,4),(5,7,3),
            (31,33,35),(33,31,17),(3,65,67)]
    cases=[]
    for m,k,n in shapes:
        a=[rng.randrange(-32,33)/16 for _ in range(m*k)]
        b=[rng.randrange(-32,33)/16 for _ in range(k*n)]
        cases.append(('dyadic %dx%dx%d'%(m,k,n),m,k,n,a,b))
    cases += [('decimal cancellation',2,3,2,[.1,1e12,-1e12,-.3,.2,.7],[.3,-.2,.1,.4,.1,.4]),
              ('mixed scales',2,3,2,[1e100,1e-100,3,-1e100,2,1e-100],[1e-100,1e-100,1e100,-1e100,2,3])]
    records=[]
    for name,m,k,n,a,b in cases:
        wanted=reference(a,b,m,k,n)
        magnitude=[sum(abs(Fraction(a[i*k+p])*Fraction(b[p*n+j])) for p in range(k)) for i in range(m) for j in range(n)]
        for variant in variants:
            result=json.loads(execute([runner,'check',variant,m,k,n], ' '.join(format(x,'.17g') for x in a+b)))
            if (result['rows'],result['cols'])!=(m,n) or len(result['values'])!=m*n: raise AssertionError('incorrect output shape')
            for actual,expected,bound in zip(result['values'],wanted,magnitude):
                # Absolute forward error scaled by sum |a*b|, robust to cancellation.
                tolerance=32*max(1,k)*sys.float_info.epsilon*float(bound)+1e-300
                if not math.isfinite(actual) or abs(Fraction(actual)-expected)>Fraction(tolerance): raise AssertionError(variant+': '+name)
            records.append({'variant':variant,'case':name,'passed':True})
    for variant in variants:
        for mode in ('reused','allocated'):
            invoke(runner,variant,(3,5,7),3,19,mode)
    return records


def numerator(index,seed): return (index*17+seed*13)%33-16


def expected_checksum(shape,seed):
    m,k,n=shape
    if not m*n:return 0.0
    total=0
    for index in (0,m*n//2,m*n-1):
        i,j=divmod(index,n)
        total+=sum(numerator(i*k+p,seed)*numerator(p*n+j,seed+1) for p in range(k))
    return total/256


def invoke(runner,variant,shape,iterations,seed,mode):
    for attempt in range(5):
        result=json.loads(execute([runner,'bench',variant,*shape,iterations,seed,mode]))
        if result.get('elapsed_ns') != 0: break
        if attempt == 4 or iterations > 625000: raise AssertionError('timer resolution insufficient after bounded retries')
        iterations *= 16
    if result.get('iterations')!=iterations or result.get('checksum')!=expected_checksum(shape,seed)*iterations or result.get('full_result_verified') is not True:
        raise AssertionError('invalid benchmark checksum or verification')
    elapsed=result.get('elapsed_ns')
    if not isinstance(elapsed,int) or elapsed<=0:raise AssertionError('invalid elapsed time')
    allocations=iterations if mode=='allocated' else 0
    if result.get('output_allocations')!=allocations or result.get('requested_output_bytes')!=allocations*shape[0]*shape[2]*8:
        raise AssertionError('allocation accounting differs from requested mode')
    result['ns_per_op']=elapsed/iterations
    return result


def collect_profile(runner,variant,iterations,output):
    profile={'variant':variant,'status':'unavailable','chronological':False,'shape':[256,256,256],'mode':'reused'}
    if platform.system()!='Darwin':
        return dict(profile,note='No native sampler configured on this platform.')
    raw=BUILD/(variant+'-sample.txt')
    raw.unlink(missing_ok=True)
    process=subprocess.Popen([str(runner),'bench',variant,'256','256','256',str(iterations),'17','reused'],cwd=ROOT,
        env=dict(os.environ,VECLIB_MAXIMUM_THREADS='1'),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        sampled=subprocess.run(['/usr/bin/sample',str(process.pid),'2','1','-file',str(raw)],capture_output=True,text=True,timeout=15)
        stdout,stderr=process.communicate(timeout=45)
        if process.returncode:raise RuntimeError(stderr)
        result=json.loads(stdout)
        if result.get('checksum')!=expected_checksum((256,256,256),17)*iterations or result.get('full_result_verified') is not True:
            raise AssertionError('profile result verification failed')
        if sampled.returncode:raise RuntimeError(sampled.stderr or sampled.stdout)
        text=sanitize_macos_sample(raw.read_text(),SANITIZER)
        stacks,weights=parse_macos_sample(text)
        if not stacks or not all(math.isfinite(w) and w>0 for w in weights):raise ValueError('no usable samples')
        destination=output/'profiles'/(variant+'.txt');destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(text)
        profile.update(status='available',stacks=stacks,weights=weights,raw_file='profiles/'+variant+'.txt',iterations=iterations,
            note='Separate process; main-thread aggregate stack counts × nominal 1 ms. Includes startup, warmup and verification. Not a chronological trace; worker threads are not measured.')
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        profile['note']=SANITIZER.text(str(error))[-1200:]
    finally:
        if process.poll() is None:
            process.terminate()
            try:process.communicate(timeout=3)
            except subprocess.TimeoutExpired:process.kill();process.communicate()
        raw.unlink(missing_ok=True)
    return public(profile)


def render(data,output,*,published=False):
    template=(HERE/'report.html').read_text()
    sanitizer = PublicSanitizer.for_published() if published else SANITIZER
    encoded=json.dumps(public(data, sanitizer),allow_nan=False).replace('<','\\u003c')
    (output/'index.html').write_text(apply_report_design(template.replace('__DATA__',encoded)))


def shape(text):
    try:values=tuple(int(v) for v in text.split('x'))
    except ValueError:raise argparse.ArgumentTypeError('use MxKxN')
    if len(values)!=3 or any(v<1 or v>1024 for v in values):raise argparse.ArgumentTypeError('three dimensions from 1 to 1024 required')
    return values


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--blas',action='store_true')
    parser.add_argument('--profiles',action='store_true')
    parser.add_argument('--verify-only',action='store_true')
    parser.add_argument('--sanitize',action='store_true',help='address/undefined sanitizers; correctness only')
    parser.add_argument('--samples',type=int,default=5)
    parser.add_argument('--shapes',type=shape,nargs='+',default=list(map(shape,['32x32x32','128x128x128','256x256x256','512x512x512','96x257x65'])))
    parser.add_argument('--output',type=Path,default=HERE/'results')
    parser.add_argument('--render',type=Path,help='render existing JSON without building or measuring')
    args=parser.parse_args()
    if args.render:
        args.output.mkdir(parents=True,exist_ok=True);render(json.loads(args.render.read_text()),args.output);return
    if args.samples<3 or args.samples>30:parser.error('samples must be between 3 and 30')
    if args.sanitize and not args.verify_only:parser.error('--sanitize requires --verify-only')
    runner,variants,build_info=build(args.blas,args.sanitize)
    checks=verify(runner,variants);print('%d exact/rational correctness checks passed'%len(checks),flush=True)
    if args.verify_only:return
    output=args.output;output.mkdir(parents=True,exist_ok=True)
    records=[];rng=random.Random(9211)
    for dimensions in args.shapes:
        for mode in ('reused','allocated'):
            calibrated={}
            for variant in variants:
                probe=invoke(runner,variant,dimensions,1,17,mode)
                iterations=max(1,min(10000000,int(50_000_000/probe['ns_per_op'])))
                probe=invoke(runner,variant,dimensions,iterations,17,mode)
                calibrated[variant]=max(1,min(10000000,int(probe['iterations']*50_000_000/probe['elapsed_ns'])))
            samples={variant:[] for variant in variants}
            for round_index in range(args.samples):
                order=variants[:];rng.shuffle(order)
                for position,variant in enumerate(order):
                    sample=invoke(runner,variant,dimensions,calibrated[variant],17,mode)
                    sample.update(round=round_index,order=position);samples[variant].append(sample)
            for variant in variants:
                times=[s['ns_per_op'] for s in samples[variant]];median=statistics.median(times)
                records.append(dict(variant=variant,shape=dimensions,mode=mode,median_ns=median,mad_ns=statistics.median(abs(t-median) for t in times),min_ns=min(times),max_ns=max(times),samples=samples[variant]))
            print('%s %s measured'%('x'.join(map(str,dimensions)),mode),flush=True)
    profiles=[]
    if args.profiles:
        for variant in variants:
            probe=invoke(runner,variant,(256,256,256),5,17,'reused')
            iterations=max(1,min(10000000,int(6e9/probe['ns_per_op'])))
            profiles.append(collect_profile(runner,variant,iterations,output))
            print('Profile',variant,profiles[-1]['status'],flush=True)
    # Disassembly is evidence about emitted instructions, not a timing source.
    disassembly=None
    if shutil.which('xcrun'):
        text=execute(['xcrun','llvm-objdump','--disassemble','--no-show-raw-insn',runner])
        disassembly='disassembly.txt';(output/disassembly).write_text(SANITIZER.text(text))
    source_hash=hashlib.sha256()
    for path in sorted((HERE/'src').glob('*')):
        source_hash.update(path.name.encode());source_hash.update(path.read_bytes())
    source_hash.update((HERE/'run.py').read_bytes())
    data=public(dict(schema_version=1,created_at=datetime.now(timezone.utc).isoformat(),
        revision=execute(['git','rev-parse','HEAD']).strip(),dirty=bool(execute(['git','status','--porcelain']).strip()),source_sha256=source_hash.hexdigest(),
        machine=dict(os=platform.system(),release=platform.release(),architecture=platform.machine(),logical_cpus=os.cpu_count()),
        build=build_info,variants=[dict(id=v,label=VARIANTS[v][0],hypothesis=VARIANTS[v][1]) for v in variants],
        checks=checks,results=records,profiles=profiles,disassembly=disassembly,
        methodology=dict(samples=args.samples,seed=17,order_seed=9211,tile=32,blas_threads_requested=1 if args.blas else None,
            timing='Fresh process per sample; three warmup calls; calibrated ~50 ms batches; monotonic timer. Initialization, kernel, checksum and optional output allocation/free included. Startup and full-result verification excluded.',
            allocations='Output malloc calls and requested bytes are counted by the harness. Reused mode has zero output allocations inside timing; both modes have three setup buffers. No handwritten kernel allocates. BLAS internal allocation is unknown; counts are not resident memory or allocator overhead.',
            controls='One native language/runtime; all scalar auto-vectorization and FP contraction disabled; explicit NEON uses separate multiply/add. BLAS is a vendor baseline with independent compiler, reduction, threading, and internal allocation choices.',
            limitations='Same input buffers reused: warm-cache steady state, not streaming. Five processes are not five machines. Serial shuffled order reduces but does not eliminate frequency, thermal, and scheduling effects. Tiles are fixed at 32, not tuned. No cache misses or bandwidth measured; locality explanations are hypotheses. Profiles are separate runs, not timing evidence.')))
    (output/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');render(data,output)

if __name__=='__main__':main()
