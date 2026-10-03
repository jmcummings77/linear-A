#!/usr/bin/env python3
"""Compile and run a bounded libFuzzer/ASan/UBSan campaign (or replay one input)."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',type=int,default=10000)
    parser.add_argument('--seed',type=int,default=20261002)
    parser.add_argument('--replay',type=Path)
    parser.add_argument('--smoke',action='store_true',help='deterministic sanitizer driver, without coverage feedback')
    parser.add_argument('--assembly',action='store_true',help='exercise ARM64 kernels with C validation')
    args=parser.parse_args()
    if not 1<=args.runs<=1000000:parser.error('runs must be 1..1000000')
    compiler=os.environ.get('CC') or shutil.which('clang')
    if not compiler:parser.error('Clang with libFuzzer is required')
    build=ROOT/'.build/fuzz';build.mkdir(parents=True,exist_ok=True)
    corpus=build/('assembly-corpus' if args.assembly else 'c-corpus');corpus.mkdir(exist_ok=True)
    # Seed dimensions, SIMD/tile tails, empty shapes and extreme byte values.
    for index,data in enumerate([bytes([0,0,0,0]),bytes([1,1,1,1,128,128]),bytes([5,5,5,5,0,255,127,129]),bytes(range(64))]):
        (corpus/('seed-%d'%index)).write_bytes(data)
    target=build/('assembly-fuzz' if args.assembly else 'matrix-fuzz')
    command=[compiler,'-std=c11','-O1','-g','-ffp-contract=off','-fno-omit-frame-pointer',
             '-fsanitize=address,undefined' if args.smoke else '-fsanitize=fuzzer,address,undefined','-fno-sanitize-recover=all',
             '-ffile-prefix-map='+str(ROOT)+'=.','-fdebug-prefix-map='+str(ROOT)+'=.',
             str(ROOT/'ports/c/matrix.c'),str(ROOT/'tests/fuzz/matrix_fuzz.c'),'-lm','-o',str(target)]
    if args.assembly:command += ['-DMATRIX_USE_ASM',str(ROOT/'ports/assembly/kernels.S')]
    if args.smoke:command.append(str(ROOT/'tests/fuzz/smoke.c'))
    if args.smoke and args.replay:parser.error('--smoke cannot replay libFuzzer artifacts')
    subprocess.run(command,check=True,cwd=ROOT)
    runtime=[str(target),'-seed='+str(args.seed),'-runs='+str(args.runs),'-max_len=128','-timeout=5',
             '-rss_limit_mb=1024','-artifact_prefix='+str(build)+'/']
    runtime.append(str(args.replay) if args.replay else str(corpus))
    # Raw sanitizer diagnostics stay local; they can contain runtime host paths.
    if args.smoke:runtime=[str(target),str(args.runs),str(args.seed)]
    return subprocess.run(runtime,cwd=ROOT,timeout=180).returncode

if __name__=='__main__':raise SystemExit(main())
