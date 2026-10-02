#!/usr/bin/env python3
"""Isolated ARM64 scalar dot-product experiment. Python standard library only."""
import argparse
import html
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import struct
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / '.build'
SRC = ROOT / 'src'
RESULTS = ROOT / 'results'
LABELS = {'machine-code': 'Hand-encoded machine code', 'assembly': 'Assembly (identical instructions)',
          'c': 'C', 'cpp': 'C++', 'rust': 'Rust', 'go': 'Go', 'csharp': 'C#', 'fsharp': 'F#',
          'typescript': 'TypeScript / Node', 'python': 'Python', 'julia': 'Julia'}
ENV = dict(os.environ, DOTNET_CLI_HOME=str(BUILD / 'cli'), NUGET_PACKAGES=str(BUILD / 'nuget'),
           DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_GENERATE_ASPNET_CERTIFICATE='false',
           DOTNET_TieredCompilation='0', DOTNET_NOLOGO='1', GOCACHE=str(BUILD / 'go-cache'),
           GOPATH=str(BUILD / 'go'), GOTOOLCHAIN='local', JULIA_NUM_THREADS='1',
           JULIA_DEPOT_PATH=str(BUILD / 'julia-depot'))

def command(args, timeout=120):
    result = subprocess.run([str(x) for x in args], cwd=ROOT, env=ENV,
                            text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError('Command failed: ' + ' '.join(map(str, args)) + '\n' + result.stdout + result.stderr)
    return result.stdout.strip()

def executable(name, override=None):
    found = override or shutil.which(name)
    if not found:
        raise RuntimeError('Missing tool: ' + name)
    return str(Path(found).resolve())

def build(args):
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        raise RuntimeError('This experiment targets ARM64 macOS.')
    BUILD.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    for directory in ('cli', 'nuget', 'go-cache', 'go', 'julia-depot'):
        (BUILD / directory).mkdir(exist_ok=True)
    clang, cpp, rust, go, node = [executable(n) for n in ('clang', 'clang++', 'rustc', 'go', 'node')]
    dotnet = executable('dotnet', args.dotnet)
    julia = executable('julia', args.julia)
    versions = {key: command(cmd).splitlines()[0] for key, cmd in {
        'clang': [clang, '--version'], 'rust': [rust, '--version'], 'go': [go, 'version'],
        'node': [node, '--version'], 'python': [sys.executable, '--version'],
        'dotnet': [dotnet, '--version'], 'julia': [julia, '--version']}.items()}
    flags = ['-O3', '-ffp-contract=off', '-Wall', '-Wextra']
    if not args.no_build:
        command([clang, *flags, '-DRAW_MACHINE_CODE', SRC / 'native_runner.c', '-o', BUILD / 'machine-code'])
        command([clang, *flags, '-c', SRC / 'native_runner.c', '-o', BUILD / 'harness.o'])
        command([clang, *flags, '-c', SRC / 'dot.c', '-o', BUILD / 'c.o'])
        command([cpp, *flags, '-std=c++17', '-c', SRC / 'dot.cpp', '-o', BUILD / 'cpp.o'])
        command([clang, '-c', SRC / 'assembly.S', '-o', BUILD / 'assembly.o'])
        for key in ('c', 'cpp', 'assembly'):
            command([clang, BUILD / 'harness.o', BUILD / (key + '.o'), '-o', BUILD / key])
        command([rust, '-C', 'opt-level=3', '-C', 'target-cpu=native', SRC / 'dot.rs', '-o', BUILD / 'rust'])
        command([go, 'build', '-o', BUILD / 'go-runner', SRC / 'dot.go'])
        sdk_version = versions['dotnet']
        sdk_source = Path(dotnet).parent / 'sdk' / sdk_version / 'FSharp' / 'library-packs'
        for language, extension in (('csharp', 'csproj'), ('fsharp', 'fsproj')):
            cmd = [dotnet, 'build', ROOT / language / ('Dot.' + extension), '-c', 'Release',
                   '--disable-build-servers', '-m:1', '-p:NuGetAudit=false']
            # FSharp.Core ships with the SDK; this avoids an unnecessary network dependency.
            if sdk_source.is_dir():
                cmd.extend(['--source', sdk_source])
            command(cmd)
    repository_prefix = str(ROOT.parents[1]) + os.sep
    disassembly = command(['xcrun', 'llvm-objdump', '--disassemble', BUILD / 'assembly.o']).replace(repository_prefix, '')
    assembled_words = [int(x, 16) for x in re.findall(r'^\s*[0-9a-f]+:\s+([0-9a-f]{8})\s', disassembly, re.M)]
    manual_words = [int(x, 16) for x in re.findall(r'0x([0-9a-f]{8})', (SRC / 'raw_dot.h').read_text())]
    if assembled_words != manual_words or len(manual_words) != 9:
        raise RuntimeError('Assembly control does not match the hand-encoded 36-byte function.')
    (RESULTS / 'raw_dot.bin').write_bytes(struct.pack('<9I', *manual_words))
    (RESULTS / 'disassembly.txt').write_text(disassembly)
    (RESULTS / 'c-disassembly.txt').write_text(command(['xcrun', 'llvm-objdump', '--disassemble', BUILD / 'c.o']).replace(repository_prefix, ''))
    runners = {key: [str(BUILD / key)] for key in ('machine-code', 'assembly', 'c', 'cpp', 'rust')}
    runners.update({'go': [str(BUILD / 'go-runner')], 'typescript': [node, str(SRC / 'dot.ts')],
                    'python': [sys.executable, str(SRC / 'dot.py')],
                    'julia': [julia, '--startup-file=no', '--history-file=no', str(SRC / 'dot.jl')]})
    for key in ('csharp', 'fsharp'):
        runners[key] = [dotnet, str(ROOT / key / 'bin' / 'Release' / 'net10.0' / 'Dot.dll')]
    return runners, versions, manual_words

_REFERENCE = {}
def numerators(n, seed):
    key = n, seed
    if key not in _REFERENCE:
        # Exact integer oracle: every input is an integer divided by 16.
        _REFERENCE[key] = [sum(((i * 17 + seed * 13) % 101 - 50) *
                                  ((i * 29 + (seed + 1) * 7) % 103 - 51)
                                  for i in range(batch*n, (batch+1)*n)) for batch in range(32)]
    return _REFERENCE[key]

def invoke(runner, n, iterations, seed):
    data = json.loads(command(runner + [str(n), str(iterations), str(seed)]))
    parts = numerators(n, seed)
    expected = ((iterations // 32) * sum(parts) + sum(parts[:iterations % 32])) / 256
    if data['iterations'] != iterations or not math.isfinite(data['checksum']):
        raise RuntimeError('Invalid benchmark output: ' + repr(data))
    if data['checksum'] != expected:
        raise RuntimeError('Exact checksum mismatch: expected %r, got %r' % (expected, data['checksum']))
    if not math.isfinite(data['elapsed_ns']) or data['elapsed_ns'] < 0:
        raise RuntimeError('Invalid elapsed time')
    return data

def report(result):
    data = json.dumps(result, allow_nan=False).replace('<', '\\u003c')
    page = r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>36 bytes vs 10 languages</title><style>
:root{font-family:system-ui,sans-serif;color:#e9ecf4;background:#101522;color-scheme:dark}body{max-width:1180px;margin:40px auto;padding:0 24px}h1{font-size:42px;letter-spacing:-1.8px;margin-bottom:8px}p{color:#b6bed0;line-height:1.6}.tag{color:#7ce4bb;letter-spacing:2px;font-size:12px;text-transform:uppercase}.cards{display:flex;gap:12px;margin:25px 0;flex-wrap:wrap}.card{background:#1a2336;padding:16px 22px;border-radius:10px;flex:1}.card b{display:block;font-size:25px}.card span{font-size:13px;color:#b6bed0}.controls{display:flex;gap:15px;align-items:center;margin:28px 0}select{padding:8px;border-radius:6px}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:12px 10px;border-bottom:1px solid #29334a}th{font-size:12px;color:#aab5ce}.raw{color:#7ce4bb}.bar{height:10px;background:#719ff1;border-radius:3px;min-width:2px}.raw .bar{background:#7ce4bb}td small{display:block;color:#919cb5;margin-top:5px}code,pre{background:#1a2336;border-radius:6px}pre{padding:18px;overflow:auto;font-size:13px;line-height:1.6}details{margin:24px 0}a{color:#91b8ff}.note{border-left:3px solid #719ff1;padding-left:15px}footer{font-size:12px;color:#919cb5;margin:28px 0}.units{white-space:nowrap}</style>
<div class="tag">ARM64 · scalar float64 dot product</div><h1>36 bytes vs 10 languages</h1>
<p>A hand-encoded machine-code function, compared with the same scalar operation in every project language. Assembly is a byte-for-byte control.</p>
<div class="cards"><div class="card"><b>36 bytes</b><span>9 hand-encoded instructions</span></div><div class="card"><b id="passed"></b><span>exact correctness checks passed</span></div><div class="card"><b id="samples"></b><span>independent timing samples / row</span></div></div>
<div class="controls"><label>Vector length <select id="size"></select></label><span id="working"></span></div>
<table><thead><tr><th>Implementation</th><th>Median / call</th><th>Sample range</th><th>Time ÷ machine code</th><th>Relative latency (log scale)</th></tr></thead><tbody id="rows"></tbody></table>
<p class="note">Lower is faster. Ratios compare this kernel and these implementations on this machine; they are not general language rankings. Min–max ranges show run-to-run variation. Profiling, startup, ahead-of-time compilation, and input generation are excluded. JIT optimization may still occur during a short timed run. This was not run on an isolated host.</p>
<details open><summary>How this experiment works</summary><p>The kernel computes Σ a[i] × b[i], left to right, using separate double-precision multiply and add. No BLAS, SIMD intrinsics, parallelism, reassociation, or fast-math. 32 input pairs rotate to discourage loop-invariant elimination. Compilers may still auto-vectorize products, unroll loops, and optimize bounds checks. Inspection of Clang’s C output found eight-element blocks with SIMD multiplication and ordered scalar additions; the manual kernel has one scalar multiply per iteration. Julia uses @inbounds; raw machine code, C/C++, and assembly trust the caller's lengths. Other runtimes retain their normal bounds-checking semantics.</p><p>Each process warms the function 8–128 times. .NET tiered compilation is disabled to avoid mixing tiers in a short trial. Timed loops consume every result in a verified checksum. Multiple samples are run in shuffled, serial order. Calibration targets roughly 30 ms per sample, capped at ten million calls. Checksums use an independent exact-integer oracle.</p><p>The C program only provides allocation, timing, and Apple's JIT-memory interface. It copies the literal instruction words into executable memory and calls them; it does not compile the dot-product function. All native controls use an indirect function call. The assembly control is assembled separately, then its instruction words are checked against the manual encoding.</p></details>
<details><summary>Hand-encoded instructions</summary><pre id="code"></pre></details><details><summary>Machine, toolchains and run metadata</summary><pre id="metadata"></pre></details>
<footer><a href="results.json">Raw results and all samples</a> · <a href="raw_dot.bin">36-byte function</a> · <a href="disassembly.txt">Assembly control disassembly</a> · <a href="c-disassembly.txt">Clang C disassembly</a><p>© 2026 J.M. Cummings.</p></footer>
<script>const data=__DATA__;const names=__NAMES__;document.querySelector('#passed').textContent=data.correctness_checks;document.querySelector('#samples').textContent=data.samples;const sizes=[...new Set(data.timings.map(r=>r.n))];const select=document.querySelector('#size');for(const n of sizes){const o=document.createElement('option');o.value=n;o.textContent=n.toLocaleString();select.append(o)}select.value=sizes[Math.min(2,sizes.length-1)];const fmt=x=>x<1000?x.toFixed(1)+' ns':x<1e6?(x/1000).toFixed(2)+' µs':(x/1e6).toFixed(2)+' ms';function draw(){const n=+select.value;const rows=data.timings.filter(r=>r.n===n).sort((a,b)=>a.median_ns-b.median_ns);const base=rows.find(r=>r.implementation==='machine-code').median_ns;const min=Math.min(...rows.map(r=>r.median_ns));const max=Math.max(...rows.map(r=>r.median_ns));document.querySelector('#working').textContent=(n*32*16/1048576).toFixed(2)+' MiB input working set';const body=document.querySelector('#rows');body.replaceChildren();for(const r of rows){const tr=document.createElement('tr');if(r.implementation==='machine-code')tr.className='raw';const values=[names[r.implementation],fmt(r.median_ns),fmt(Math.min(...r.ns_per_call))+' – '+fmt(Math.max(...r.ns_per_call)),(r.median_ns/base).toFixed(2)+'×'];for(const v of values){const td=document.createElement('td');td.textContent=v;tr.append(td)}const td=document.createElement('td');const bar=document.createElement('div');bar.className='bar';bar.style.width=(10+90*Math.log(r.median_ns/min)/Math.max(1e-9,Math.log(max/min)))+'%';td.append(bar);tr.append(td);body.append(tr)}}select.onchange=draw;draw();document.querySelector('#code').textContent=data.words.map((w,i)=>String(i*4).padStart(2,'0')+'  '+w).join('\n');document.querySelector('#metadata').textContent=JSON.stringify(data.metadata,null,2);</script></html>'''
    page = page.replace('__DATA__', data).replace('__NAMES__', json.dumps(LABELS))
    (RESULTS / 'index.html').write_text(page)
    (RESULTS / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dotnet')
    parser.add_argument('--julia')
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--sizes', type=int, nargs='+', default=[16, 256, 4096, 65536])
    args = parser.parse_args()
    if args.samples < 3 or any(n < 1 or n > 1048576 for n in args.sizes): parser.error('At least 3 samples and lengths 1..1048576 required')
    runners, versions, words = build(args)
    checks = 0
    for key, runner in runners.items():
        for n, count, seed in ((0,33,7),(1,35,0),(3,67,7),(16,97,13),(255,65,991),(4096,33,7)):
            invoke(runner,n,count,seed); checks += 1
        print('Verified ' + LABELS[key], flush=True)
    timings = []
    rng = random.Random(7312)
    for n in args.sizes:
        iterations, samples = {}, {key: [] for key in runners}
        order = list(runners); rng.shuffle(order)
        for key in order:
            initial = invoke(runners[key],n,32,7)
            estimate = max(1, initial['elapsed_ns']/32)
            iterations[key] = max(32, min(10000000, int(30000000/estimate)))
        for sample in range(args.samples):
            rng.shuffle(order)
            for key in order:
                value = invoke(runners[key],n,iterations[key],7)
                checks += 1
                samples[key].append(value['elapsed_ns']/iterations[key])
            print('Length %d: sample %d/%d complete' % (n,sample+1,args.samples),flush=True)
        for key in runners:
            timings.append(dict(implementation=key,n=n,iterations=iterations[key],
                                ns_per_call=samples[key],median_ns=statistics.median(samples[key])))
        (RESULTS / 'timing-checkpoint.json').write_text(json.dumps(timings, indent=2))
    result = dict(metadata=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
                               os=platform.platform(),architecture=platform.machine(),versions=versions,
                               cpu=platform.processor() or platform.machine(),
                               compiler_flags='Native: -O3 -ffp-contract=off; Rust: opt-level=3 target-cpu=native',
                               dotnet_tiered_compilation=False,seed=7,order_seed=7312,
                               scope='Separate scalar dot kernels, not full matrix-library comparisons'),
                  correctness_checks=checks,samples=args.samples,
                  words=['0x%08x'%word for word in words],timings=timings)
    report(result)
    print('Report: '+str(RESULTS/'index.html'),flush=True)
    for n in args.sizes:
        print('\nLength %d'%n)
        for row in sorted((r for r in timings if r['n']==n),key=lambda r:r['median_ns']):
            print('%-32s %12.2f ns/call' % (LABELS[row['implementation']], row['median_ns']))

if __name__ == '__main__':
    main()
