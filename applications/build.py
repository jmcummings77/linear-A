#!/usr/bin/env python3
"""Build a portable numerical playground from a sanitized, verified WASM bundle."""
import hashlib
import json
from pathlib import Path
import re
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'benchmarks'))
from report import live_bundle, json_for_html
from report_design import apply_report_design
from run import build_one, verify

def render(bundle,destination=None):
    """Render the interface around an existing bundle without rebuilding kernels."""
    destination=Path(destination) if destination is not None else ROOT/'applications/index.html'
    replacements={'__BUNDLE__':json_for_html(bundle),'__APP__':(ROOT/'applications/app.mjs').read_text()}
    template=(ROOT/'applications/template.html').read_text()
    html=re.sub('|'.join(replacements),lambda m:replacements[m[0]],template)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(apply_report_design(html))

def main():
    implementation=build_one('wasm',no_build=True)
    if not verify(implementation):raise RuntimeError('WebAssembly shared checks failed')
    bundle=live_bundle({'implementations':[implementation]})
    if not bundle['available']:raise RuntimeError(bundle['reason'])
    # Keep only runtime assets needed by these applications; no measured report data.
    bundle={k:bundle[k] for k in ('module_source','wrapper_source','wasm_base64')}
    bundle['math_source']=(ROOT/'applications/math.mjs').read_text()
    bundle['worker_source']=(ROOT/'applications/worker.mjs').read_text()
    bundle['sha256']=hashlib.sha256(json.dumps(bundle,sort_keys=True).encode()).hexdigest()
    render(bundle)
    print('Generated applications/index.html after %d shared WASM checks'%len(implementation['checks']))
if __name__=='__main__':main()
