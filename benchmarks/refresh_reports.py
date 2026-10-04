#!/usr/bin/env python3
"""Refresh saved report presentation, preserving measurements and executable bundles."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil

import compare
import report
import sparse
import ordering_bench
import cholesky_bench
import amd_bench
from report_design import render_directory

ROOT = Path(__file__).resolve().parents[1]


def embedded(path, identity):
    content = path.read_text(encoding='utf-8')
    pattern = r'<script\b(?=[^>]*\bid="' + re.escape(identity) + r'")(?=[^>]*\btype="application/json")[^>]*>(.*?)</script>'
    match = re.search(pattern, content, re.S)
    if not match:
        raise ValueError(f'{path.relative_to(ROOT)} has no saved {identity} bundle; refusing to replace it with a new build')
    return json.loads(match[1])


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def measurements():
    paths = list((ROOT / 'benchmarks/reports').glob('*/results.json'))
    paths += list((ROOT / 'experiments').glob('*/results/results.json'))
    return {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def refresh():
    before = measurements()
    count = 0
    for path in sorted((ROOT / 'benchmarks/reports').glob('*/index.html')):
        data = embedded(path, 'data')
        if path.parent.name == 'amd':
            live = embedded(path, 'live')
            amd_bench.render(data, path, live_override=live)
            assert embedded(path, 'live') == live
        elif path.parent.name == 'cholesky':
            live = embedded(path, 'live')
            cholesky_bench.render(data, path, live_override=live)
            assert embedded(path, 'live') == live
        elif path.parent.name == 'ordering':
            live = embedded(path, 'live')
            ordering_bench.render(data, path, live_override=live)
            assert embedded(path, 'live') == live
        elif path.parent.name in ('sparse', 'gmres', 'ilu'):
            live = embedded(path, 'live')
            sparse.render(data, path, live_override=live)
            assert embedded(path, 'live') == live
        else:
            live = embedded(path, 'live-data')
            report.render(data, path, live_override=live)
            assert embedded(path, 'live-data') == live
        assert embedded(path, 'data') == data, f'Recorded data changed in {path}'
        count += 1
    path = ROOT / 'applications/index.html'
    bundle = embedded(path, 'bundle')
    module('playground_build', 'applications/build.py').render(bundle, path)
    assert embedded(path, 'bundle') == bundle
    path = ROOT / 'benchmarks/comparison/index.html'
    snapshots = embedded(path, 'snapshots')
    compare.render([(item['label'], item['data']) for item in snapshots], path)
    assert embedded(path, 'snapshots') == snapshots
    for name in ('matmul-locality', 'machine-code-dot'):
        directory = ROOT / 'experiments' / name / 'results'
        study = module(name.replace('-', '_'), f'experiments/{name}/run.py')
        data = json.loads((directory / 'results.json').read_text())
        study.render(data, directory if name == 'matmul-locality' else directory / 'index.html')
        assert embedded(directory / 'index.html', 'data') == data
    render_directory()
    if before != measurements():
        raise AssertionError('A measurement JSON file changed during presentation refresh')
    print(f'Refreshed {count + 5} HTML pages. Measurement files and embedded live bundles are unchanged.')


def assemble_preview(destination):
    """Match the published directory layout, including relative report links."""
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copytree(ROOT / 'benchmarks/reports', destination, dirs_exist_ok=True)
    shutil.copy2(ROOT / 'benchmarks/pages/index.html', destination / 'index.html')
    for name in ('matmul-locality', 'machine-code-dot'):
        shutil.copytree(ROOT / 'experiments' / name / 'results', destination / name, dirs_exist_ok=True)
    for name, source in [('applications', 'applications/index.html'), ('compare', 'benchmarks/comparison/index.html')]:
        (destination / name).mkdir(exist_ok=True)
        shutil.copy2(ROOT / source, destination / name / 'index.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir', type=Path, help='Also assemble a local directory matching the published site')
    args = parser.parse_args()
    refresh()
    if args.preview_dir:
        assemble_preview(args.preview_dir)
        print(f'Preview directory: {args.preview_dir}')
