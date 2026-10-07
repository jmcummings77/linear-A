"""Render a portable report from saved experiment data without rerunning solvers."""
import html
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from report_design import apply_report_design


def number(value):
    return '0' if value == 0 else format(value, '.6g')


def render_html(data):
    """Keep the saved values intact; all dynamic content comes from this payload."""
    runs = data['runs']
    rows = []
    for row in runs:
        rows.append('<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%s</td><td>%s</td></tr>' % (
            number(row['s']), 'Right Jacobi' if row['right_jacobi'] else 'None', row['restart'],
            row['iterations'], html.escape(row['reason'].replace('_', ' ')),
            number(row['history'][-1]['true_residual'])))
    summary = []
    for polynomial in data['polynomials']:
        pair = [next(r for r in runs if r['s'] == polynomial['s'] and r['restart'] == restart
                     and not r['right_jacobi']) for restart in (1, 2)]
        one, two = pair
        summary.append('<tr><th scope="row">%s</th><td>%d · %s</td><td>%s</td><td>%d</td><td>%s</td></tr>' % (
            number(polynomial['s']), one['iterations'], html.escape(one['reason'].replace('_', ' ')),
            number(one['history'][-1]['true_residual']), two['iterations'],
            number(two['history'][-1]['true_residual'])))
    polynomials = ''.join('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (
        number(p['s']), number(p['polynomial_norm']), number(p['polynomial_supremum']),
        '—' if not p['s'] else number(p['polynomial_norm']/p['polynomial_supremum']))
        for p in data['polynomials'])
    graph_rows = []
    for case in data.get('graphs', {}).get('cases', []):
        solves = case['solves']
        graph_rows.append('<tr><th scope="row">%s</th><td>%d</td><td>%d</td><td>%s</td><td>%s</td><td>%s / %s</td><td>%s</td></tr>' % (
            html.escape(case['name']), case['n'], case['nnz'], number(case['spectral_gap']),
            number(case['kappa2']), solves[0]['iterations'], solves[1]['iterations'],
            number(max(s['relative_true_residual'] for s in solves))))
    source = data['source']
    revision = html.escape(str(source['revision']))
    runtime = data.get('runtime', {})
    replacements = {
        'DATA': json.dumps(data, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026'),
        'SCRIPT': (HERE / 'report.js').read_text(),
        'SUMMARY_ROWS': ''.join(summary), 'RUN_ROWS': ''.join(rows), 'POLYNOMIAL_ROWS': polynomials,
        'GRAPH_ROWS': ''.join(graph_rows), 'REVISION': revision,
        'CLEAN_STATE': 'clean' if source['dirty'] is False else 'dirty or unavailable',
        'LIMIT': str(data['max_iterations']), 'RUN_COUNT': str(len(runs)),
        'COUPLING_COUNT': str(len(data['polynomials'])),
        'RESTART2_MAX': str(max(r['iterations'] for r in runs if r['restart'] == 2)),
        'RTOL': number(data['rtol']),
        'RUNTIME': html.escape(' · '.join(str(runtime.get(key, 'unspecified')) for key in ('python', 'system', 'architecture'))),
        'CREATED_AT': html.escape(str(data.get('created_at', 'not recorded'))),
    }
    template = (HERE / 'report.html').read_text()
    rendered = re.sub(r'@@([A-Z0-9_]+)@@', lambda m: replacements[m[1]], template)
    return apply_report_design(rendered)
