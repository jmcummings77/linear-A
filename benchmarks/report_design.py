"""Inline the shared report design so saved HTML remains portable offline."""
from pathlib import Path
import html as markup
import json
import re

HERE = Path(__file__).resolve().parent


def apply_report_design(html):
    """Apply common tokens, navigation, controls and theme without network assets."""
    html = re.sub(r'<style id="report-design">.*?</style>', '', html, flags=re.S)
    css = (HERE / 'report-ui.css').read_text(encoding='utf-8')
    theme = (HERE / 'report-theme.js').read_text(encoding='utf-8')
    style = '<style id="report-design">\n' + css + '\n</style>'
    # Put shared rules after template styles, including compact implicit-head pages.
    position = html.rfind('</style>')
    if position >= 0:
        position += len('</style>')
        html = html[:position] + '\n' + style + html[position:]
    else:
        html = html.replace('</head>', style + '\n</head>', 1)
    if theme not in html:
        html = html.replace('<style', '<script id="report-appearance">' + theme + '</script>\n<style', 1)
    return html


def implementation_summary(name, data):
    """Distinguish algorithm variants and external references from project ports."""
    passed = [item for item in data.get('implementations', []) if item.get('status') == 'passed']
    if name == 'determinants':
        ports = {item.get('name', item.get('id', '')).split(' / ')[0] for item in passed}
        return f'{len(passed)} verified algorithm variants across {len(ports)} ports'
    references = [item for item in passed if item.get('id') == 'suitesparse']
    if references:
        return f'{len(passed) - len(references)} verified ports + SuiteSparse reference'
    return f'{len(passed)} verified implementations'


def render_directory(destination=None):
    """Render the report directory with dates and counts from saved measurements."""
    template = (HERE / 'pages/template.html').read_text(encoding='utf-8')
    replacements = {}
    for marker in set(re.findall(r'__([A-Z0-9_]+)_META__', template)):
        name = marker.lower().replace('_', '-')
        data = json.loads((HERE / 'reports' / name / 'results.json').read_text())
        date = str(data.get('created_at', 'Undated run')).split('T')[0]
        count = sum(row.get('status') == 'passed' for row in data.get('results', []))
        label = f'{date} · {implementation_summary(name, data)} · {count} timing results'
        replacements['__' + marker + '_META__'] = markup.escape(label)
    html = re.sub(r'__[A-Z0-9_]+_META__', lambda match: replacements[match[0]], template)
    Path(destination or HERE / 'pages/index.html').write_text(apply_report_design(html))
