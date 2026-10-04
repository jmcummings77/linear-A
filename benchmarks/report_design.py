"""Inline the shared report design so saved HTML remains portable offline."""
from pathlib import Path
import html as markup
import json
import re

HERE = Path(__file__).resolve().parent


def apply_report_design(html):
    """Apply common tokens, navigation, controls and theme without network assets."""
    # Source templates may mention a report whose first capture is not published yet.
    catalog_path = HERE / 'pages/catalog.json'
    if catalog_path.exists():
        for item in json.loads(catalog_path.read_text())['items']:
            if item['kind'] == 'benchmark' and not (HERE / 'reports' / item['id'] / 'index.html').exists():
                html = re.sub(r'<a\b(?=[^>]*href="\.\./' + re.escape(item['id']) + r'/")[^>]*>.*?</a>', '', html, flags=re.S)
    html = re.sub(r'<style id="report-design">.*?</style>', '', html, flags=re.S)
    css = (HERE / 'report-ui.css').read_text(encoding='utf-8') + '\n' + (HERE / 'report-navigation.css').read_text(encoding='utf-8')
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
    html = re.sub(r'<script id="report-controls">.*?</script>', '', html, flags=re.S)
    controls = '<script id="report-controls">\n' + (HERE / 'report-ui.js').read_text(encoding='utf-8') + '\n</script>'
    if '</body>' in html:
        html = html.replace('</body>', controls + '\n</body>', 1)
    else:
        html += '\n' + controls
    html = re.sub(r'<script id="report-navigation">.*?</script>', '', html, flags=re.S)
    navigation = '<script id="report-navigation">\n' + (HERE / 'report-navigation.js').read_text(encoding='utf-8') + '\n</script>'
    if '</body>' in html:
        html = html.replace('</body>', navigation + '\n</body>', 1)
    else:
        html += '\n' + navigation
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


def directory_catalog():
    """One topic identity with explicit destinations for each reader task."""
    catalog = json.loads((HERE / 'pages/catalog.json').read_text())
    modes = {'results': ('View results', 'chart'), 'demo': ('Try demo', 'play'),
             'evidence': ('Inspect evidence', 'check')}
    groups = []
    options = []
    for group in catalog['groups']:
        entries = []
        for item in catalog['items']:
            if item['group'] != group['id']:
                continue
            metadata = ''
            if item['kind'] == 'benchmark':
                source = HERE / 'reports' / item['id']
                # A new algorithm can land before its first captured report.
                if not (source / 'results.json').exists() or not (source / 'index.html').exists():
                    continue
                data = json.loads((source / 'results.json').read_text())
                date = str(data.get('created_at', 'Undated run')).split('T')[0]
                count = sum(row.get('status') == 'passed' for row in data.get('results', []))
                metadata = f'{date} · {implementation_summary(item["id"], data)} · {count} timing results'
            search = ' '.join([item['title'], group['label'], *item['questions'].values()])
            actions, questions = [], []
            for mode, href in item['links'].items():
                label, icon = modes[mode]
                actions.append(f'<a class="entry-action" data-entry-view="{mode}" data-report-icon="{icon}" aria-label="{label}: {markup.escape(item["title"], quote=True)}" href="{markup.escape(href, quote=True)}">{label}</a>')
                questions.append(f'<p data-entry-question="{mode}">{markup.escape(item["questions"][mode])}</p>')
            title = markup.escape(item['title'])
            meta = f'<span class="run-meta">{markup.escape(metadata)}</span>' if metadata else ''
            entries.append(f'<article class="report-entry" data-entry="{item["id"]}" data-search="{markup.escape(search, quote=True)}" aria-labelledby="topic-{item["id"]}"><div><h4 id="topic-{item["id"]}">{title}</h4>{"".join(questions)}{meta}</div><div class="entry-actions">{"".join(actions)}</div></article>')
        if entries:
            options.append(f'<option value="{group["id"]}">{markup.escape(group["label"])}</option>')
            groups.append(f'<section class="report-family" id="{group["id"]}" data-topic="{group["id"]}"><h3>{markup.escape(group["label"])}</h3><p>{markup.escape(group["description"])}</p><div class="report-list">{"".join(entries)}</div></section>')
    return ''.join(groups), ''.join(options)


def render_directory(destination=None):
    """Render the report directory with dates and counts from saved measurements."""
    template = (HERE / 'pages/template.html').read_text(encoding='utf-8')
    replacements = {}
    if '__REPORT_CATALOG__' in template:
        catalog, options = directory_catalog()
        template = template.replace('__REPORT_CATALOG__', catalog).replace('__TOPIC_OPTIONS__', options)
        template = template.replace('__DIRECTORY_SCRIPT__', (HERE / 'pages/directory.js').read_text())
    for marker in set(re.findall(r'__([A-Z0-9_]+)_META__', template)):
        name = marker.lower().replace('_', '-')
        data = json.loads((HERE / 'reports' / name / 'results.json').read_text())
        date = str(data.get('created_at', 'Undated run')).split('T')[0]
        count = sum(row.get('status') == 'passed' for row in data.get('results', []))
        label = f'{date} · {implementation_summary(name, data)} · {count} timing results'
        replacements['__' + marker + '_META__'] = markup.escape(label)
    html = re.sub(r'__[A-Z0-9_]+_META__', lambda match: replacements[match[0]], template)
    Path(destination or HERE / 'pages/index.html').write_text(apply_report_design(html))
