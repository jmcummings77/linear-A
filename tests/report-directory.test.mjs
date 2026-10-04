import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import { readFileSync } from 'node:fs';

const source = readFileSync(new URL('../benchmarks/pages/directory.js', import.meta.url), 'utf8');

function element({ id = '', dataset = {}, className = '', children = [], textContent = '', href } = {}) {
  const attributes = new Map(href ? [['href', href]] : []), events = new Map();
  const node = {
    id, dataset, className, children, textContent, hidden: false, disabled: false, value: '', focused: false,
    setAttribute: (name, value) => attributes.set(name, String(value)),
    getAttribute: name => attributes.get(name) ?? null,
    addEventListener(name, fn) { events.set(name, fn); },
    dispatch(name) { events.get(name)?.({ target: node }); },
    focus() { node.focused = true; },
    events,
    querySelectorAll(selector) {
      const attr = /^\[data-([a-z-]+)\]$/.exec(selector);
      const key = attr?.[1].replace(/-([a-z])/g, (_, letter) => letter.toUpperCase());
      assert.ok(key || selector === '.run-meta', `Unimplemented fixture selector: ${selector}`);
      const matches = item => key ? Object.hasOwn(item.dataset, key) : item.className.split(' ').includes('run-meta');
      return children.flatMap(child => [...(matches(child) ? [child] : []), ...child.querySelectorAll(selector)]);
    },
    querySelector(selector) { return this.querySelectorAll(selector)[0] ?? null; },
  };
  return node;
}

function entry(id, views, search, { metadata = true } = {}) {
  const actions = views.map(view => element({ dataset: { entryView: view }, textContent: `${view}: ${id}`, href: `${id}/#${view}` }));
  const questions = views.map(view => element({ dataset: { entryQuestion: view }, textContent: `${view} question for ${id}` }));
  const meta = metadata ? element({ className: 'run-meta', textContent: '2026-10-03 · 12 saved measurements' }) : null;
  return element({ id, dataset: { entry: id, search }, children: [...actions, ...questions, ...(meta ? [meta] : [])] });
}

function library({ hash = '', absent = false } = {}) {
  const entries = [
    entry('dense', ['results', 'demo', 'evidence'], 'Dense matrix performance multiplication mémoire'),
    entry('solver', ['results', 'demo', 'evidence'], 'Sparse solver preconditioning runtime memory'),
    entry('machine-code', ['results', 'evidence'], 'Machine code dot product runtime'),
    entry('applications', ['demo', 'evidence'], 'Applications rotations image compression', { metadata: false }),
  ];
  const groups = [
    element({ dataset: { topic: 'matrix' }, children: [entries[0]] }),
    element({ dataset: { topic: 'solvers' }, children: [entries[1]] }),
    element({ dataset: { topic: 'hardware' }, children: [entries[2]] }),
    element({ dataset: { topic: 'applications' }, children: [entries[3]] }),
  ];
  const buttons = ['results', 'demo', 'evidence'].map(mode => element({ id: mode, dataset: { directoryMode: mode } }));
  const root = element({ id: 'report-library', children: [...buttons, ...groups] });
  const ids = new Map([[root.id, root], ...entries.map(node => [node.id, node]), ...buttons.map(node => [node.id, node])]);
  for (const id of ['library-search', 'library-topic', 'library-title', 'task-context', 'library-count',
    'library-empty', 'compare', 'clear-filters', 'empty-clear', 'task-choices', 'library-filters']) ids.set(id, element({ id }));
  const get = id => ids.get(id);
  get('library-topic').value = 'all'; get('task-choices').hidden = get('library-filters').hidden = true;
  const location = { hash }, historyEntries = [hash], events = new Map();
  let cursor = 0, networkCalls = 0;
  const network = () => { networkCalls++; throw new Error('The portable directory must not use the network'); };
  const history = { pushState(state, unused, url) {
    assert.equal(state, null); assert.equal(unused, ''); assert.match(url, /^#/);
    historyEntries.splice(++cursor, Infinity, url); location.hash = url;
  } };
  vm.runInNewContext(source, {
    document: { getElementById: id => absent && id === 'report-library' ? null : get(id) },
    location, history, fetch: network, XMLHttpRequest: network,
    window: { addEventListener: (name, fn) => events.set(name, fn), fetch: network },
    navigator: { sendBeacon: network },
  });
  return {
    get, buttons, entries, groups, location, historyEntries,
    get networkCalls() { return networkCalls; },
    click: id => get(id).dispatch('click'),
    search(value) { get('library-search').value = value; get('library-search').dispatch('input'); },
    topic(value) { get('library-topic').value = value; get('library-topic').dispatch('change'); },
    hash(value) { location.hash = value; events.get('hashchange')?.(); },
    back() { assert.ok(cursor > 0); location.hash = historyEntries[--cursor]; events.get('hashchange')?.(); },
    visible() { return entries.filter(node => !node.hidden).map(node => node.id); },
  };
}

function assertView(page, mode) {
  assert.deepEqual(page.buttons.filter(button => button.getAttribute('aria-pressed') === 'true').map(button => button.id), [mode]);
  for (const entry of page.entries.filter(node => !node.hidden)) {
    const actions = entry.querySelectorAll('[data-entry-view]');
    assert.deepEqual(actions.filter(node => !node.hidden).map(node => node.dataset.entryView), [mode]);
    for (const action of actions) assert.equal(action.dataset.selected, String(!action.hidden));
    assert.deepEqual(entry.querySelectorAll('[data-entry-question]').filter(node => !node.hidden).map(node => node.dataset.entryQuestion), [mode]);
    const metadata = entry.querySelector('.run-meta');
    if (metadata) assert.equal(metadata.hidden, mode === 'demo');
  }
  assert.equal(page.get('compare').hidden, mode !== 'results');
}

test('each mode shows its supported topics, question, action, and recorded-data context', () => {
  const p = library();
  assert.deepEqual(p.visible(), ['dense', 'solver', 'machine-code']);
  assert.equal(p.get('library-title').textContent, 'Compare recorded performance');
  assert.equal(p.get('library-count').textContent, '3 reports available');
  assert.match(p.get('task-context').textContent, /saved measurements/);
  assert.equal(p.get('task-choices').hidden, false); assert.equal(p.get('library-filters').hidden, false);
  assertView(p, 'results');
  p.click('demo');
  assert.deepEqual(p.visible(), ['dense', 'solver', 'applications']);
  assert.equal(p.groups.find(group => group.dataset.topic === 'hardware').hidden, true);
  assert.equal(p.get('library-count').textContent, '3 interactive topics available');
  assert.match(p.get('task-context').textContent, /on your device/);
  assertView(p, 'demo');
  p.click('evidence');
  assert.deepEqual(p.visible(), ['dense', 'solver', 'machine-code', 'applications']);
  assert.equal(p.get('library-count').textContent, '4 reports available');
  assert.match(p.get('task-context').textContent, /verification/);
  assertView(p, 'evidence');
});

test('search matches every normalized word and combines with topic filtering', () => {
  const p = library();
  p.search('  MEMOIRE   matrix ');
  assert.deepEqual(p.visible(), ['dense']);
  assert.equal(p.get('library-count').textContent, '1 report matches your filters');
  p.search('runtime');
  assert.deepEqual(p.visible(), ['solver', 'machine-code']);
  p.topic('solvers');
  assert.deepEqual(p.visible(), ['solver']);
  assert.equal(p.groups.filter(group => !group.hidden).length, 1);
  p.click('demo');
  assert.equal(p.get('library-count').textContent, '1 interactive topic matches your filters');
  p.click('results');
  p.search('runtime rotations');
  assert.deepEqual(p.visible(), []);
  assert.equal(p.get('library-count').textContent, '0 reports match your filters');
  assert.equal(p.get('library-empty').hidden, false);
  assert.equal(p.get('clear-filters').disabled, false);
});

test('empty-state clear resets both filters, retains the chosen mode, and restores search focus', () => {
  const p = library({ hash: '#demo' });
  p.topic('hardware');
  assert.equal(p.get('library-empty').hidden, false, 'hardware has no interactive action');
  p.search('missing');
  p.click('empty-clear');
  assert.equal(p.get('library-search').value, ''); assert.equal(p.get('library-topic').value, 'all');
  assert.equal(p.get('library-search').focused, true);
  assert.equal(p.get('library-empty').hidden, true);
  assert.equal(p.get('clear-filters').disabled, true);
  assert.deepEqual(p.visible(), ['dense', 'solver', 'applications']);
  assert.equal(p.location.hash, '#demo'); assert.equal(p.historyEntries.length, 1);
  assertView(p, 'demo');
});

test('clear availability follows active filters and mode changes preserve a deliberate filter', () => {
  const p = library();
  assert.equal(p.get('clear-filters').disabled, true);
  p.search('   '); assert.equal(p.get('clear-filters').disabled, true);
  p.search('runtime'); p.topic('hardware'); p.click('demo');
  assert.equal(p.get('library-search').value, 'runtime'); assert.equal(p.get('library-topic').value, 'hardware');
  assert.equal(p.get('library-empty').hidden, false); assert.equal(p.get('clear-filters').disabled, false);
  p.click('clear-filters');
  assert.equal(p.get('clear-filters').disabled, true);
  assert.equal(p.get('library-search').focused, true);
  assertView(p, 'demo');
});

test('deep views, legacy hashes, and browser history restore the matching mode', () => {
  for (const [hash, mode] of [['#demo', 'demo'], ['#evidence', 'evidence'], ['#results', 'results'],
    ['#explore', 'demo'], ['#benchmarks', 'results'], ['#verification', 'evidence'], ['#unknown', 'results'],
    ['#constructor', 'results'], ['#toString', 'results'], ['#__proto__', 'results']]) {
    const p = library({ hash }); assertView(p, mode); assert.equal(p.historyEntries.length, 1);
  }
  const p = library({ hash: '#results' });
  p.click('demo'); p.click('demo'); p.click('evidence');
  assert.deepEqual(p.historyEntries, ['#results', '#demo', '#evidence']);
  p.back(); assertView(p, 'demo'); p.back(); assertView(p, 'results');
  p.hash('#verification'); assertView(p, 'evidence');
});

test('filtering preserves native report links and makes no external calls', () => {
  const p = library();
  const links = p.entries.flatMap(entry => entry.querySelectorAll('[data-entry-view]'));
  const hrefs = links.map(link => link.getAttribute('href'));
  p.click('demo'); p.search('rotations'); p.click('clear-filters'); p.click('evidence');
  assert.deepEqual(links.map(link => link.getAttribute('href')), hrefs);
  assert.ok(links.every(link => link.events.size === 0), 'the directory must not intercept native report navigation');
  assert.equal(p.networkCalls, 0);
  assert.doesNotThrow(() => library({ absent: true }), 'the shared asset is harmless outside its directory');
});
