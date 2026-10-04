import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import { readFileSync } from 'node:fs';

const source = readFileSync(new URL('../benchmarks/report-ui.js', import.meta.url), 'utf8');

function control(id, textContent, attrs = {}, tagName = 'BUTTON') {
  const attributes = new Map(Object.entries(attrs));
  return {
    id, tagName, textContent, disabled: false, onclick() {},
    getAttribute: name => attributes.get(name) ?? null,
    hasAttribute: name => attributes.has(name),
    setAttribute: (name, value) => attributes.set(name, String(value)),
    removeAttribute: name => attributes.delete(name),
    contains: node => node.parentControl === id,
  };
}

function page(controls, { ready = 'complete', withoutObserver = false } = {}) {
  const events = {}, observations = [];
  let callback, disconnected = false;
  const document = {
    readyState: ready,
    querySelectorAll: () => controls,
    addEventListener(name, listener, options) {
      assert.equal(name, 'DOMContentLoaded'); assert.equal(options.once, true); events[name] = listener;
    },
  };
  class MutationObserver {
    constructor(fn) { callback = fn; }
    observe(target, options) { observations.push({ target, options }); }
    disconnect() { disconnected = true; }
  }
  vm.runInNewContext(source, {
    document,
    window: { addEventListener: (name, listener) => { events[name] = listener; } },
    ...(withoutObserver ? {} : { MutationObserver }),
  });
  return {
    observations, events,
    changed(target) { callback?.([{ target }]); },
    get disconnected() { return disconnected; },
  };
}

test('common controls gain decorative icons without changing labels, handlers, or availability', () => {
  const play = control('play', 'Animate elimination'), next = control('next', 'Next pivot');
  play.disabled = true;
  const handler = play.onclick;
  const download = control('', 'Download recorded JSON', { download: '', href: 'results.json' }, 'A');
  const arbitrary = control('profile-sample-7', '7');
  const optedOut = control('reset', 'Reset zoom', { 'data-report-icon': 'none' });
  page([play, next, download, arbitrary, optedOut]);
  assert.equal(play.getAttribute('data-report-icon'), 'play');
  assert.equal(next.getAttribute('data-report-icon'), 'next');
  assert.equal(download.getAttribute('data-report-icon'), 'download');
  assert.equal(play.textContent, 'Animate elimination');
  assert.equal(play.disabled, true); assert.equal(play.onclick, handler);
  assert.equal(play.getAttribute('aria-pressed'), null, 'a changing-label action is not a toggle');
  assert.equal(arbitrary.getAttribute('data-report-icon'), null);
  assert.equal(optedOut.getAttribute('data-report-icon'), 'none');
});

test('playback icons follow textContent replacements and nested text changes without observing plots', () => {
  const play = control('play', 'Play'), geometry = control('geometry-play', 'Play transformation');
  const run = control('live-run', 'Run live benchmark');
  const p = page([play, geometry, run]);
  assert.deepEqual(p.observations.map(({ target }) => target.id), ['play', 'geometry-play']);
  for (const { options } of p.observations) {
    assert.equal(options.childList, true); assert.equal(options.characterData, true);
    assert.equal(options.subtree, true); assert.equal(options.attributes, undefined);
  }
  play.textContent = 'Pause'; p.changed(play);
  assert.equal(play.getAttribute('data-report-icon'), 'pause');
  assert.equal(play.textContent, 'Pause');
  geometry.textContent = 'Pause animation'; p.changed({ parentControl: 'geometry-play' });
  assert.equal(geometry.getAttribute('data-report-icon'), 'pause');
  play.textContent = 'Animate ordering'; p.changed(play);
  assert.equal(play.getAttribute('data-report-icon'), 'play');
  assert.equal(run.getAttribute('data-report-icon'), 'play');
});

test('settings has a visible label while retaining its dialog semantics', () => {
  const settings = control('geometry-settings-open', '', {
    'aria-label': 'Visualization settings', 'aria-haspopup': 'dialog',
    'aria-controls': 'geometry-settings', title: 'Visualization settings',
  });
  page([settings]);
  assert.equal(settings.textContent, 'Settings');
  assert.equal(settings.getAttribute('data-report-icon'), 'settings');
  assert.equal(settings.getAttribute('aria-controls'), 'geometry-settings');
  assert.equal(settings.getAttribute('aria-haspopup'), 'dialog');
  assert.equal(settings.getAttribute('aria-label'), 'Visualization settings');
  assert.equal(settings.getAttribute('title'), null);
});

test('explicit icons support new controls and cached-page navigation preserves live updates', () => {
  const button = control('new-player', 'Pause replay', { 'data-report-icon': 'play' });
  const link = control('', 'Report guide', { 'data-report-icon': 'info' }, 'A');
  const p = page([button, link], { ready: 'loading' });
  assert.equal(p.observations.length, 0);
  p.events.DOMContentLoaded();
  assert.equal(button.getAttribute('data-report-icon'), 'pause');
  assert.equal(link.getAttribute('data-report-icon'), 'info');
  p.events.pagehide({ persisted: true }); assert.equal(p.disconnected, false);
  button.textContent = 'Play replay'; p.changed(button);
  assert.equal(button.getAttribute('data-report-icon'), 'play');
  p.events.pagehide({ persisted: false }); assert.equal(p.disconnected, true);
});

test('lack of mutation observation leaves working native labels and initial icons', () => {
  const button = control('play', 'Pause');
  page([button], { withoutObserver: true });
  assert.equal(button.textContent, 'Pause');
  assert.equal(button.getAttribute('data-report-icon'), 'pause');
});
