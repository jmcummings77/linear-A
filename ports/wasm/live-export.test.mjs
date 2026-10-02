import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import vm from 'node:vm';

const source = await readFile(new URL('../../benchmarks/live-report.mjs', import.meta.url), 'utf8');

function page(userAgent) {
  const elements = new Map();
  const element = () => ({ value: '', textContent: '', append() {}, replaceChildren() {},
    setAttribute() {}, click() {}, remove() {} });
  const get = id => {
    if (!elements.has(id)) elements.set(id, element());
    return elements.get(id);
  };
  get('live-data').textContent = JSON.stringify({ available: true, sha256: 'ab'.repeat(32),
    byte_length: 80, fixtures: [{ name: 'rectangular multiplication' }], worker_source: '' });
  get('live-operation').value = 'multiply';
  get('live-samples').value = '3';
  let worker;
  class Worker {
    constructor() { worker = this; }
    postMessage() {}
    terminate() {}
  }
  vm.runInNewContext(source, { document: { getElementById: get, createElement: element,
    body: element() }, navigator: { userAgent }, Worker, WebAssembly: {}, Blob,
    URL: { createObjectURL: () => 'blob:test', revokeObjectURL() {} },
    setTimeout: () => 1, clearTimeout() {} });
  return { get, deliver(message) { worker.onmessage({ data: message }); } };
}

test('live exports retain measurements but omit private worker metadata and full user agents', () => {
  const ui = page('Mozilla/5.0 Chrome/145.0.1.2 private-machine /Users/private-person API_TOKEN=secret-value');
  ui.get('live-run').onclick();
  const samples = [{ elapsed_ms: 0.25, ns_per_op: 125000, checksum: 17.25,
    environment: { API_TOKEN: 'secret-value' } }];
  ui.deliver({ type: 'done', kind: 'benchmark', checks: { passed: 0, total: 1,
    checks: [{ name: '/Users/private-person/test', passed: false,
      error: 'API_TOKEN=secret-value private-machine', environment: { token: 'secret-value' } }],
    hostname: 'private-machine' }, result: { operation: 'multiply', size: 48, seed: 17,
    iterations: 2, samples, median_ns: 125000, mad_ns: 0, min_ns: 125000,
    max_ns: 125000, checksum_verified: true, hostname: 'private-machine',
    environment: { API_TOKEN: 'secret-value' } } });
  ui.get('live-export').onclick();
  const serialized = ui.get('live-export-json').textContent;
  const payload = JSON.parse(serialized);
  assert.deepEqual(payload.browser, { name: 'Chrome', version: '145.0.1.2' });
  assert.equal(payload.benchmarks[0].median_ns, 125000);
  assert.deepEqual(payload.benchmarks[0].samples, [{ elapsed_ms: 0.25, ns_per_op: 125000, checksum: 17.25 }]);
  assert.equal(payload.checks.checks[0].name, 'rectangular multiplication');
  assert.equal(payload.checks.checks[0].passed, false);
  assert.match(payload.checks.completed_at, /^\d{4}-\d\d-/);
  for (const value of ['private-person', 'private-machine', 'secret-value', 'API_TOKEN', '/Users/']) {
    assert.ok(!serialized.includes(value), value);
  }
});

test('unrecognized browser metadata is not copied into exports', () => {
  const ui = page('Custom /home/private-person machine-uuid=private-machine');
  ui.get('live-check').onclick();
  ui.deliver({ type: 'done', kind: 'checks', result: { passed: 1, total: 1,
    checks: [{ name: 'rectangular multiplication', passed: true }] } });
  ui.get('live-export').onclick();
  const payload = JSON.parse(ui.get('live-export-json').textContent);
  assert.deepEqual(payload.browser, { name: 'Unknown' });
  assert.deepEqual(payload.benchmarks, []);
});
