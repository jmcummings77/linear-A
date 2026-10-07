import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {resolve} from 'node:path';
import {fileURLToPath} from 'node:url';
import {spawnSync} from 'node:child_process';
import test from 'node:test';
import vm from 'node:vm';
const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const sandbox = {module: {exports: {}}};
vm.runInNewContext(readFileSync(resolve(root, 'experiments/nonnormal-gmres/report.js'), 'utf8'), sandbox);
const {format, selectedRuns, chartMarkup, rangeMarkup, historyRows} = sandbox.module.exports;
const generated = spawnSync('python3', ['-c', `import importlib.util,json; s=importlib.util.spec_from_file_location('experiment','experiments/nonnormal-gmres/run.py'); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); print(json.dumps(m.experiment()))`], {cwd: root, encoding: 'utf8'});
assert.equal(generated.status, 0, generated.stderr);
const data = JSON.parse(generated.stdout);

test('recorded controls select every coupling and scaling without changing evidence', () => {
  const before = JSON.stringify(data);
  for (const coupling of [0, .25, .75, 2, 4]) for (const jacobi of [false, true]) {
    const runs = selectedRuns(data, coupling, jacobi);
    assert.deepEqual(Array.from(runs, run => run.restart), [1, 2]);
    assert.ok(runs.every(run => run.s === coupling && run.right_jacobi === jacobi));
    const svg = chartMarkup(runs, data.max_iterations, data.rtol, true, 1, 0);
    assert.equal((svg.match(/<circle/g) || []).length, runs.reduce((sum, run) => sum + run.history.length, 0));
    assert.equal((svg.match(/<polyline/g) || []).length, 4);
    for (const run of runs) {
      const rows = historyRows(run, run.iterations);
      assert.equal((rows.match(/<tr/g) || []).length, run.history.length);
      assert.equal((rows.match(/class="selected-frame"/g) || []).length, 1);
      assert.ok(rows.includes(format(run.history.at(-1).true_residual)));
    }
  }
  assert.equal(JSON.stringify(data), before);
});

test('zero residuals have finite display coordinates and preserve zero values', () => {
  const runs = selectedRuns(data, 0, false);
  const svg = chartMarkup(runs, 24, data.rtol, false, 2, 1);
  assert.ok(!/NaN|Infinity/.test(svg));
  assert.ok(svg.includes('true residual 0 (displayed at chart floor)'));
  assert.equal((svg.match(/<polyline/g) || []).length, 2);
  assert.equal(format(0), '0');
  assert.equal(format(100000), '100000');
  assert.equal(format(1e-10), '1e-10');
});

test('numerical range preserves a point disk and finite geometry at the largest coupling', () => {
  for (const coupling of [0, .25, .75, 2, 4]) {
    const svg = rangeMarkup(coupling);
    assert.ok(svg.includes(`radius ${coupling}`));
    assert.ok(!/NaN|Infinity/.test(svg));
  }
  assert.ok(rangeMarkup(0).includes('r="0"'));
});


test('short captures have distinct iteration ticks', () => {
  const svg = chartMarkup(selectedRuns(data, .25, false), 2, data.rtol, false, 1, 0);
  const ticks = [...svg.matchAll(/<text x="[^"]+" y="318" text-anchor="middle">(.*?)<\/text>/g)].map(match => match[1]);
  assert.deepEqual(ticks, ['0', '1', '2']);
});
