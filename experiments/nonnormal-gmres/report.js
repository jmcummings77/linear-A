/* Inspect recorded trajectories; never recompute solver output in the browser. */
(function () {
  'use strict';
  const floor = 1e-16;
  const format = value => value === 0 ? '0' : String(Number(Number(value).toPrecision(6)));
  const stop = run => `${run.iterations} steps, ${run.reason.replaceAll('_', ' ')}, residual ${format(run.history.at(-1).true_residual)}`;
  function selectedRuns(data, coupling, jacobi) {
    return [1, 2].map(restart => data.runs.find(run => run.s === coupling && run.right_jacobi === jacobi && run.restart === restart));
  }
  function chartMarkup(runs, limit, rtol, showEnvelope, selectedRestart, frame, width = 640) {
    const left = 64, right = width - 16, top = 24, bottom = 296;
    const x = k => left + k / limit * (right - left);
    const y = value => bottom - (Math.log10(Math.max(value, floor)) + 16) / 16 * (bottom - top);
    const line = (x1, y1, x2, y2, extra = '') => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" ${extra}/>`;
    let markup = '<title>Recorded true residuals for restart lengths 1 and 2</title>';
    for (const exponent of [0, -4, -8, -12, -16]) {
      const py = y(10 ** exponent);
      markup += line(left, py, right, py, 'stroke="var(--line)"');
      markup += `<text x="${left - 10}" y="${py + 4}" text-anchor="end">1e${exponent}</text>`;
    }
    const ticks = [...new Set(Array.from({length: 5}, (_, i) => Math.round(i * limit / 4)))];
    for (const k of ticks) {
      markup += `<text x="${x(k)}" y="${bottom + 22}" text-anchor="middle">${k}</text>`;
    }
    markup += `<text x="${(left + right) / 2}" y="342" text-anchor="middle">Iteration</text>`;
    markup += line(left, y(rtol), right, y(rtol), 'stroke="var(--muted)" stroke-dasharray="3 5" opacity="0.65"');
    markup += `<text x="${right - 3}" y="${y(rtol) - 6}" text-anchor="end">Stopping tolerance ${format(rtol)}</text>`;
    for (const run of runs) {
      const color = run.restart === 1 ? 'var(--green)' : 'var(--orange)';
      const dash = run.restart === 2 ? 'stroke-dasharray="8 5"' : '';
      const points = field => run.history.map(h => `${x(h.iteration)},${y(h[field])}`).join(' ');
      if (showEnvelope) markup += `<polyline points="${points('cycle_envelope')}" fill="none" stroke="${color}" stroke-width="1.5" stroke-dasharray="2 5" opacity="0.7"/>`;
      markup += `<polyline points="${points('true_residual')}" fill="none" stroke="${color}" stroke-width="2.6" ${dash}/>`;
      for (const h of run.history) {
        const selected = run.restart === selectedRestart && h.iteration === frame;
        markup += `<circle cx="${x(h.iteration)}" cy="${y(h.true_residual)}" r="${selected ? 5 : 2.4}" fill="${color}" stroke="var(--card)" stroke-width="${selected ? 2 : 0}"><title>Restart ${run.restart}, iteration ${h.iteration}, true residual ${h.true_residual}${h.true_residual < floor ? ' (displayed at chart floor)' : ''}</title></circle>`;
      }
    }
    return markup;
  }
  function rangeMarkup(s) {
    const radius = Math.abs(s), extent = Math.max(1.5, radius + 0.4), scale = 118 / extent;
    const origin = 160 - scale;
    return `<title>Numerical range: disk centered at 1 with radius ${radius}</title>
      <line x1="18" y1="165" x2="304" y2="165" stroke="var(--line)"/>
      <line x1="${origin}" y1="24" x2="${origin}" y2="300" stroke="var(--line)"/>
      <circle cx="160" cy="165" r="${radius * scale}" fill="var(--tint)" fill-opacity="0.6" stroke="var(--green)" stroke-width="2"/>
      <line x1="160" y1="165" x2="${160 + radius * scale}" y2="165" stroke="var(--green)" stroke-dasharray="3 3"/>
      <circle cx="160" cy="165" r="4" fill="var(--green)"/><circle cx="${origin}" cy="165" r="3" fill="var(--orange)"/>
      <text x="160" y="185" text-anchor="middle">1 (both eigenvalues)</text><text x="${origin}" y="151" text-anchor="middle">0</text>
      <text x="298" y="151" text-anchor="end">Re z</text><text x="${origin + 8}" y="32">Im z</text>
      <text x="160" y="327" text-anchor="middle">Radius |s| = ${radius}</text>`;
  }
  function historyRows(run, selected) {
    return run.history.map(h => `<tr${h.iteration === selected ? ' class="selected-frame"' : ''}><th scope="row">${h.iteration}</th>${[...h.x, h.true_residual, h.estimated_residual, h.cycle_envelope, h.rounding_allowance].map(value => `<td>${format(value)}</td>`).join('')}</tr>`).join('');
  }
  if (typeof module !== 'undefined' && module.exports) module.exports = {format, selectedRuns, chartMarkup, rangeMarkup, historyRows};
  if (typeof document === 'undefined') return;
  const el = id => document.getElementById(id);
  const data = JSON.parse(el('data').textContent);
  function update(reset = false) {
    const coupling = Number(el('coupling').value), jacobi = el('scaling').value === 'jacobi';
    const runs = selectedRuns(data, coupling, jacobi), restart = Number(el('inspect-restart').value);
    const run = runs[restart - 1];
    const frame = Math.min(reset ? 0 : Number(el('frame').value), run.iterations);
    el('frame').max = String(run.iterations); el('frame').value = String(frame);
    el('frame-number').textContent = `Iteration ${frame} of ${run.iterations}`;
    const width = Math.max(280, Math.round(el('residual-plot').getBoundingClientRect().width));
    el('residual-plot').setAttribute('viewBox', `0 0 ${width} 350`);
    el('residual-plot').innerHTML = chartMarkup(runs, data.max_iterations, data.rtol, el('show-envelope').checked, restart, frame, width);
    el('range-plot').innerHTML = rangeMarkup(coupling);
    el('range-note').textContent = coupling >= 1 ? 'This disk contains the origin. The candidate numerical-range envelope gives no guaranteed decrease.' : 'The disk excludes the origin. Its radius describes information absent from the unchanged eigenvalues.';
    el('trajectory-status').textContent = `Restart 1: ${stop(runs[0])}. Restart 2: ${stop(runs[1])}.`;
    const h = run.history[frame];
    el('frame-readout').textContent = `x = (${format(h.x[0])}, ${format(h.x[1])}); true residual ${format(h.true_residual)}; cycle starts at ${h.cycle_start}, local step ${h.local_step}.`;
    el('history-rows').innerHTML = historyRows(run, frame);
  }
  for (const id of ['coupling', 'scaling', 'inspect-restart']) el(id).addEventListener('change', () => update(true));
  el('show-envelope').addEventListener('change', () => update());
  el('frame').addEventListener('input', () => update());
  window.addEventListener('resize', () => update());
  update();
  el('trajectory-controls').hidden = false;
  el('frame-controls').hidden = false;
})();
