const $ = id => document.getElementById(id);
const data = JSON.parse($('data').textContent), names = JSON.parse($('names').textContent);
const add = (parent, tag, text, className) => {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  parent.append(node);
  return node;
};
const fmt = x => x < 1000 ? x.toFixed(1) + ' ns' : x < 1e6 ? (x / 1000).toFixed(2) + ' µs' : (x / 1e6).toFixed(2) + ' ms';
$('run-context').textContent = ((data.metadata.timestamp_utc || 'Saved run').split('T')[0]) + ' · ARM64 · scalar float64 · ' + data.samples + ' samples per row';
$('passed').textContent = data.correctness_checks;
$('samples').textContent = data.samples;
const sizes = [...new Set(data.timings.map(r => r.n))], select = $('size');
for (const n of sizes) { const option = add(select, 'option', n.toLocaleString()); option.value = n; }
select.value = sizes[Math.min(2, sizes.length - 1)];

function draw() {
  const n = Number(select.value), rows = data.timings.filter(r => r.n === n).sort((a, b) => a.median_ns - b.median_ns);
  const base = rows.find(r => r.implementation === 'machine-code')?.median_ns;
  const hasReference = Number.isFinite(base) && base > 0;
  $('working').textContent = (n * 32 * 16 / 1048576).toFixed(2) + ' MiB input working set';
  $('rows').replaceChildren();
  for (const r of rows) {
    const tr = add($('rows'), 'tr', undefined, r.implementation === 'machine-code' ? 'raw' : '');
    for (const value of [names[r.implementation] || r.implementation, fmt(r.median_ns),
      fmt(Math.min(...r.ns_per_call)) + ' – ' + fmt(Math.max(...r.ns_per_call)), hasReference ? (r.median_ns / base).toFixed(2) + '×' : '—']) add(tr, 'td', value);
  }
  const chart = $('latency-chart'); chart.replaceChildren();
  if (!hasReference) {
    $('latency-note').textContent = 'Relative latency needs a positive machine-code timing. Exact measurements remain available below.';
    return;
  }
  const positive = rows.flatMap(r => [r.median_ns, ...r.ns_per_call]).filter(v => Number.isFinite(v) && v > 0).map(v => v / base);
  let low = Math.floor(Math.log10(Math.min(1, ...positive))), high = Math.ceil(Math.log10(Math.max(1, ...positive)));
  if (low === high) { low--; high++; }
  const position = ratio => 100 * (Math.log10(ratio) - low) / (high - low);
  const ticks = Array.from({length: high - low + 1}, (_, i) => 10 ** (low + i));
  const tickLabel = ratio => ratio >= .001 && ratio < 100000 ? Number(ratio.toPrecision(4)) + '×' : ratio.toExponential(0) + '×';
  const axis = add(chart, 'div', undefined, 'latency-row latency-axis');
  add(axis, 'span', 'Time ÷ machine code', 'latency-label');
  const axisTrack = add(axis, 'div', undefined, 'latency-track axis-track');
  ticks.forEach((tick, i) => {
    const label = add(axisTrack, 'span', tickLabel(tick), 'axis-tick' + (i === 0 ? ' first' : i === ticks.length - 1 ? ' last' : ''));
    label.style.left = position(tick) + '%';
  });
  for (const r of rows) {
    const row = add(chart, 'div', undefined, 'latency-row' + (r.implementation === 'machine-code' ? ' raw' : ''));
    const label = add(row, 'div', undefined, 'latency-label');
    add(label, 'span', names[r.implementation] || r.implementation);
    add(label, 'span', `${fmt(r.median_ns)} · ${(r.median_ns / base).toFixed(2)}×`, 'latency-value');
    const track = add(row, 'div', undefined, 'latency-track');
    for (const tick of ticks) { const grid = add(track, 'span', undefined, 'latency-grid'); grid.style.left = position(tick) + '%'; }
    const reference = add(track, 'span', undefined, 'latency-reference'); reference.style.left = position(1) + '%';
    const samples = r.ns_per_call.filter(v => Number.isFinite(v) && v > 0);
    if (samples.length) {
      const whisker = add(track, 'span', undefined, 'latency-range');
      const left = position(Math.min(...samples) / base), right = position(Math.max(...samples) / base);
      whisker.style.left = left + '%'; whisker.style.width = (right - left) + '%';
    }
    if (r.median_ns > 0) { const dot = add(track, 'span', undefined, 'latency-dot'); dot.style.left = position(r.median_ns / base) + '%'; }
    else add(track, 'span', '0 · outside log scale', 'latency-zero');
  }
  const zero = rows.some(r => r.median_ns === 0 || r.ns_per_call.some(v => v === 0));
  $('latency-note').textContent = 'Dots: medians. Whiskers: observed sample ranges. Orange dashed line: hand-encoded machine code (1×). Left is faster; right is slower. All rows share one logarithmic ratio axis.' + (zero ? ' Zero timings are excluded from the logarithmic positions and remain in the table.' : '');
}
select.onchange = draw;
draw();
$('code').textContent = data.words.map((word, i) => String(i * 4).padStart(2, '0') + '  ' + word).join('\n');
$('metadata').textContent = JSON.stringify(data.metadata, null, 2);
