/** Chart-local language selection; recorded measurements and exports stay intact. */
(() => {
  const languageId = id => id.replace(/-(?:lu|cholesky|cofactor)$/, '');
  const languages = new Map();
  for (const implementation of passed) {
    const id = languageId(implementation.id);
    if (!languages.has(id)) languages.set(id, implementation.name.split(' / ')[0]);
  }
  const selected = new Set(languages.keys()), checkboxes = [];
  const dropdown = $('timing-languages'), summary = $('timing-language-summary');
  for (const [id, name] of languages) {
    const label = el('label'), checkbox = el('input');
    checkbox.type = 'checkbox';
    checkbox.checked = true;
    checkbox.dataset.language = id;
    checkbox.onchange = () => {
      if (checkbox.checked) selected.add(id); else selected.delete(id);
      renderChart();
    };
    label.append(checkbox, el('span', name));
    $('timing-language-options').append(label);
    checkboxes.push(checkbox);
  }
  function selectAll(checked) {
    selected.clear();
    for (const checkbox of checkboxes) {
      checkbox.checked = checked;
      if (checked) selected.add(checkbox.dataset.language);
    }
    renderChart();
  }
  $('timing-language-all').onclick = () => selectAll(true);
  $('timing-language-none').onclick = () => selectAll(false);
  document.addEventListener('click', event => {
    if (!dropdown.contains(event.target)) dropdown.open = false;
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape' && dropdown.open) {
      dropdown.open = false;
      summary.focus();
      event.preventDefault();
    }
  });

  for (const op of [...new Set(data.results.map(row => row.operation))]) addOption($('operation'), op, op);
  for (const implementation of passed) addOption($('baseline'), implementation.id, implementation.name);
  if (passed.some(implementation => implementation.id === 'csharp')) $('baseline').value = 'csharp';

  function sizes() {
    const old = $('size').value;
    $('size').replaceChildren();
    const sizes = [...new Set(data.results.filter(row => row.operation === $('operation').value).map(row => row.size))].sort((a, b) => a - b);
    for (const size of sizes) addOption($('size'), size, $('operation').value === 'cross' ? '3D vectors' : `${size} × ${size}`);
    if ([...$('size').options].some(option => option.value === old)) $('size').value = old;
    renderChart();
  }

  function renderChart() {
    const matching = valid.filter(row => row.operation === $('operation').value && row.size === Number($('size').value));
    const rows = matching.filter(row => selected.has(languageId(row.implementation))).sort((a, b) => a.median_ns - b.median_ns);
    summary.textContent = `Languages (${selected.size}/${languages.size})`;
    const previous = $('baseline').value, available = new Set(rows.map(row => row.implementation));
    $('baseline').replaceChildren();
    for (const implementation of passed.filter(item => available.has(item.id))) addOption($('baseline'), implementation.id, implementation.name);
    if (available.has(previous)) $('baseline').value = previous;
    else if (available.has(previous + '-lu')) $('baseline').value = previous + '-lu';
    $('baseline').disabled = rows.length === 0;
    const base = rows.find(row => row.implementation === $('baseline').value);
    $('chart').replaceChildren();
    $('samples').replaceChildren();
    $('chart-note').textContent = '';
    if (!rows.length) {
      const message = selected.size === 0 ? 'Select at least one language to show timings.'
        : matching.length ? 'No successful timing samples for the selected languages.'
          : 'No successful timing samples for this selection.';
      $('chart').append(el('div', message, 'empty'));
      return;
    }
    const max = Math.max(...rows.map(row => row.median_ns));
    for (const result of rows) {
      const row = el('div', undefined, 'bar-row');
      row.dataset.base = result.implementation === $('baseline').value;
      const track = el('div', undefined, 'bar-track'), bar = el('div', undefined, 'bar');
      bar.style.width = (max > 0 ? 100 * result.median_ns / max : 0) + '%';
      track.append(bar);
      const ratio = base && result.median_ns > 0 ? `${speed(base.median_ns / result.median_ns)}× speed` : '—';
      row.append(el('span', names[result.implementation]), track, el('span', time(result.median_ns), 'timing'), el('span', ratio, 'ratio'));
      $('chart').append(row);
    }
    $('chart-note').textContent = `${rows.length}/${matching.length} implementations shown · scale fits shown results · ${rows[0].samples.length} samples per implementation · ratio = baseline median / implementation median · output allocation included`;
    const header = el('tr');
    for (const title of ['Implementation', 'Iterations / sample', 'Median', 'MAD', 'Range', 'Samples']) header.append(el('th', title));
    $('samples').append(header);
    for (const result of rows) {
      const row = el('tr');
      for (const value of [names[result.implementation], result.iterations, time(result.median_ns), time(result.mad_ns), `${time(result.min_ns)} – ${time(result.max_ns)}`, result.samples.map(sample => time(sample.ns_per_op)).join(', ')]) row.append(el('td', value));
      $('samples').append(row);
    }
  }
  $('operation').onchange = sizes;
  $('size').onchange = renderChart;
  $('baseline').onchange = renderChart;
  sizes();
})();
