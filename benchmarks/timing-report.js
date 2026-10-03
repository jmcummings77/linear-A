/** Chart-local language selection; recorded measurements and exports stay intact. */
(() => {
  const operationLabels = {
    add: 'Matrix addition', subtract: 'Matrix subtraction', scale: 'Scalar multiplication',
    transpose: 'Matrix transpose', multiply: 'Matrix multiplication', trace: 'Matrix trace',
    determinant: 'Determinant · general matrices', determinant_small: 'Determinant · small matrices',
    determinant_spd: 'Determinant · positive-definite matrices',
    eigen_symmetric: 'Eigenvalues & eigenvectors · symmetric', eigen_general: 'Eigenvalues & eigenvectors · general real',
    cross: 'Cross product · 3D vectors', rotation2d: 'Rotation matrix · 2D', rotation3d: 'Rotation matrix · 3D'
  };
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

  for (const op of [...new Set(data.results.map(row => row.operation))]) addOption($('operation'), op, operationLabels[op] || op.replaceAll('_', ' '));
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
    $('timing-axis').replaceChildren();
    $('timing-workload').textContent = '';
    if (!rows.length) {
      const message = selected.size === 0 ? 'Select at least one language to show timings.'
        : matching.length ? 'No successful timing samples for the selected languages.'
          : 'No successful timing samples for this selection.';
      $('chart').append(el('div', message, 'empty'));
      return;
    }
    const max = Math.max(...rows.map(row => row.median_ns));
    const logarithmic = $('timing-scale').value === 'log';
    const positive = rows.flatMap(row => [row.min_ns, row.median_ns, row.max_ns]).filter(value => Number.isFinite(value) && value > 0);
    const low = positive.length ? Math.log10(Math.min(...positive)) : 0;
    const high = positive.length ? Math.log10(Math.max(...positive)) : 1;
    const padding = Math.max((high - low) * .06, .1), domainLow = low - padding, domainHigh = high + padding;
    const position = value => Math.max(0, Math.min(100, 100 * (Math.log10(value) - domainLow) / (domainHigh - domainLow)));
    $('timing-workload').textContent = `${operationLabels[$('operation').value] || $('operation').value} · ${$('size').options[$('size').selectedIndex]?.textContent || ''}. Median per operation; lower is faster. MAD describes sample variation, not uncertainty in the mean.`;
    for (const result of rows) {
      const row = el('div', undefined, 'bar-row');
      row.dataset.base = result.implementation === $('baseline').value;
      const track = el('div', undefined, logarithmic ? 'bar-track log-track' : 'bar-track');
      track.setAttribute('aria-hidden', 'true');
      if (logarithmic) {
        if (result.median_ns > 0) {
          if (result.min_ns > 0 && result.max_ns > 0) {
            const range = el('span', undefined, 'timing-range');
            range.style.left = position(result.min_ns) + '%';
            range.style.width = (position(result.max_ns) - position(result.min_ns)) + '%';
            track.append(range);
          }
          const dot = el('span', undefined, 'timing-dot');
          dot.style.left = position(result.median_ns) + '%';
          track.append(dot);
        } else track.append(el('span', 'Zero is not on a log scale', 'small'));
      } else {
        const bar = el('div', undefined, 'bar');
        bar.style.width = (max > 0 ? 100 * result.median_ns / max : 0) + '%';
        track.append(bar);
      }
      const ratio = base && result.median_ns > 0 ? `${speed(base.median_ns / result.median_ns)}× baseline speed` : '—';
      const metric = el('span', undefined, 'timing');
      metric.append(el('span', time(result.median_ns), 'timing-value'));
      if (Number.isFinite(result.mad_ns)) metric.append(el('span', `MAD ${time(result.mad_ns)}`, 'timing-spread'));
      if (Number.isFinite(result.min_ns) && Number.isFinite(result.max_ns)) metric.append(el('span', `${time(result.min_ns)} – ${time(result.max_ns)}`, 'timing-spread'));
      row.append(el('span', names[result.implementation], 'bar-label'), track, metric, el('span', ratio, 'ratio'));
      $('chart').append(row);
    }
    const labels = el('div', undefined, 'timing-axis-labels');
    if (logarithmic) {
      labels.append(el('span', time(10 ** domainLow)), el('span', time(10 ** ((domainLow + domainHigh) / 2))), el('span', time(10 ** domainHigh)));
    } else labels.append(el('span', '0'), el('span', time(max)));
    $('timing-axis').append(labels);
    const sampleCounts = [...new Set(rows.map(row => row.samples.length))].sort((a, b) => a - b);
    $('chart-note').textContent = `${rows.length}/${matching.length} implementations shown · ${logarithmic ? 'logarithmic dots; whiskers show observed min–max' : 'linear bars from zero'} · scale fits shown results · ${sampleCounts.join('–')} samples per implementation · ratio = baseline median / implementation median · output allocation included`;
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
  $('timing-scale').onchange = renderChart;
  sizes();
})();
