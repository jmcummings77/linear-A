/** A single topic catalog, filtered by the reader's question. No network requests. */
(() => {
  const $ = id => document.getElementById(id);
  const root = $('report-library');
  if (!root) return;
  const tasks = {
    results: {title: 'Compare recorded performance', noun: 'report', note: 'These are saved measurements. Compare implementations or methods for the stated workload; use Compare two runs for changes between snapshots.'},
    demo: {title: 'Understand an algorithm', noun: 'interactive topic', note: 'These demos run WebAssembly on your device. Their live calculations are separate from the recorded benchmarks.'},
    evidence: {title: 'Inspect the evidence', noun: 'report', note: 'Go directly to the measurement method, verification, or limitations. Original data downloads remain available within each report.'},
  };
  const buttons = [...root.querySelectorAll('[data-directory-mode]')];
  const groups = [...root.querySelectorAll('[data-topic]')];
  const entries = [...root.querySelectorAll('[data-entry]')];
  const normalize = value => value.toLocaleLowerCase().normalize('NFKD').replace(/[\u0300-\u036f]/g, '');
  const legacyViews = {explore:'demo', benchmarks:'results', verification:'evidence'};
  const fromHash = () => {
    const hash = location.hash.slice(1);
    return Object.hasOwn(legacyViews, hash) ? legacyViews[hash] : Object.hasOwn(tasks, hash) ? hash : 'results';
  };
  let mode = fromHash();
  function render() {
    const query = normalize($('library-search').value.trim()), topic = $('library-topic').value;
    let count = 0;
    for (const group of groups) {
      let visible = 0;
      for (const entry of group.querySelectorAll('[data-entry]')) {
        const actions = [...entry.querySelectorAll('[data-entry-view]')];
        const supported = actions.some(action => action.dataset.entryView === mode);
        const match = supported && (topic === 'all' || group.dataset.topic === topic) && query.split(/\s+/).every(word => normalize(entry.dataset.search).includes(word));
        entry.hidden = !match;
        for (const action of actions) { action.hidden = action.dataset.entryView !== mode; action.dataset.selected = String(!action.hidden); }
        for (const question of entry.querySelectorAll('[data-entry-question]')) question.hidden = question.dataset.entryQuestion !== mode;
        const metadata = entry.querySelector('.run-meta');
        if (metadata) metadata.hidden = mode === 'demo';
        if (match) { visible++; count++; }
      }
      group.hidden = visible === 0;
    }
    for (const button of buttons) button.setAttribute('aria-pressed', String(button.dataset.directoryMode === mode));
    $('library-title').textContent = tasks[mode].title;
    $('task-context').textContent = tasks[mode].note;
    $('library-count').textContent = `${count} ${tasks[mode].noun}${count === 1 ? '' : 's'}${query || topic !== 'all' ? (count === 1 ? ' matches your filters' : ' match your filters') : ' available'}`;
    $('library-empty').hidden = count > 0;
    $('compare').hidden = mode !== 'results';
    $('clear-filters').disabled = !query && topic === 'all';
  }
  function clear() { $('library-search').value = ''; $('library-topic').value = 'all'; render(); $('library-search').focus(); }
  for (const button of buttons) button.addEventListener('click', () => {
    mode = button.dataset.directoryMode;
    if (location.hash !== '#' + mode) history.pushState(null, '', '#' + mode);
    render();
  });
  $('library-search').addEventListener('input', render);
  $('library-topic').addEventListener('change', render);
  $('clear-filters').addEventListener('click', clear);
  $('empty-clear').addEventListener('click', clear);
  window.addEventListener('hashchange', () => { mode = fromHash(); render(); });
  $('task-choices').hidden = false;
  $('library-filters').hidden = false;
  render();
})();
