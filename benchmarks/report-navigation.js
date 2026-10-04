/* Progressive report views: all original content remains available without JavaScript. */
(() => {
  'use strict';
  const boot = () => {
    const main = document.querySelector('main'), originalNav = main?.querySelector('.report-nav');
    if (!main || !originalNav || main.classList.contains('report-modes-enabled')) return;
    const win = window, find = id => document.getElementById(id);
    const decode = hash => { try { return decodeURIComponent(hash.replace(/^#/, '')); } catch { return hash.replace(/^#/, ''); } };
    const targetFor = hash => find(decode(hash));
    const internalTarget = link => {
      if (!link) return null;
      try {
        const url = new URL(link.getAttribute('href'), win.location.href);
        return url.origin === win.location.origin && url.pathname === win.location.pathname && url.search === win.location.search ? targetFor(url.hash) : null;
      } catch { return null; }
    };
    const originalLinks = [...originalNav.querySelectorAll('a[href]')];
    const labels = new Map(originalLinks.map(link => [internalTarget(link)?.id, link.textContent.trim()]));
    const topLevel = node => {
      while (node && node.parentElement !== main) node = node.parentElement;
      return node && node.tagName !== 'HEADER' ? node : null;
    };
    const available = node => !!node && !node.closest('[hidden]');
    let definitions;
    if (find('timings-panel') && find('method-panel')) definitions = {results:['timings-panel','profiles-panel'], demo:['live-panel','accuracy-panel','geometry-panel'], evidence:['correctness-panel','method-panel']};
    else if (find('runtime') && find('storage-section') && find('explore')) definitions = {results:['runtime','reuse-section','storage-section'], demo:['explore'], evidence:['method-section']};
    else if (find('results') && find('explore') && find('method-section')) definitions = {results:['results'], demo:['explore'], evidence:['method-section']};
    const back = main.querySelector('.report-topbar a');
    if (back && !find('report-library')) { back.textContent = 'Report library'; back.setAttribute('data-report-icon','back'); }
    // Playgrounds, comparisons, studies and the directory retain the complete reading order.
    if (!definitions || find('report-library')) {
      const links = originalLinks.filter(link => internalTarget(link));
      const mark = target => {
        const active = links.find(link => { const section = internalTarget(link); return section === target || section?.contains(target); }) || links[0];
        for (const link of links) link === active ? link.setAttribute('aria-current','location') : link.removeAttribute('aria-current');
      };
      mark(targetFor(win.location.hash));
      win.addEventListener('hashchange', () => mark(targetFor(win.location.hash)));
      if (typeof IntersectionObserver === 'function') {
        const visible = new Set();
        const observer = new IntersectionObserver(entries => {
          for (const entry of entries) entry.isIntersecting ? visible.add(entry.target) : visible.delete(entry.target);
          const first = links.map(internalTarget).find(target => visible.has(target));
          if (first) mark(first);
        }, {rootMargin:'-80px 0px -55% 0px'});
        links.map(internalTarget).forEach(target => observer.observe(target));
      }
      return;
    }
    const names = {results:'Results', demo:'Interactive demo', evidence:'Method & data'};
    const scopes = {
      results:'Saved measurements from the recorded run. These results do not change when you use the interactive demo.',
      demo:'Interactive calculations on your device. They are separate from the saved benchmark measurements.',
      evidence:'Correctness checks, measurement methods and run details for interpreting the saved results.'
    };
    const sections = Object.fromEntries(Object.entries(definitions).map(([mode, ids]) => [mode, ids.map(find).filter(node => topLevel(node) === node)]));
    const assignments = new Map();
    for (const [mode, nodes] of Object.entries(sections)) for (const node of nodes) assignments.set(node, mode);
    const assign = (node, mode) => { const block = topLevel(node); if (block) assignments.set(block, mode); };
    if (find('timings-panel')) {
      find('explore')?.querySelector('nav')?.classList.add('report-view-redundant-nav');
      assign(find('explore'), 'demo'); assign(find('stats'), 'results'); assign(find('saved-context'), 'results');
      main.querySelectorAll('.availability-note').forEach(node => assign(node, 'evidence'));
    } else {
      assign(find('overview'), 'results');
      main.querySelectorAll('.explore-label').forEach(node => assign(node, 'demo'));
    }
    const make = (tag, className, text) => {
      const node = document.createElement(tag);
      if (className) node.className = className;
      if (text) node.textContent = text;
      return node;
    };
    // Legacy observers keep their old detached anchors and cannot overwrite view selection.
    const nav = make('nav','report-view-nav'); nav.setAttribute('aria-label','Report views');
    const utilities = make('div','report-view-utilities');
    for (const link of originalLinks) if (!internalTarget(link)) utilities.append(link);
    const scope = make('p','report-view-scope'), localNav = make('nav','report-view-sections');
    localNav.setAttribute('aria-label','In this view');
    originalNav.replaceWith(nav, utilities, scope, localNav);
    // Related topics stay outside modes and out of the primary navigation hierarchy.
    const related = [...main.querySelectorAll('.related-reports, .ic0-related')];
    if (related.length) {
      const details = make('details','report-view-related'); details.append(make('summary','','Related reports'));
      for (const row of related) details.append(row);
      utilities.append(details);
    }
    utilities.hidden = !utilities.children.length;
    main.classList.add('report-modes-enabled');
    let activeMode, modes = [], modeLinks = new Map(), localLinks = [], resizeFrame;
    const redraw = () => { win.cancelAnimationFrame(resizeFrame); resizeFrame = win.requestAnimationFrame(() => win.dispatchEvent(new Event('resize'))); };
    const activeSections = mode => sections[mode].filter(available);
    const defaultMode = () => modes.includes('results') ? 'results' : modes.includes('evidence') ? 'evidence' : modes[0];
    const modeFor = target => { const mode = assignments.get(topLevel(target)); return modes.includes(mode) ? mode : null; };
    const localLabel = node => ({'live-panel':'Run a benchmark','accuracy-panel':'Sensitivity','geometry-panel':'Matrix geometry'}[node.id] || labels.get(node.id) || node.querySelector('h2')?.textContent.trim() || node.id);
    const markLocal = target => {
      for (const link of localLinks) {
        const section = internalTarget(link);
        if (section === target || section?.contains(target)) link.setAttribute('aria-current','location');
        else link.removeAttribute('aria-current');
      }
    };
    const activate = (mode, target) => {
      mode = modes.includes(mode) ? mode : defaultMode(); if (!mode) return;
      const changed = activeMode !== mode; activeMode = mode; main.setAttribute('data-report-view', mode);
      for (const [node, assigned] of assignments) {
        // Verification-only captures show their run summary with evidence, without an empty Results view.
        const actual = assigned === 'results' && !modes.includes('results') ? 'evidence' : assigned;
        node.classList.toggle('report-mode-inactive', actual !== mode);
      }
      for (const [key, link] of modeLinks) key === mode ? link.setAttribute('aria-current','page') : link.removeAttribute('aria-current');
      scope.textContent = scopes[mode];
      const nodes = activeSections(mode); localNav.replaceChildren(); localLinks = []; localNav.hidden = nodes.length < 2;
      if (nodes.length > 1) {
        localNav.append(make('span','','In this view:'));
        for (const node of nodes) { const link = make('a','',localLabel(node)); link.setAttribute('href',`#${node.id}`); localNav.append(link); localLinks.push(link); }
      }
      markLocal(target); if (changed) redraw();
    };
    const sync = () => {
      modes = Object.keys(sections).filter(mode => activeSections(mode).length); modeLinks = new Map(); nav.replaceChildren();
      for (const mode of modes) { const link = make('a','',names[mode]); link.setAttribute('data-report-icon', {results:'chart',demo:'play',evidence:'info'}[mode]); link.setAttribute('href',`#${activeSections(mode)[0].id}`); nav.append(link); modeLinks.set(mode,link); }
      const target = targetFor(win.location.hash); activate(modeFor(target) || activeMode || defaultMode(), target); redraw();
    };
    const openAncestors = target => { for (let node = target?.parentElement; node && node !== main; node = node.parentElement) if (node.tagName === 'DETAILS') node.open = true; };
    const revealHash = () => {
      const target = targetFor(win.location.hash); activate(modeFor(target) || defaultMode(), target); openAncestors(target);
      // Native history may scroll before a hidden target is revealed.
      if (target && available(target) && !target.closest('.report-mode-inactive')) win.requestAnimationFrame(() => target.scrollIntoView({block:'start',behavior:'auto'}));
    };
    sync(); if (win.location.hash) revealHash();
    win.addEventListener('hashchange', revealHash); win.addEventListener('popstate', revealHash);
    document.addEventListener('click', event => {
      if (event.defaultPrevented || event.button > 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
      const link = event.target.closest('a[href]'), target = internalTarget(link), mode = modeFor(target);
      // Reveal before native link scrolling, without replacing browser history or controls.
      if (mode) { activate(mode, target); openAncestors(target); }
    });
    if (typeof MutationObserver === 'function') {
      const observer = new MutationObserver(sync);
      for (const nodes of Object.values(sections)) for (const node of nodes) observer.observe(node,{attributes:true,attributeFilter:['hidden']});
    }
    win.addEventListener('matrix-profile-focus', () => {
      if (modes.includes('results') && available(find('profiles-panel'))) { activate('results',find('profiles-panel')); win.location.hash = 'profiles-panel'; }
    }, {capture:true});
  };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded',boot,{once:true}); else boot();
})();
