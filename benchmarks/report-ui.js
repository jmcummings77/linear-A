/** Decorative action icons; native controls and visible labels remain authoritative. */
(() => {
  const actionIcons = {
    csv: "download", download: "download", "live-export": "download",
    previous: "previous", "geometry-prev": "previous",
    next: "next", "geometry-next": "next",
    reset: "reset", "geometry-reset": "reset", "pca-reset": "reset", "fit-reset": "reset",
    "live-check": "check", "live-run": "play", "live-stop": "stop",
    swap: "swap", "profile-more": "more",
    "timing-language-all": "check", "timing-language-none": "clear",
    "geometry-settings-open": "settings", "geometry-settings-close": "check",
    "geometry-randomize": "shuffle", "geometry-profile-link": "search",
    "pca-add": "add", "fit-add": "add", "pca-remove": "remove", "fit-remove": "remove",
  };
  const playbackIds = new Set(["play", "geometry-play"]);

  function decorate() {
    const playback = [];
    for (const control of document.querySelectorAll("button[id], a[download], [data-report-icon]")) {
      // An explicit icon is an opt-in for new controls; "none" opts out.
      const explicit = control.getAttribute("data-report-icon");
      if (explicit === "none") continue;
      const isPlayback = playbackIds.has(control.id) ||
        (control.tagName === "BUTTON" && (explicit === "play" || explicit === "pause"));
      const icon = explicit || (isPlayback ? "play" : actionIcons[control.id]) ||
        (control.tagName === "A" && control.hasAttribute("download") ? "download" : null);
      if (!icon) continue;
      // The existing settings button has only an SVG. Give it a visible label too.
      if (control.id === "geometry-settings-open" && !control.textContent.trim()) {
        control.textContent = "Settings";
        control.removeAttribute("title");
      }
      control.setAttribute("data-report-icon", icon);
      if (isPlayback) playback.push(control);
    }

    function sync(control) {
      const icon = /^\s*pause\b/i.test(control.textContent) ? "pause" : "play";
      if (control.getAttribute("data-report-icon") !== icon) control.setAttribute("data-report-icon", icon);
    }
    playback.forEach(sync);
    // Observe just the playback labels, not plots or continuously updated frame data.
    // The CSS pseudo-element survives all existing textContent replacements.
    if (typeof MutationObserver === "function" && playback.length) {
      const observer = new MutationObserver(records => {
        for (const control of playback) {
          if (records.some(record => record.target === control || control.contains(record.target))) sync(control);
        }
      });
      for (const control of playback) observer.observe(control, { childList: true, characterData: true, subtree: true });
      window.addEventListener("pagehide", event => { if (!event.persisted) observer.disconnect(); });
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", decorate, { once: true });
  else decorate();
})();
