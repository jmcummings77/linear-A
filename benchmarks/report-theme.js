/** Apply the saved or system theme before the report's first paint. */
(() => {
  const key = "linear-a-report-theme";
  const preferences = new Set(["system", "light", "dark"]);
  let preference = "system";
  try {
    const saved = window.localStorage.getItem(key);
    if (preferences.has(saved)) preference = saved;
  } catch { /* Storage can be unavailable in private or local-file contexts. */ }

  const media = typeof window.matchMedia === "function"
    ? window.matchMedia("(prefers-color-scheme: dark)") : null;

  function apply() {
    const theme = preference === "system" ? (media?.matches ? "dark" : "light") : preference;
    if (document.documentElement.dataset.theme === theme) return;
    document.documentElement.dataset.theme = theme;
    window.dispatchEvent(new Event("matrix-theme-change"));
  }

  function systemChanged() {
    if (preference === "system") apply();
  }
  if (typeof media?.addEventListener === "function") media.addEventListener("change", systemChanged);
  else if (typeof media?.addListener === "function") media.addListener(systemChanged);
  apply();

  function bindSelect() {
    const select = document.getElementById("report-theme");
    if (!select) return;
    select.value = preference;
    select.addEventListener("change", () => {
      if (!preferences.has(select.value)) { select.value = preference; return; }
      preference = select.value;
      try {
        if (preference === "system") window.localStorage.removeItem(key);
        else window.localStorage.setItem(key, preference);
      } catch { /* The current report still changes theme when persistence fails. */ }
      apply();
    });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", bindSelect, { once: true });
  else bindSelect();
})();
