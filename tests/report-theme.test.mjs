import assert from "node:assert/strict";
import test from "node:test";
import { readFileSync } from "node:fs";
import vm from "node:vm";

const source = readFileSync(new URL("../benchmarks/report-theme.js", import.meta.url), "utf8");
const key = "linear-a-report-theme";

function report({ saved, dark = false, ready = "loading", storageFailure, legacyMedia = false,
  withoutMedia = false, withoutSelect = false, initialTheme } = {}) {
  const values = new Map(saved === undefined ? [] : [[key, saved]]);
  const writes = [], events = [], documentListeners = new Map(), mediaListeners = [], selectListeners = [];
  const select = {
    value: "uninitialized",
    addEventListener(type, handler) { assert.equal(type, "change"); selectListeners.push(handler); },
  };
  const document = {
    readyState: ready,
    documentElement: { dataset: initialTheme ? { theme: initialTheme } : {} },
    getElementById(id) { assert.equal(id, "report-theme"); return withoutSelect ? null : select; },
    addEventListener(type, handler, options) {
      assert.equal(type, "DOMContentLoaded");
      assert.equal(options.once, true);
      documentListeners.set(type, handler);
    },
  };
  const media = { matches: dark };
  if (legacyMedia) media.addListener = handler => mediaListeners.push(handler);
  else media.addEventListener = (type, handler) => {
    assert.equal(type, "change"); mediaListeners.push(handler);
  };
  const storage = {
    getItem(name) {
      assert.equal(name, key);
      if (storageFailure === "read") throw new Error("storage read unavailable");
      return values.get(name) ?? null;
    },
    setItem(name, value) {
      if (storageFailure === "write") throw new Error("storage write unavailable");
      writes.push(["set", name, value]); values.set(name, value);
    },
    removeItem(name) {
      if (storageFailure === "write") throw new Error("storage removal unavailable");
      writes.push(["remove", name]); values.delete(name);
    },
  };
  const window = {
    dispatchEvent(event) {
      events.push({ type: event.type, theme: document.documentElement.dataset.theme });
      return true;
    },
  };
  Object.defineProperty(window, "localStorage", { get() {
    if (storageFailure === "access") throw new Error("storage access unavailable");
    return storage;
  } });
  if (!withoutMedia) window.matchMedia = query => {
    assert.equal(query, "(prefers-color-scheme: dark)"); return media;
  };
  class Event { constructor(type) { this.type = type; } }
  vm.runInNewContext(source, { window, document, Event });
  return {
    document, select, values, writes, events,
    get theme() { return document.documentElement.dataset.theme; },
    ready() {
      document.readyState = "interactive";
      const callback = documentListeners.get("DOMContentLoaded");
      documentListeners.delete("DOMContentLoaded");
      callback?.();
    },
    system(dark) {
      media.matches = dark;
      for (const listener of mediaListeners) listener({ matches: dark });
    },
    choose(value) {
      select.value = value;
      for (const listener of selectListeners) listener();
    },
  };
}

test("system preference applies before DOM readiness and initializes the select later", () => {
  for (const dark of [false, true]) {
    const page = report({ dark });
    assert.equal(page.theme, dark ? "dark" : "light");
    assert.equal(page.select.value, "uninitialized");
    assert.deepEqual(page.events, [{ type: "matrix-theme-change", theme: page.theme }]);
    assert.deepEqual(page.writes, []);
    page.ready();
    assert.equal(page.select.value, "system");
    assert.equal(page.events.length, 1);
  }
});

test("saved manual preferences override the system and persist in the select", () => {
  for (const [saved, dark] of [["light", true], ["dark", false]]) {
    const page = report({ saved, dark });
    assert.equal(page.theme, saved);
    page.ready();
    assert.equal(page.select.value, saved);
    page.system(!dark);
    page.system(dark);
    assert.equal(page.theme, saved);
    assert.equal(page.events.length, 1, "OS changes must not redraw a manual theme");
    assert.deepEqual(page.writes, []);
  }
});

test("system changes update immediately and emit only effective transitions", () => {
  const page = report();
  page.system(true);
  assert.equal(page.theme, "dark", "an OS change before DOMContentLoaded still applies");
  page.system(true);
  page.ready();
  page.system(false);
  assert.equal(page.select.value, "system");
  assert.deepEqual(page.events.map(event => event.theme), ["light", "dark", "light"]);
  assert.deepEqual(page.writes, []);
});

test("manual changes persist and selecting system resumes the current OS preference", () => {
  const page = report({ dark: true, ready: "complete" });
  assert.equal(page.select.value, "system");
  page.choose("light");
  assert.equal(page.theme, "light");
  assert.equal(page.values.get(key), "light");
  page.system(false);
  page.system(true);
  assert.equal(page.theme, "light");
  page.choose("dark");
  assert.equal(page.values.get(key), "dark");
  assert.equal(page.theme, "dark");
  page.choose("system");
  assert.equal(page.values.has(key), false);
  assert.equal(page.theme, "dark");
  page.system(false);
  assert.equal(page.theme, "light");
  assert.equal(page.select.value, "system");
  assert.deepEqual(page.events.map(event => event.theme), ["dark", "light", "dark", "light"]);
  assert.deepEqual(page.writes, [["set", key, "light"], ["set", key, "dark"], ["remove", key]]);
});

test("selecting an override with the same effective theme saves without redrawing", () => {
  const page = report({ dark: true, ready: "interactive", initialTheme: "dark" });
  assert.equal(page.events.length, 0);
  page.choose("dark");
  page.choose("dark");
  assert.equal(page.events.length, 0);
  assert.equal(page.values.get(key), "dark");
  page.choose("system");
  assert.equal(page.events.length, 0);
  page.system(false);
  assert.deepEqual(page.events, [{ type: "matrix-theme-change", theme: "light" }]);
});

test("unavailable storage never prevents current-page theme changes", () => {
  for (const storageFailure of ["read", "write", "access"]) {
    const page = report({ storageFailure, ready: "complete" });
    assert.equal(page.theme, "light");
    page.choose("dark");
    assert.equal(page.theme, "dark");
    assert.equal(page.select.value, "dark");
    page.system(true);
    page.system(false);
    assert.equal(page.theme, "dark");
    page.choose("system");
    assert.equal(page.theme, "light");
    page.system(true);
    assert.equal(page.theme, "dark");
  }
});

test("invalid saved or selected preferences cannot become theme attributes", () => {
  for (const saved of ["", "unknown", "DARK", "<script>", "system"]) {
    const page = report({ saved, dark: true, ready: "complete" });
    assert.equal(page.theme, "dark");
    assert.equal(page.select.value, "system");
    page.choose("unknown");
    assert.equal(page.select.value, "system");
    assert.equal(page.theme, "dark");
    assert.deepEqual(page.writes, []);
  }
});

test("missing controls or media APIs retain a usable theme, including legacy media listeners", () => {
  const absentSelect = report({ withoutSelect: true });
  absentSelect.ready();
  absentSelect.system(true);
  assert.equal(absentSelect.theme, "dark");
  const absentMedia = report({ withoutMedia: true, ready: "complete" });
  assert.equal(absentMedia.theme, "light");
  absentMedia.choose("dark");
  assert.equal(absentMedia.theme, "dark");
  const legacy = report({ legacyMedia: true });
  legacy.system(true);
  assert.equal(legacy.theme, "dark");
});
