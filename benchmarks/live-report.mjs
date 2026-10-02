// Browser controls for the embedded WebAssembly test worker.
(() => {
  const get = id => document.getElementById(id);
  const bundle = JSON.parse(get("live-data").textContent);
  const controls = ["live-check", "live-run", "live-operation", "live-size", "live-samples"].map(get);
  const history = [];
  let active = null, workerUrl = null, timeout = null, lastChecks = null;
  const format = ns => ns >= 1e6 ? `${(ns / 1e6).toFixed(3)} ms` : `${(ns / 1e3).toFixed(2)} µs`;
  const dimensions = (operation, size) => operation === "cross" ? "3D vectors" : `${size} × ${size}`;
  const text = (tag, value, className) => {
    const node = document.createElement(tag);
    node.textContent = value;
    if (className) node.className = className;
    return node;
  };

  function publicBrowser() {
    // Export only recognized product versions, never the full user-agent string.
    for (const [name, pattern] of [["Edge", /Edg\/([\d.]+)/], ["Firefox", /Firefox\/([\d.]+)/],
      ["Chrome", /Chrome\/([\d.]+)/], ["Safari", /Version\/([\d.]+).*Safari\//]]) {
      const match = navigator.userAgent.match(pattern);
      if (match) return { name, version: match[1] };
    }
    return { name: "Unknown" };
  }

  function publicChecks(result) {
    if (!result) return null;
    return { completed_at: result.completed_at, passed: result.passed, total: result.total,
      checks: result.checks.map((check, index) => ({
        name: bundle.fixtures[index]?.name || `Check ${index + 1}`, passed: check.passed === true,
        ...(!check.passed ? { error: "Correctness check failed." } : {}),
      })) };
  }

  function publicBenchmark(result) {
    const keys = ["size", "seed", "iterations", "median_ns", "mad_ns", "min_ns", "max_ns"];
    const operation = ["add", "subtract", "scale", "transpose", "multiply", "cross", "rotation2d", "rotation3d", "trace", "eigen_symmetric", "eigen_general", "determinant", "determinant_lu", "determinant_spd_lu", "determinant_cholesky"]
      .includes(result.operation) ? result.operation : "unknown";
    return { operation, completed_at: result.completed_at,
      ...Object.fromEntries(keys.map(key => [key, result[key]])),
      checksum_verified: result.checksum_verified === true,
      samples: result.samples.map(sample => ({ elapsed_ms: sample.elapsed_ms,
        ns_per_op: sample.ns_per_op, checksum: sample.checksum })) };
  }

  function status(message, state = "") {
    get("live-status").textContent = message;
    get("live-status").className = `live-status ${state}`;
  }

  function finish() {
    if (active) active.terminate();
    active = null;
    if (workerUrl) URL.revokeObjectURL(workerUrl);
    workerUrl = null;
    clearTimeout(timeout);
    timeout = null;
    controls.forEach(control => { control.disabled = false; });
    get("live-stop").disabled = true;
    get("live-panel").setAttribute("aria-busy", "false");
  }

  function showCheck(check) {
    const row = document.createElement("tr");
    row.append(text("td", check.name), text("td", check.passed ? "PASS" : "FAIL", check.passed ? "pass" : "failed"),
      text("td", check.error || ""));
    get("live-check-rows").append(row);
  }

  function showChecks(result) {
    lastChecks = { ...result, completed_at: new Date().toISOString() };
    get("live-check-rows").replaceChildren();
    for (const check of result.checks) showCheck(check);
    get("live-check-summary").textContent = `${result.passed}/${result.total} correctness checks passed`;
    get("live-export").disabled = false;
    if (result.passed !== result.total) get("live-check-details").open = true;
  }

  function showBenchmark(result) {
    if (!result.checksum_verified) throw new Error("The live result did not pass checksum verification.");
    history.unshift({ ...result, completed_at: new Date().toISOString() });
    // Keep an on-page session small even if the report is left open for hours.
    history.splice(20);
    get("live-results").hidden = false;
    get("live-results-empty").hidden = true;
    get("live-result-rows").replaceChildren();
    for (const run of history) {
      const row = document.createElement("tr");
      for (const value of [run.operation, dimensions(run.operation, run.size), format(run.median_ns),
        format(run.mad_ns), `${format(run.min_ns)} – ${format(run.max_ns)}`,
        `${run.samples.length} × ${run.iterations}`, "Verified"]) row.append(text("td", value));
      get("live-result-rows").append(row);
    }
    get("live-export").disabled = false;
  }

  function sizes() {
    const old = get("live-size").value;
    get("live-size").replaceChildren();
    const operation = get("live-operation").value;
    const allowed = operation === "rotation2d" ? [2] : operation === "cross" || operation === "rotation3d" ? [3]
      : operation.startsWith("determinant") ? [4, 8, 16, 24, 48]
      : operation.startsWith("eigen_") ? [4, 8, 16, 32, 48] : [16, 48, 64, 128, 256];
    for (const size of allowed) {
      const option = text("option", dimensions(operation, size));
      option.value = String(size);
      get("live-size").append(option);
    }
    get("live-size").value = allowed.includes(Number(old)) ? old : String(allowed.includes(48) ? 48 : allowed[0]);
    get("live-operation-note").textContent = operation === "cross" ? "Cross two 3D column vectors; the result is a vector, not a matrix product."
      : operation === "rotation2d" ? "Construct a 2D rotation matrix at 0.5 radians."
      : operation === "rotation3d" ? "Construct a 3D rotation matrix at 0.5 radians about the axis (1, 2, 3)." : "";
  }

  function start(type) {
    if (active) return;
    let worker;
    try {
      workerUrl = URL.createObjectURL(new Blob([bundle.worker_source], { type: "text/javascript" }));
      worker = new Worker(workerUrl, { type: "module", name: "matrix-live-tests" });
      active = worker;
      controls.forEach(control => { control.disabled = true; });
      get("live-stop").disabled = false;
      get("live-panel").setAttribute("aria-busy", "true");
      get("live-check-rows").replaceChildren();
      get("live-check-summary").textContent = "Checking the embedded implementation…";
      status("Loading WebAssembly in this browser…");
      worker.onmessage = ({ data: message }) => {
        if (active !== worker) return;
        try {
          if (message.type === "progress") {
            status(message.message || "Running…");
            if (message.check) showCheck(message.check);
          } else if (message.type === "done") {
            if (message.kind === "checks") {
              showChecks(message.result);
              const passed = message.result.passed === message.result.total;
              status(passed ? `All ${message.result.total} correctness checks passed in this browser.` : "Some correctness checks failed. See the case results below.", passed ? "pass" : "failed");
            } else {
              if (message.checks) showChecks(message.checks);
              showBenchmark(message.result);
              status(`${message.result.operation} ${dimensions(message.result.operation, message.result.size)}: ${format(message.result.median_ns)} median per operation. Checksums verified.`, "pass");
            }
            finish();
          } else if (message.type === "error") {
            if (message.checks) showChecks(message.checks);
            else get("live-check-summary").textContent = "Run incomplete — see status above";
            status(`Live run failed: ${message.message}`, "failed");
            finish();
          }
        } catch (error) {
          status(`Live run failed: ${error.message}`, "failed");
          finish();
        }
      };
      worker.onerror = event => {
        if (active !== worker) return;
        event.preventDefault();
        status(`Live tests could not run: ${event.message || "This browser or its content policy blocked the worker."}`, "failed");
        get("live-check-summary").textContent = "Run incomplete";
        finish();
      };
      timeout = setTimeout(() => {
        if (active !== worker) return;
        finish();
        status("Stopped after 45 seconds. Try a smaller matrix or fewer samples.", "failed");
        get("live-check-summary").textContent = "Run stopped before completion";
      }, 45000);
      worker.postMessage({ type, bundle, config: { operation: get("live-operation").value,
        size: Number(get("live-size").value), samples: Number(get("live-samples").value), seed: 17 } });
    } catch (error) {
      finish();
      status(`Live tests could not start: ${error.message}`, "failed");
    }
  }

  get("live-operation").onchange = sizes;
  sizes();
  get("live-check").onclick = () => start("checks");
  get("live-run").onclick = () => start("benchmark");
  get("live-stop").onclick = () => {
    finish();
    status("Stopped. You can start another live run.");
    get("live-check-summary").textContent = "Run stopped before completion";
  };
  get("live-export").onclick = () => {
    const payload = { schema_version: 1, kind: "browser-webassembly-session", exported_at: new Date().toISOString(),
      bundle_sha256: bundle.sha256, browser: publicBrowser(), checks: publicChecks(lastChecks),
      benchmarks: history.map(publicBenchmark),
      timing: "Browser worker performance.now(); allocation, JS/WASM calls, checksum, and disposal included; module loading, inputs, checks, and warmup excluded." };
    const serialized = JSON.stringify(payload, null, 2);
    get("live-export-json").textContent = serialized;
    get("live-export-details").hidden = false;
    const url = URL.createObjectURL(new Blob([serialized], { type: "application/json" }));
    const link = document.createElement("a");
    link.href = url;
    link.download = "browser-wasm-results.json";
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  if (!bundle.available || typeof Worker === "undefined" || typeof WebAssembly === "undefined") {
    controls.forEach(control => { control.disabled = true; });
    status(bundle.reason || "This browser does not support the WebAssembly and worker features needed for live tests.", "unavailable");
    get("live-check-summary").textContent = "Live tests unavailable";
  } else {
    status(`Ready to run ${bundle.fixtures.length} correctness checks or a browser benchmark.`);
    get("live-bundle").textContent = `Embedded WASM: ${(bundle.byte_length / 1024).toFixed(1)} KiB · bundle ${bundle.sha256.slice(0, 16)}`;
  }
})();
