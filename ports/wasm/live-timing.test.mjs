import assert from "node:assert/strict";
import test from "node:test";
import { runBenchmark } from "../../benchmarks/live-worker.mjs";

const MAX_ITERATIONS = 1_048_576;
const config = (samples = 3) => ({ operation: "trace", size: 16, samples, seed: 17 });

// Advance simulated time only when the workload runs. Quantization therefore
// models a real coarse clock, rather than supplying a canned list of readings.
async function withClock(options, body) {
  const performanceDescriptor = Object.getOwnPropertyDescriptor(globalThis, "performance");
  const dateDescriptor = Object.getOwnPropertyDescriptor(Date, "now");
  const state = { elapsed: 0, wall: 0, cost: .01, calls: 0, created: 0, disposed: 0,
    live: 0, peakLive: 0, reads: [], progress: [], ...options };
  const now = () => {
    const value = state.readPerformance ? state.readPerformance(state)
      : state.frozen ? 0 : state.quantum ? Math.floor(state.elapsed / state.quantum) * state.quantum : state.elapsed;
    state.reads.push({ calls: state.calls, value });
    return value;
  };
  const advance = () => {
    state.calls++;
    state.elapsed += state.cost;
    state.onConsume?.(state);
  };
  class Matrix {
    constructor(rows, cols, values) {
      this.rows = rows; this.cols = cols; this.values = values;
      this.disposed = false;
      this.diagonal = 0;
      for (let i = 0; i < Math.min(rows, cols); i++) this.diagonal += values[i * cols + i];
      state.created++; state.live++;
      state.peakLive = Math.max(state.peakLive, state.live);
    }
    trace() { assert.equal(this.disposed, false); advance(); return this.diagonal; }
    toArray() { assert.equal(this.disposed, false); return this.values; }
    checksum() {
      assert.equal(this.disposed, false);
      return this.values[0] + this.values[Math.floor(this.values.length / 2)] + this.values.at(-1);
    }
    dispose() {
      assert.equal(this.disposed, false, "each owned matrix must be disposed exactly once");
      this.disposed = true; state.disposed++; state.live--;
    }
    static rotationAxisAngle(axis, angle) {
      assert.equal(axis.disposed, false);
      advance();
      const norm = Math.hypot(...axis.values), [x, y, z] = axis.values.map(value => value / norm);
      const c = Math.cos(angle), s = Math.sin(angle), t = 2 * Math.sin(angle / 2) ** 2;
      return new Matrix(3, 3, [c+x*x*t, x*y*t-z*s, x*z*t+y*s,
        y*x*t+z*s, c+y*y*t, y*z*t-x*s, z*x*t-y*s, z*y*t+x*s, c+z*z*t]);
    }
  }
  Object.defineProperty(globalThis, "performance", { configurable: true, value: { now } });
  Object.defineProperty(Date, "now", { configurable: true, value: () => state.readWall ? state.readWall(state) : state.wall });
  try {
    return await body({ state, Matrix, run: (settings = config(), onProgress = () => {}) =>
      runBenchmark(Matrix, settings, message => {
        state.progress.push({ ...message });
        onProgress(message, state);
      }) });
  } finally {
    Object.defineProperty(globalThis, "performance", performanceDescriptor);
    Object.defineProperty(Date, "now", dateDescriptor);
    assert.equal(state.live, 0, "successful, discarded, and failed batches must release every matrix");
    assert.equal(state.disposed, state.created);
  }
}

function assertSamples(result, samples) {
  assert.equal(result.samples.length, samples);
  assert.equal(result.checksum_verified, true);
  assert.ok(Number.isSafeInteger(result.iterations) && result.iterations > 0 && result.iterations <= MAX_ITERATIONS);
  // The fixed 16×16, seed-17 trace workload sums to −55/16. The 3D rotation
  // checksum samples its three diagonal entries, whose sum is 1 + 2 cos θ.
  const expected = result.operation === "rotation3d" ? 1 + 2 * Math.cos(.5) : -55 / 16;
  for (const sample of result.samples) {
    assert.ok(Number.isFinite(sample.elapsed_ms) && sample.elapsed_ms >= 15, "every retained sample must be measurable");
    assert.equal(sample.ns_per_op, sample.elapsed_ms * 1e6 / result.iterations, "one common iteration count must explain every timing");
    assert.ok(Math.abs(sample.checksum / result.iterations - expected) <= 1e-8 * Math.max(1, Math.abs(expected)),
      "one common iteration count must explain every checksum");
  }
  for (const key of ["median_ns", "mad_ns", "min_ns", "max_ns"]) assert.ok(Number.isFinite(result[key]), `${key} must be finite`);
  assert.ok(result.min_ns > 0 && result.mad_ns >= 0 && result.min_ns <= result.median_ns && result.median_ns <= result.max_ns);
}

// These tests replace global clocks and must remain serial within this process.
const serial = (name, body) => test(name, { concurrency: false }, body);

serial("coarse clocks calibrate tiny workloads beyond the former 4096-operation limit", async () => {
  await withClock({ cost: .001, quantum: 2 }, async ({ run, state }) => {
    const result = await run();
    assertSamples(result, 3);
    assert.ok(result.iterations > 4096);
    assert.ok(state.reads.some((read, i) => i > 0 && read.calls > state.reads[i-1].calls && read.value === state.reads[i-1].value),
      "the simulated clock must actually produce zero-duration work");
    assert.equal(state.created, 2);
  });
});

serial("a speedup that makes the second calibration read zero grows another batch", async () => {
  await withClock({ cost: 1, quantum: 1, onConsume(state) { if (state.calls === 6) state.cost = .0005; } }, async ({ run, state }) => {
    const result = await run();
    assertSamples(result, 3);
    assert.ok(result.iterations > 4096);
    assert.ok(state.reads.some((read, i) => i > 0 && read.calls > 6 && read.calls < 100
      && read.calls > state.reads[i-1].calls && read.value === state.reads[i-1].value),
    "a later calibration batch must have run inside a single clock tick");
  });
});

serial("rotation3d remains usable across fresh 3, 5, and 3 sample runs after a late speedup", async () => {
  const runs = [];
  for (const samples of [3, 5, 3]) {
    await withClock({ cost: .02, quantum: 1 }, async ({ run, state }) => {
      let spedUp = false;
      const result = await run({ operation: "rotation3d", size: 3, samples, seed: 17 }, message => {
        if (samples === 5 && message.completed === 3 && !spedUp) {
          spedUp = true; state.cost = .0002;
          // Keep time monotonic and put the next small batch wholly inside a tick.
          state.elapsed = Math.ceil(state.elapsed) + .1;
        }
      });
      assertSamples(result, samples);
      if (samples === 5) {
        assert.equal(spedUp, true);
        assert.ok(state.progress.filter(message => message.completed === 1).length >= 2,
          "earlier samples must be discarded and the full sample set restarted");
        assert.equal(state.progress.at(-1).completed, 5);
      }
      assert.ok(state.peakLive <= 2, "rotation results must not accumulate across retries");
      runs.push(result);
    });
  }
  assert.ok(runs[1].iterations > 10 * runs[0].iterations);
  assert.ok(runs[2].iterations < runs[1].iterations, "a fresh run must calibrate independently");
  assert.ok(runs[1].median_ns < runs[0].median_ns / 10);
});

serial("a frozen performance clock fails at a finite work limit without fabricated timings", async () => {
  await withClock({ frozen: true }, async ({ run, state }) => {
    await assert.rejects(run(), /timer|clock|measur|resolution/i);
    assert.ok(state.calls >= MAX_ITERATIONS, "the final allowed batch must be attempted");
    assert.ok(state.calls < 4 * MAX_ITERATIONS, "a stopped clock must not create an endless retry loop");
    assert.ok(!state.progress.some(message => message.completed > 0));
    assert.equal(state.created, 2);
  });
});

serial("the maximum batch can retain measurable samples below the calibration target without inventing time", async () => {
  await withClock({ cost: .00002, quantum: 1 }, async ({ run }) => {
    const result = await run();
    assertSamples(result, 3);
    assert.equal(result.iterations, MAX_ITERATIONS);
    assert.ok(result.samples.every(sample => sample.elapsed_ms < 30));
  });
});

serial("an independent wall-clock budget interrupts a large batch even when the performance clock is frozen", async () => {
  await withClock({ frozen: true, onConsume(state) { if (state.calls === 1500) state.wall = 15_001; } }, async ({ run, state }) => {
    await assert.rejects(run(), /time budget/i);
    assert.ok(state.calls >= 1500 && state.calls <= 1500 + 1024,
      "the worker must notice the budget within one checkpoint interval, not finish the oversized batch");
    assert.equal(state.created, 2);
  });
});

serial("warmup time is included in the independent safety budget", async () => {
  await withClock({ frozen: true, onConsume(state) { state.wall = 15_001; } }, async ({ run, state }) => {
    await assert.rejects(run(), /time budget/i);
    assert.equal(state.calls, 5, "the budget must be checked after warmup, before the timed operation");
    assert.equal(state.created, 2);
  });
});

serial("invalid clock readings reject the run and clean up owned inputs", async () => {
  for (const invalid of [NaN, Infinity, -Infinity]) {
    for (const clock of ["performance", "wall"]) {
      await withClock({
        readPerformance: state => clock === "performance" && state.calls >= 20 ? invalid : state.elapsed,
        readWall: state => clock === "wall" && state.calls >= 20 ? invalid : 0,
      }, async ({ run, state }) => {
        await assert.rejects(run(), /clock|elapsed|timer/i, `${clock} must reject ${String(invalid)}`);
        assert.equal(state.created, 2);
        assert.ok(state.calls < 4096, "invalid clocks must fail promptly");
      });
    }
  }
});

serial("invalid initial clocks reject before allocating inputs or creating an unusable deadline", async () => {
  for (const invalid of [NaN, Infinity, -Infinity]) {
    for (const clock of ["performance", "wall"]) {
      await withClock({
        readPerformance: () => clock === "performance" ? invalid : 0,
        readWall: () => clock === "wall" ? invalid : 0,
      }, async ({ run, state }) => {
        await assert.rejects(run(), /clock|elapsed|timer/i);
        assert.equal(state.created, 0);
        assert.equal(state.calls, 0);
      });
    }
  }
});

serial("a clock that goes backwards during calibration rejects negative elapsed time", async () => {
  await withClock({ readPerformance: state => state.calls >= 20 ? -1 : state.elapsed }, async ({ run, state }) => {
    await assert.rejects(run(), /clock|elapsed|timer/i);
    assert.equal(state.created, 2);
    assert.ok(state.calls < 4096);
  });
});
