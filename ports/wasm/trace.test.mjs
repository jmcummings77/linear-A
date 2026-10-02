import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";
import { test } from "node:test";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";

const traceURL = process.env.LINEAR_A_WASM_TRACE_MODULE
  ? pathToFileURL(resolve(process.env.LINEAR_A_WASM_TRACE_MODULE))
  : new URL("../../.build/wasm-trace/matrix.mjs", import.meta.url);
const hasTrace = existsSync(traceURL);
const options = { moduleUrl: traceURL };
const Matrix = hasTrace ? (await createMatrixAPI(options)).Matrix : undefined;
const traceTest = (name, run) => test(name, { skip: hasTrace ? false : "build ports/wasm/build.py --trace first" }, run);
const sourceSnapshot = () => JSON.parse(readFileSync(new URL("trace-source.json", traceURL), "utf8"));

function use(matrices, run) {
  try { return run(...matrices); } finally { for (const matrix of matrices) matrix.dispose(); }
}

traceTest("trace captures the actual 27 multiply updates, including signed fractions and source lines", () => {
  const left = [1, -2, 0.5, 0, 3, -1, 2, 0.25, -4];
  const right = [0.5, 2, -1, 3, -0.5, 0, -2, 1, 4];
  const expected = [-6.5, 3.5, 1, 11, -2.5, -4, 9.75, -0.125, -18];
  use([new Matrix(3, 3, left), new Matrix(3, 3, right)], (a, b) => {
    const captured = a.multiplyWithTrace(b);
    use([captured.matrix, a.multiply(b)], (result, plain) => {
      assert.equal(captured.source, "ports/c/matrix.c");
      assert.deepEqual([...result.toArray()], expected);
      assert.deepEqual(result.toArray(), plain.toArray());
      assert.equal(captured.steps.length, 27);
      const partial = Array(9).fill(0);
      const source = sourceSnapshot();
      captured.steps.forEach((step, index) => {
        const row = Math.floor(index / 9), col = Math.floor(index / 3) % 3, k = index % 3;
        assert.deepEqual([step.row, step.col, step.k], [row, col, k]);
        assert.equal(step.left, left[row * 3 + k]);
        assert.equal(step.right, right[k * 3 + col]);
        assert.equal(step.product, step.left * step.right);
        partial[row * 3 + col] += step.product;
        assert.equal(step.sum, partial[row * 3 + col]);
        assert.ok(Number.isSafeInteger(step.source_line) && step.source_line > 0);
        assert.match(source.lines[step.source_line - source.start_line], /m_trace_record\(.*__LINE__\)/);
      });
      assert.deepEqual(partial, expected);
    });
  });
});

traceTest("trace source manifest contains the compiled function and full-source hash without local paths", () => {
  const snapshot = sourceSnapshot();
  const source = readFileSync(new URL("../c/matrix.c", import.meta.url));
  const lines = source.toString("utf8").split(/\r?\n/);
  assert.deepEqual(Object.keys(snapshot).sort(), ["lines", "path", "source_sha256", "start_line"]);
  assert.equal(snapshot.path, "ports/c/matrix.c");
  assert.equal(snapshot.source_sha256, createHash("sha256").update(source).digest("hex"));
  assert.ok(Number.isSafeInteger(snapshot.start_line) && snapshot.start_line > 0);
  assert.deepEqual(snapshot.lines, lines.slice(snapshot.start_line - 1, snapshot.start_line - 1 + snapshot.lines.length));
  assert.match(snapshot.lines[0], /^matrix_status m_multiply\(/);
  assert.equal(snapshot.lines.at(-1), "}");
  assert.ok(snapshot.lines.some(line => line.includes("m_trace_record(")));
  assert.doesNotMatch(JSON.stringify(snapshot), /\/Users\/|\/private\/|\/tmp\/|[A-Z]:\\\\/);
});

traceTest("trace preserves left-to-right floating-point rounding and rectangular multiplication", () => {
  use([new Matrix(1, 3, [1e16, 1, -1e16]), new Matrix(3, 2, [1, 0.5, 1, 0.5, 1, 0.5])], (a, b) => {
    const { matrix, steps } = a.multiplyWithTrace(b);
    use([matrix], result => {
      assert.deepEqual([result.rows, result.cols], [1, 2]);
      assert.deepEqual([...result.toArray()], [0, 0]);
      assert.deepEqual(steps.map(step => step.sum), [1e16, 1e16, 0, 5e15, 5e15, 0]);
    });
  });
});

traceTest("trace snapshots survive subsequent calls, input edits, and memory growth", () => {
  use([new Matrix(1, 1, [3]), new Matrix(1, 1, [-2])], (a, b) => {
    const first = a.multiplyWithTrace(b);
    use([first.matrix], result => {
      const saved = structuredClone(first.steps);
      a.set(0, 0, 9);
      const second = a.multiplyWithTrace(b);
      use([second.matrix, new Matrix(2048, 2048)], () => {
        assert.deepEqual(first.steps, saved);
        assert.equal(result.get(0, 0), -6);
        assert.equal(second.steps[0].sum, -18);
      });
    });
  });
});

traceTest("trace enforces its 8-by-8 bound and accepts exactly 512 terms", () => {
  use([Matrix.identity(8), Matrix.identity(9), new Matrix(2, 3), new Matrix(0, 9)], (eight, nine, rectangle, emptyWide) => {
    const captured = eight.multiplyWithTrace(eight);
    use([captured.matrix], result => {
      assert.equal(captured.steps.length, 512);
      assert.deepEqual(result.toArray(), eight.toArray());
    });
    assert.throws(() => nine.multiplyWithTrace(nine), /0 to 8/);
    assert.throws(() => rectangle.multiplyWithTrace(rectangle), /dimensions/);
    assert.throws(() => emptyWide.multiplyWithTrace(nine), /0 to 8/);
  });
});

traceTest("zero inner dimensions produce a zero result and no invented trace events", () => {
  use([new Matrix(2, 0), new Matrix(0, 3)], (a, b) => {
    const { matrix, steps } = a.multiplyWithTrace(b);
    use([matrix], result => {
      assert.deepEqual([result.rows, result.cols], [2, 3]);
      assert.deepEqual([...result.toArray()], [0, 0, 0, 0, 0, 0]);
      assert.deepEqual(steps, []);
    });
  });
});

traceTest("failed arithmetic cannot publish partial results and later captures still work", () => {
  use([new Matrix(1, 1, [1e308]), new Matrix(1, 1, [2]), new Matrix(1, 2, [1e308, 1e308]), new Matrix(2, 1, [1, 1])], (large, two, row, col) => {
    assert.throws(() => large.multiplyWithTrace(two), /nonfinite/);
    assert.throws(() => row.multiplyWithTrace(col), /nonfinite/);
    assert.equal(large.get(0, 0), 1e308);
    const { matrix, steps } = two.multiplyWithTrace(two);
    use([matrix], result => {
      assert.equal(result.get(0, 0), 4);
      assert.equal(steps.length, 1);
      assert.equal(steps[0].sum, 4);
    });
  });
});

traceTest("tracing rejects disposed and foreign-instance handles", async () => {
  const { Matrix: OtherMatrix } = await createMatrixAPI(options);
  use([new Matrix(1, 1), new OtherMatrix(1, 1), new Matrix(1, 1)], (a, foreign, disposed) => {
    disposed.dispose();
    assert.throws(() => a.multiplyWithTrace(foreign), /same WebAssembly instance/);
    assert.throws(() => a.multiplyWithTrace(disposed), /disposed/);
    assert.throws(() => disposed.multiplyWithTrace(a), /disposed/);
  });
});

test("ordinary WebAssembly builds do not expose trace instrumentation", async () => {
  const normalURL = new URL("../../.build/wasm/matrix.mjs", import.meta.url);
  if (!existsSync(normalURL)) return;
  const { Matrix: NormalMatrix } = await createMatrixAPI({ moduleUrl: normalURL });
  assert.equal("multiplyWithTrace" in NormalMatrix.prototype, false);
});
