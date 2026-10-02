import assert from "node:assert/strict";
import { test } from "node:test";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";

const options = process.env.LINEAR_A_WASM_MODULE ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);

function use(matrices, operation) {
  try { operation(...matrices); } finally { for (const matrix of matrices) matrix.dispose(); }
}

test("zero, identity, and copied values own independent WASM storage", () => {
  const values = [1, 2, 3, 4];
  use([new Matrix(2, 2, values), new Matrix(2, 3), Matrix.identity(3)], (a, zero, identity) => {
    values[0] = 99;
    assert.equal(a.get(0, 0), 1);
    assert.deepEqual([...zero.toArray()], [0, 0, 0, 0, 0, 0]);
    assert.equal(identity.trace(), 3);
    assert.equal(identity.determinant(), 1);
    use([a.copy()], copy => {
      copy.set(0, 0, -5);
      assert.equal(a.get(0, 0), 1);
    });
    const exported = a.toArray();
    exported[0] = 88;
    assert.equal(a.get(0, 0), 1);
  });
});

test("rectangular arithmetic, transpose, and extracted copies", () => {
  use([new Matrix(2, 3, [1, 2, 3, 4, 5, 6]), new Matrix(3, 2, [7, 8, 9, 10, 11, 12])], (a, b) => {
    use([a.multiply(b), a.transpose(), a.add(a), a.subtract(a), a.scale(0.5)], (product, transpose, sum, difference, scaled) => {
      assert.deepEqual([...product.toArray()], [58, 64, 139, 154]);
      assert.deepEqual([transpose.rows, transpose.cols], [3, 2]);
      assert.deepEqual([...transpose.toArray()], [1, 4, 2, 5, 3, 6]);
      assert.deepEqual([...sum.toArray()], [2, 4, 6, 8, 10, 12]);
      assert.deepEqual([...difference.toArray()], [0, 0, 0, 0, 0, 0]);
      assert.deepEqual([...scaled.toArray()], [0.5, 1, 1.5, 2, 2.5, 3]);
      assert.equal(a.checksum(), 11);
    });
    assert.deepEqual([...a.row(1)], [4, 5, 6]);
    assert.deepEqual([...a.column(1)], [2, 5]);
    const row = a.row(0); row[0] = -100;
    assert.equal(a.get(0, 0), 1);
    assert.throws(() => a.add(b), /dimensions/);
    assert.throws(() => a.multiply(a), /dimensions/);
    assert.throws(() => a.trace(), /dimensions/);
    assert.throws(() => a.determinant(), /dimensions/);
    assert.deepEqual(a.triangular(), { upper: false, lower: false });
  });
});

test("determinant base cases, row pivoting, singularity, and triangular predicates", () => {
  use([
    new Matrix(), new Matrix(1, 1, [-7]), new Matrix(2, 2, [0, 1, 1, 0]),
    new Matrix(2, 2, [1, 2, 2, 4]), new Matrix(3, 3, [6, 1, 1, 4, -2, 5, 2, 8, 7]),
    new Matrix(3, 3, [2, 3, 4, 0, -3, 5, 0, 0, 7]),
  ], (empty, singleton, swapped, singular, three, upper) => {
    assert.equal(empty.determinant(), 1);
    assert.equal(empty.trace(), 0);
    assert.equal(singleton.determinant(), -7);
    assert.equal(swapped.determinant(), -1);
    assert.equal(singular.determinant(), 0);
    assert.ok(Math.abs(three.determinant() + 306) < 1e-10);
    assert.equal(three.trace(), 11);
    assert.deepEqual(upper.triangular(), { upper: true, lower: false });
    use([upper.transpose()], lower => assert.deepEqual(lower.triangular(), { upper: false, lower: true }));
  });
});

test("selectable determinant methods agree with exact products and preserve operands", () => {
  // This SPD matrix is L L^T for L = [[2,0,0],[1,3,0],[-1,2,4]].
  const values = [4, 2, -2, 2, 10, 5, -2, 5, 21];
  use([new Matrix(), new Matrix(3, 3, values), new Matrix(3, 3, [6, 1, 1, 4, -2, 5, 2, 8, 7])], (empty, spd, general) => {
    for (const algorithm of ["auto", "lu", "cholesky"]) {
      assert.equal(empty.determinant(algorithm), 1);
      assert.ok(Math.abs(spd.determinant(algorithm) - 576) < 1e-10);
      assert.deepEqual([...spd.toArray()], values);
    }
    assert.equal(general.determinant("auto"), -306);
    assert.ok(Math.abs(general.determinant("lu") + 306) < 1e-10);
    assert.throws(() => general.determinant("cholesky"), /symmetric positive definite/);
    for (const algorithm of ["LU", "unknown", 1, null, true, {}, "toString"]) {
      assert.throws(() => spd.determinant(algorithm), /unknown determinant algorithm/);
    }
  });
});

test("Cholesky rejects nonsymmetric, indefinite, and singular matrices explicitly", () => {
  for (const values of [[2, 100, 1, 2], [1, 2, 2, 1], [-1, 0, 0, -1], [1, 1, 1, 1], [2, 1, 1 + Number.EPSILON, 2]]) {
    use([new Matrix(2, 2, values)], a => {
      assert.throws(() => a.determinant("cholesky"), /symmetric positive definite/);
      assert.deepEqual([...a.toArray()], values);
      assert.ok(Number.isFinite(a.determinant("auto")));
    });
  }
  use([new Matrix(2, 3)], a => {
    for (const algorithm of ["lu", "cholesky"]) assert.throws(() => a.determinant(algorithm), /dimensions/);
  });
});

test("scaled determinant products retain range independently of pivot order", () => {
  for (const exponents of [[800, 800, -800, -800], [-800, -800, 800, 800], [800, -800, 800, -800]]) {
    const values = Array(16).fill(0);
    exponents.forEach((exponent, i) => values[i * 4 + i] = 2 ** exponent);
    use([new Matrix(4, 4, values)], a => {
      for (const algorithm of ["auto", "lu", "cholesky"]) assert.equal(a.determinant(algorithm), 1);
      assert.deepEqual([...a.toArray()], values);
    });
  }
  use([new Matrix(1, 1, [Number.MIN_VALUE]), new Matrix(2, 2, [Number.MAX_VALUE, 0, 0, 2])], (tiny, overflowing) => {
    for (const algorithm of ["auto", "lu", "cholesky"]) {
      assert.equal(tiny.determinant(algorithm), Number.MIN_VALUE);
      assert.throws(() => overflowing.determinant(algorithm), /nonfinite/);
    }
  });
  use([new Matrix(2, 2, [4e300, 2, 2, 3e-300])], a => assert.ok(Math.abs(a.determinant("cholesky") - 8) < 1e-13));
});

test("tiny fast paths fall back on overflow, cancellation, and subnormal products", () => {
  for (const [size, values, expected] of [
    [2, [1e200, 1e200, 1e200, 1e200], 0],
    [2, [1, 1, 1, 1 + Number.EPSILON], Number.EPSILON],
    [2, [1e-200, 1e-200, 3e-124, 6e-124], Number.MIN_VALUE],
    [3, [1e-200, 1e-200, 0, 3e-124, 6e-124, 0, 0, 0, 1], Number.MIN_VALUE],
    [3, [1e200, 1e200, 0, 1e200, 1e200, 0, 0, 0, 1], 0],
  ]) use([new Matrix(size, size, values)], a => {
    for (const algorithm of ["auto", "lu"]) assert.equal(a.determinant(algorithm), expected);
  });
  for (const values of [[1e308, 1e308, 1e-308, 2e-308], [1e308, 1e-308, 1e308, 2e-308]]) {
    use([new Matrix(2, 2, values)], a => {
      for (const algorithm of ["auto", "lu"]) assert.ok(Math.abs(a.determinant(algorithm) - 1) < 1e-14);
      assert.deepEqual([...a.toArray()], values);
    });
  }
});

test("empty rectangular dimensions retain shapes and zero inner products", () => {
  use([new Matrix(2, 0), new Matrix(0, 3)], (a, b) => {
    use([a.multiply(b), b.transpose()], (product, transpose) => {
      assert.deepEqual([product.rows, product.cols], [2, 3]);
      assert.deepEqual([...product.toArray()], [0, 0, 0, 0, 0, 0]);
      assert.deepEqual([transpose.rows, transpose.cols], [3, 0]);
      assert.deepEqual([...transpose.row(2)], []);
    });
  });
});

test("invalid dimensions, indices, values, and overflow fail explicitly", () => {
  for (const value of [-1, 0.5, NaN, Infinity, "2", true, 0x100000000]) assert.throws(() => new Matrix(value, 1));
  assert.throws(() => new Matrix(0xffffffff, 2), /overflow/);
  assert.throws(() => new Matrix(2, 2, [1]), /count/);
  for (const value of [NaN, Infinity, -Infinity, true, "1"]) assert.throws(() => new Matrix(1, 1, [value]), /finite/);
  use([new Matrix(1, 1, [1e308])], a => {
    assert.throws(() => a.scale(2), /nonfinite/);
    assert.throws(() => a.scale("1"), /finite/);
    assert.throws(() => a.set(0, 0, NaN), /finite/);
    assert.equal(a.get(0, 0), 1e308);
    for (const index of [-1, 1, 0.5, true, "0"]) {
      assert.throws(() => a.get(index, 0), /bounds/);
      assert.throws(() => a.row(index), /bounds/);
      assert.throws(() => a.column(index), /bounds/);
    }
  });
});

test("disposed handles cannot be reused and double disposal is harmless", () => {
  const a = new Matrix(1, 1, [3]);
  a.dispose(); a.dispose();
  assert.equal(a.disposed, true);
  for (const operation of [() => a.rows, () => a.get(0, 0), () => a.set(0, 0, 1), () => a.copy(), () => a.toArray(), () => a.trace()]) {
    assert.throws(operation, /disposed/);
  }
  use([new Matrix(1, 1)], b => assert.throws(() => b.add(a), /disposed/));
});

test("matrices from separate WebAssembly instances cannot mix handles", async () => {
  const { Matrix: OtherMatrix } = await createMatrixAPI(options);
  use([new Matrix(1, 1), new OtherMatrix(1, 1)], (a, b) => {
    assert.throws(() => a.add(b), /same WebAssembly instance/);
  });
});

test("memory growth preserves existing matrices and exported arrays remain copies", () => {
  use([new Matrix(1, 2, [3, 4])], small => {
    const copied = small.toArray();
    // More than the default initial 16 MiB forces the Emscripten heap to grow.
    use([new Matrix(1, 3_000_000)], large => {
      large.set(0, 2_999_999, 7);
      assert.equal(large.get(0, 2_999_999), 7);
      assert.equal(small.get(0, 1), 4);
      small.set(0, 1, 8);
      assert.equal(small.get(0, 1), 8);
      assert.deepEqual([...copied], [3, 4]);
    });
  });
});
