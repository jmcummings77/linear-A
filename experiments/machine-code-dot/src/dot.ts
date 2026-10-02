import { performance } from "node:perf_hooks";
function dot(a: Float64Array, b: Float64Array, offset: number, n: number): number {
    let sum = 0;
    for (let i = 0; i < n; i++) sum += a[offset + i] * b[offset + i];
    return sum;
}
const args = process.argv.slice(2).map(Number);
if (args.length !== 3 || args.some(x => !Number.isSafeInteger(x))) throw Error("expected n iterations seed");
const [n, iterations, seed] = args;
if (n < 0 || n > 1048576 || iterations < 1 || iterations > 100000000 || seed < 0 || seed > 1000000) throw Error("out of range");
const count = Math.max(n, 1) * 32;
const a = new Float64Array(count), b = new Float64Array(count);
for (let i = 0; i < count; i++) {
    a[i] = ((i * 17 + seed * 13) % 101 - 50) / 16;
    b[i] = ((i * 29 + (seed + 1) * 7) % 103 - 51) / 16;
}
let sink = 0;
for (let i = 0; i < Math.max(8, Math.min(iterations, 128)); i++) sink += dot(a, b, (i % 32) * n, n);
let checksum = 0;
const start = performance.now();
for (let i = 0; i < iterations; i++) checksum += dot(a, b, (i % 32) * n, n);
const elapsed = (performance.now() - start) * 1e6;
if (!Number.isFinite(checksum) || !Number.isFinite(sink)) throw Error("nonfinite result");
console.log(JSON.stringify({ elapsed_ns: elapsed, iterations, checksum }));
