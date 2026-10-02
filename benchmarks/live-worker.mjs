/** Isolated, terminable browser worker for real WASM checks and measurements. */
const OPERATIONS = new Set(["add", "subtract", "scale", "transpose", "multiply", "cross", "rotation2d", "rotation3d", "trace", "eigen_symmetric", "eigen_general", "determinant", "determinant_lu", "determinant_spd_lu", "determinant_cholesky"]);
const SIZES = new Set([16, 48, 64, 128, 256]);
const DETERMINANT_SIZES = new Set([4, 8, 16, 24, 48]);
const EIGEN_SIZES = new Set([4, 8, 16, 32, 48]);
const MAX_ITERATIONS = 4096;
const TARGET_MS = 30;
const MAX_RUN_MS = 15_000;

function finite(value) { return typeof value === "number" && Number.isFinite(value); }
function close(actual, expected, relative = 1e-9, absolute = 1e-10) {
  return finite(actual) && finite(expected)
    && Math.abs(actual - expected) <= Math.max(absolute, relative * Math.max(Math.abs(actual), Math.abs(expected)));
}
function yieldBatch() { return new Promise(resolve => setTimeout(resolve, 0)); }
function errorMessage(error) { return error instanceof Error ? error.message : String(error); }

function operation(name, a, b, scalar = 1.25) {
  switch (name) {
    case "add": return a.add(b);
    case "subtract": return a.subtract(b);
    case "multiply": return a.multiply(b);
    case "cross": return a.cross(b);
    case "rotation2d":
      if (a.rows !== 0 || a.cols !== 0) throw new RangeError("rotation2d expects an empty 0-by-0 input");
      return a.constructor.rotation2D(scalar);
    case "rotation3d": return a.constructor.rotationAxisAngle(a, scalar);
    case "scale": return a.scale(scalar);
    case "transpose": return a.transpose();
    case "trace": return a.trace();
    case "eigen_symmetric": return a.eigenSymmetric();
    case "eigen_general": return a.eigenGeneral();
    case "determinant": return a.determinant();
    case "determinant_lu": case "determinant_spd_lu": return a.determinant("lu");
    case "determinant_cholesky": return a.determinant("cholesky");
    case "triangular": return a.triangular();
    default: throw new RangeError("unsupported operation");
  }
}

function assertResult(actual, expected) {
  if (!actual || typeof actual !== "object" || Array.isArray(actual)) throw new Error("result must be an object");
  if (!expected || typeof expected !== "object" || Array.isArray(expected)) throw new Error("fixture expected result must be an object");
  const keys = Object.keys(expected).sort();
  if (JSON.stringify(Object.keys(actual).sort()) !== JSON.stringify(keys)) throw new Error("result keys differ from the independent fixture");
  for (const key of keys) {
    const value = actual[key], wanted = expected[key];
    if (key === "rows" || key === "cols") {
      if (!Number.isSafeInteger(value) || value < 0 || value !== wanted) throw new Error(`${key} is not the expected nonnegative integer`);
    } else if (key === "upper" || key === "lower") {
      if (typeof value !== "boolean" || value !== wanted) throw new Error(`${key} is not the expected boolean`);
    } else if (key === "values") {
      if (!Array.isArray(value) || !Array.isArray(wanted) || value.length !== wanted.length) throw new Error("wrong result value count or type");
      for (let i = 0; i < value.length; i++) {
        if (!close(value[i], wanted[i])) throw new Error(`values[${i}] differs from the independent fixture`);
      }
    } else if (key === "value") {
      if (!close(value, wanted)) throw new Error("scalar differs from the independent fixture");
    } else {
      throw new Error("unsupported expected result field");
    }
  }
}

/** Verify AQ=QΛ and QᵀQ=I; neither signs nor repeated-eigenspace bases are fixed. */
export function assertEigenResult(actual, input, expectedValues) {
  const n = input.rows, tolerance = 1e-9;
  if (!actual || Object.keys(actual).sort().join() !== "eigenvalues,eigenvectors") throw new Error("invalid eigenpair result keys");
  const values = actual.eigenvalues, vectors = actual.eigenvectors;
  if (!Array.isArray(values) || values.length !== n || !values.every(finite)) throw new Error("invalid eigenvalue array");
  if (!vectors || Object.keys(vectors).sort().join() !== "cols,rows,values"
    || !Number.isSafeInteger(vectors.rows) || !Number.isSafeInteger(vectors.cols)
    || vectors.rows !== n || vectors.cols !== n || !Array.isArray(vectors.values)
    || vectors.values.length !== n*n || !vectors.values.every(finite)) throw new Error("invalid eigenvector matrix");
  const q = vectors.values;
  if (values.some((value, i) => i > 0 && value < values[i-1])) throw new Error("eigenvalues are not sorted");
  if (q.some(value => Math.abs(value) > 1 + tolerance*n)) throw new Error("eigenvectors are not normalized");
  let orthogonality = 0, scale = 0;
  for (const value of input.values) scale = Math.max(scale, Math.abs(value));
  for (const value of values) scale = Math.max(scale, Math.abs(value));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let dot = 0;
    for (let k = 0; k < n; k++) dot += q[k*n+i]*q[k*n+j];
    orthogonality = Math.hypot(orthogonality, dot - (i === j ? 1 : 0));
  }
  if (orthogonality > tolerance*n) throw new Error("eigenvectors are not orthonormal");
  let norm = 0, residual = 0;
  if (scale) {
    for (const value of input.values) norm = Math.hypot(norm, value / scale);
    for (let row = 0; row < n; row++) for (let col = 0; col < n; col++) {
      let av = 0;
      for (let k = 0; k < n; k++) av += (input.values[row*n+k]/scale)*q[k*n+col];
      residual = Math.hypot(residual, av - q[row*n+col]*(values[col]/scale));
    }
    if (residual > tolerance*norm) throw new Error("eigenpair residual exceeds tolerance");
  }
  if (expectedValues !== undefined) {
    if (!Array.isArray(expectedValues) || expectedValues.length !== n || !expectedValues.every(finite)) throw new Error("invalid reference spectrum");
    for (let i = 0; i < n; i++) {
      const error = scale ? Math.abs(values[i]/scale - expectedValues[i]/scale) : Math.abs(values[i]-expectedValues[i]);
      if (error > tolerance*norm) throw new Error("eigenvalue differs from the independent spectrum");
    }
  }
  return { relative_residual: norm ? residual/norm : 0, orthogonality_error: orthogonality };
}

/** Validate complex right columns, including dependent columns for defective inputs. */
export function assertGeneralEigenResult(actual,input,expected,spectrumScale,componentwise=false) {
  const n=input.rows,tolerance=1e-9;
  if (!actual || input.cols!==n || Object.keys(actual).sort().join()!==
      'eigenvalues_imag,eigenvalues_real,eigenvectors_imag,eigenvectors_real') throw new Error('invalid complex eigenpair result');
  const vr=actual.eigenvalues_real,vi=actual.eigenvalues_imag;
  for(const values of [vr,vi]) if(!Array.isArray(values)||values.length!==n||!values.every(finite)) throw new Error('invalid complex eigenvalue array');
  for(const matrix of [actual.eigenvectors_real,actual.eigenvectors_imag]) {
    if(!matrix||Object.keys(matrix).sort().join()!=='cols,rows,values'||matrix.rows!==n||matrix.cols!==n
      ||!Number.isSafeInteger(matrix.rows)||!Number.isSafeInteger(matrix.cols)||!Array.isArray(matrix.values)
      ||matrix.values.length!==n*n||!matrix.values.every(finite)) throw new Error('invalid complex eigenvector matrix');
  }
  const qr=actual.eigenvectors_real.values,qi=actual.eigenvectors_imag.values;
  let scale=0,norm=0,residual=0;
  for(const x of [...input.values,...vr,...vi])scale=Math.max(scale,Math.abs(x));
  scale ||= 1;
  for(const x of input.values)norm=Math.hypot(norm,x/scale);
  for(let col=0;col<n;col++) {
    if(col && (vr[col]<vr[col-1] || vr[col]===vr[col-1]&&vi[col]<vi[col-1])) throw new Error('complex eigenvalues are not sorted');
    let length=0,error=0;
    for(let row=0;row<n;row++) {
      const re=qr[row*n+col],im=qi[row*n+col];length=Math.hypot(length,re,im);
      let ar=0,ai=0;
      for(let k=0;k<n;k++) {ar+=input.values[row*n+k]/scale*qr[k*n+col];ai+=input.values[row*n+k]/scale*qi[k*n+col];}
      error=Math.hypot(error,ar-(vr[col]/scale*re-vi[col]/scale*im),ai-(vr[col]/scale*im+vi[col]/scale*re));
    }
    if(componentwise) for(let row=0;row<n;row++) {
      let ar=0,ai=0,bound=0;
      for(let k=0;k<n;k++) {const re=input.values[row*n+k]*qr[k*n+col],im=input.values[row*n+k]*qi[k*n+col];ar+=re;ai+=im;bound+=Math.hypot(re,im);}
      const re=vr[col]*qr[row*n+col]-vi[col]*qi[row*n+col],im=vr[col]*qi[row*n+col]+vi[col]*qr[row*n+col];
      if(Math.hypot(ar-re,ai-im)>tolerance*(bound+Math.hypot(re,im)))throw new Error('componentwise eigenpair residual exceeds tolerance');
    }
    if(Math.abs(length-1)>tolerance*Math.max(n,1))throw new Error('complex eigenvectors are not normalized');
    if(error>tolerance*Math.max(norm,Math.hypot(vr[col]/scale,vi[col]/scale)))throw new Error('complex eigenpair residual exceeds tolerance');
    residual=Math.max(residual,error/(norm||1));
  }
  if(expected!==undefined) {
    if(!Array.isArray(expected)||expected.length!==n||!expected.every(x=>Array.isArray(x)&&x.length===2&&x.every(finite)))throw new Error('invalid complex reference spectrum');
    const pending=expected.map(x=>[...x]), divisor=spectrumScale??scale;
    for(let i=0;i<n;i++) {
      let closest=0,error=Infinity;
      pending.forEach(([re,im],j)=>{const d=Math.hypot(vr[i]/divisor-re/divisor,vi[i]/divisor-im/divisor);if(d<error){error=d;closest=j;}});
      if(error>tolerance*(spectrumScale===undefined?Math.max(norm,1):1))throw new Error('complex eigenvalue differs from the independent spectrum');
      pending.splice(closest,1);
    }
  }
  return {relative_residual:residual};
}

export function generalEigenObject(result,Matrix) {
  if(!result||!(result.vectorsReal instanceof Matrix)||!(result.vectorsImag instanceof Matrix))throw new Error('invalid WASM complex eigenpair result');
  return {eigenvalues_real:Array.from(result.valuesReal),eigenvalues_imag:Array.from(result.valuesImag),
    eigenvectors_real:{rows:result.vectorsReal.rows,cols:result.vectorsReal.cols,values:Array.from(result.vectorsReal.toArray())},
    eigenvectors_imag:{rows:result.vectorsImag.rows,cols:result.vectorsImag.cols,values:Array.from(result.vectorsImag.toArray())}};
}
function generalEigenSpectrum(size,seed) {
  return Array.from({length:size},(_,i)=>[1+(seed%17)/16+Math.floor(i/2)/8,
    i===size-1&&size%2?0:(.5+Math.floor(i/2)/16)*(i%2?1:-1)]);
}
function generalEigenInput(i,j,seed) {
  const k=Math.floor(i/2), re=1+(seed%17)/16+k/8,im=.5+k/16;
  return i===j?re:i%2===0&&j===i+1?-im:i%2===1&&j===i-1?im:i<j?((i*3+j*5+seed)%11-5)/32:0;
}

function eigenObject(result, Matrix) {
  if (!result || !(result.vectors instanceof Matrix)
    || (!Array.isArray(result.values) && !(result.values instanceof Float64Array))) throw new Error("invalid WASM eigenpair result");
  return { eigenvalues: Array.from(result.values), eigenvectors: {
    rows: result.vectors.rows, cols: result.vectors.cols, values: Array.from(result.vectors.toArray()),
  } };
}

function eigenSpectrum(size, seed) {
  return Array.from({length:size}, (_, i) => 2 + (seed % 17)/16 - 2*Math.cos((i+1)*Math.PI/(size+1)));
}

// Independent quaternion action on basis vectors, rather than Rodrigues' formula.
function rotationReference(axis, angle, dimensions) {
  const largest = Math.max(...axis.map(Math.abs)), scaled = axis.map(x => x / largest);
  const norm = Math.hypot(...scaled), sine = Math.sin(angle / 2);
  const q = [Math.cos(angle / 2), ...scaled.map(x => sine*x/norm)], conjugate = [q[0], -q[1], -q[2], -q[3]];
  const product = ([w,x,y,z], [r,i,j,k]) => [w*r-x*i-y*j-z*k, w*i+x*r+y*k-z*j, w*j-x*k+y*r+z*i, w*k+x*j-y*i+z*r];
  const columns = [0,1,2].map(col => product(product(q, [0, ...[0,1,2].map(row => Number(row === col))]), conjugate).slice(1));
  return { rows: dimensions, cols: dimensions,
    values: Array.from({length:dimensions*dimensions}, (_, i) => columns[i % dimensions][Math.floor(i / dimensions)]) };
}

function vectorWorkload(op, seed) {
  const a = Array.from({length:3}, (_, i) => numerator(i,seed)/16);
  const b = Array.from({length:3}, (_, i) => numerator(i,seed+1)/16);
  b[2] = -b[2];
  if (op === "cross") return {a:{rows:3,cols:1,values:a}, b:{rows:3,cols:1,values:b}, scalar:.5,
    expected:{rows:3,cols:1,values:[a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]}};
  const dimensions = op === "rotation2d" ? 2 : 3;
  return {a:dimensions===2?{rows:0,cols:0,values:[]}:{rows:3,cols:1,values:[1,2,3]}, b:null, scalar:.5,
    expected:rotationReference(dimensions===2?[0,0,1]:[1,2,3],.5,dimensions)};
}

/** Run the independently generated arithmetic and eigen fixtures against actual WASM. */
export async function runChecks(Matrix, fixtures, onProgress = () => {}) {
  if (!Array.isArray(fixtures) || fixtures.length === 0 || fixtures.length > 500) throw new RangeError("expected a bounded, nonempty fixture array");
  const checks = [];
  for (const fixture of fixtures) {
    const check = { name: String(fixture.name), passed: false };
    let a, b, result;
    try {
      if (!OPERATIONS.has(fixture.op) && fixture.op !== "triangular") throw new Error("fixture has an unsupported operation");
      if (typeof fixture.invalid !== "boolean") throw new Error("fixture invalid flag must be boolean");
      a = new Matrix(fixture.a.rows, fixture.a.cols, fixture.a.values);
      if (fixture.b !== null && fixture.b !== undefined) b = new Matrix(fixture.b.rows, fixture.b.cols, fixture.b.values);
      if (fixture.invalid) {
        let rejected = false;
        try { result = operation(fixture.op, a, b, fixture.scalar); } catch { rejected = true; }
        if (!rejected) throw new Error("invalid shape was accepted");
      } else {
        result = operation(fixture.op, a, b, fixture.scalar);
        if (fixture.op === "eigen_general") {
          assertGeneralEigenResult(generalEigenObject(result,Matrix),fixture.a,fixture.expected_complex_eigenvalues,fixture.spectrum_scale,fixture.componentwise);
        } else if (fixture.op === "eigen_symmetric") {
          assertEigenResult(eigenObject(result, Matrix), fixture.a, fixture.expected_eigenvalues);
        } else {
          const actual = result instanceof Matrix
            ? { rows: result.rows, cols: result.cols, values: Array.from(result.toArray()) }
            : typeof result === "object" && result !== null ? result : { value: result };
          assertResult(actual, fixture.expected);
        }
      }
      check.passed = true;
    } catch (error) {
      check.error = errorMessage(error);
    } finally {
      if (result instanceof Matrix) result.dispose();
      if (result?.vectors instanceof Matrix) result.vectors.dispose();
      result?.vectorsReal?.dispose(); result?.vectorsImag?.dispose();
      a?.dispose(); b?.dispose();
    }
    checks.push(check);
    onProgress({ message: `Checking ${checks.length}/${fixtures.length}: ${check.name}`, completed: checks.length, total: fixtures.length, check });
    await yieldBatch();
  }
  return { checks, passed: checks.filter(check => check.passed).length, total: checks.length };
}

function numerator(index, seed) {
  return ((index % 101) * 17 + (seed % 101) * 13) % 101 - 50;
}

// Exact fraction-free elimination on integer numerators. This independent
// oracle does not call WASM's float64 partial-pivot elimination.
function determinantChecksum(size, seed, spd = false) {
  const work = Array.from({ length: size }, (_, row) => Array.from({ length: size }, (_, col) =>
    BigInt(spd ? numerator(row * size + col, seed) + numerator(col * size + row, seed) + (row === col ? size * 128 : 0)
      : numerator(row * size + col, seed) + (row === col ? size * 64 : 0))));
  let previous = 1n, sign = 1n;
  for (let k = 0; k < size - 1; k++) {
    let pivot = k;
    while (pivot < size && work[pivot][k] === 0n) pivot++;
    if (pivot === size) return 0;
    if (pivot !== k) { [work[k], work[pivot]] = [work[pivot], work[k]]; sign = -sign; }
    const diagonal = work[k][k];
    for (let row = k + 1; row < size; row++) {
      for (let col = k + 1; col < size; col++) {
        const value = work[row][col] * diagonal - work[row][k] * work[k][col];
        if (value % previous !== 0n) throw new Error("independent determinant oracle division was not exact");
        work[row][col] = value / previous;
      }
      work[row][k] = 0n;
    }
    previous = diagonal;
  }
  return Number(sign * work[size - 1][size - 1]) / (spd ? 32 : 16) ** size;
}

function expectedChecksum(op, size, seed) {
  if (["cross", "rotation2d", "rotation3d"].includes(op)) {
    const values = vectorWorkload(op,seed).expected.values;
    return values[0]+values[Math.floor(values.length/2)]+values[values.length-1];
  }
  if (op === "eigen_general") return generalEigenSpectrum(size,seed).reduce((sum,[re,im],i)=>sum+(i+1)*(re+Math.abs(im)),size);
  if (op === "eigen_symmetric") return eigenSpectrum(size, seed).reduce((sum, value, i) => sum + (i+1)*value, size);
  if (op.startsWith("determinant")) return determinantChecksum(size, seed, op === "determinant_spd_lu" || op === "determinant_cholesky");
  if (op === "trace") {
    let total = 0;
    for (let row = 0; row < size; row++) total += numerator(row * size + row, seed);
    return total / 16;
  }
  let result = 0;
  for (const index of [0, Math.floor(size * size / 2), size * size - 1]) {
    const row = Math.floor(index / size), col = index % size;
    if (op === "add") result += (numerator(index, seed) + numerator(index, seed + 1)) / 16;
    else if (op === "subtract") result += (numerator(index, seed) - numerator(index, seed + 1)) / 16;
    else if (op === "scale") result += numerator(index, seed) * 5 / 64;
    else if (op === "transpose") result += numerator(col * size + row, seed) / 16;
    else if (op === "multiply") {
      // These small integer sums are exactly representable by JavaScript numbers.
      let product = 0;
      for (let k = 0; k < size; k++) product += numerator(row * size + k, seed) * numerator(k * size + col, seed + 1);
      result += product / 256;
    }
  }
  return result;
}

function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  return sorted[Math.floor(sorted.length / 2)];
}

/** Browser-local samples; measured separately from the saved native comparison. */
export async function runBenchmark(Matrix, config, onProgress = () => {}) {
  if (!config || typeof config !== "object" || Array.isArray(config)) throw new RangeError("benchmark configuration must be an object");
  const { operation: op, size, samples, seed } = config;
  if (!OPERATIONS.has(op)) throw new RangeError("unsupported benchmark operation");
  const allowedSizes = op === "rotation2d" ? new Set([2]) : op === "cross" || op === "rotation3d" ? new Set([3])
    : op.startsWith("determinant") ? DETERMINANT_SIZES : op.startsWith("eigen_") ? EIGEN_SIZES : SIZES;
  if (!allowedSizes.has(size)) throw new RangeError("unsupported benchmark size");
  if (samples !== 3 && samples !== 5) throw new RangeError("sample count must be 3 or 5");
  if (!Number.isSafeInteger(seed) || seed < 0 || seed > 2147483646) throw new RangeError("seed must be an integer in 0..2147483646");
  const deadline = performance.now() + MAX_RUN_MS;
  const bounded = () => { if (performance.now() > deadline) throw new Error("live benchmark exceeded its time budget; try a smaller size"); };
  onProgress({ message: "Preparing inputs and independent checksum", completed: 0, total: samples });
  await yieldBatch();
  const expected = expectedChecksum(op, size, seed);
  if (!finite(expected)) throw new Error("independent expected checksum is nonfinite");
  const left = new Float64Array(size * size), right = new Float64Array(size * size);
  for (let index = 0; index < left.length; index++) {
    left[index] = numerator(index, seed) / 16;
    right[index] = numerator(index, seed + 1) / 16;
  }
  if (op === "determinant_spd_lu" || op === "determinant_cholesky") {
    for (let row = 0; row < size; row++) for (let col = 0; col <= row; col++) {
      const value = (left[row * size + col] + left[col * size + row]) / 2;
      left[row * size + col] = left[col * size + row] = value;
    }
  }
  if (op.startsWith("determinant")) for (let i = 0; i < size; i++) left[i * size + i] += size * 4;
  if (op === "eigen_symmetric") {
    for (let row = 0; row < size; row++) for (let col = 0; col < size; col++) {
      left[row*size+col] = row === col ? 2 + (seed % 17)/16 : Math.abs(row-col) === 1 ? -1 : 0;
    }
  }
  if (op === "eigen_general") for(let i=0;i<size;i++) for(let j=0;j<size;j++) left[i*size+j]=generalEigenInput(i,j,seed);
  let a, b;
  try {
    const vector = ["cross", "rotation2d", "rotation3d"].includes(op) ? vectorWorkload(op, seed) : null;
    a = vector ? new Matrix(vector.a.rows, vector.a.cols, vector.a.values) : new Matrix(size, size, left);
    b = vector ? vector.b ? new Matrix(vector.b.rows, vector.b.cols, vector.b.values) : null : new Matrix(size, size, right);
    if (vector) {
      const probe = operation(op, a, b, vector.scalar);
      try { assertResult({rows:probe.rows, cols:probe.cols, values:Array.from(probe.toArray())}, vector.expected); }
      finally { probe?.dispose(); }
    }
    if (op === "eigen_general") {
      const probe=a.eigenGeneral();
      try {assertGeneralEigenResult(generalEigenObject(probe,Matrix),{rows:size,cols:size,values:Array.from(left)},generalEigenSpectrum(size,seed));}
      finally {probe?.vectorsReal?.dispose();probe?.vectorsImag?.dispose();}
    }
    if (op === "eigen_symmetric") {
      // Validate the complete eigenbasis of this exact timed input outside the clock.
      const probe = a.eigenSymmetric();
      try { assertEigenResult(eigenObject(probe, Matrix), {rows:size, cols:size, values:Array.from(left)}, eigenSpectrum(size,seed)); }
      finally { probe?.vectors?.dispose(); }
    }
    function consume() {
      const result = operation(op, a, b, vector ? vector.scalar : 1.25);
      if (result instanceof Matrix) {
        try {
          const checksum = result.checksum();
          if (!finite(checksum)) throw new Error("WASM returned a nonfinite or nonnumeric checksum");
          return checksum;
        } finally { result.dispose(); }
      }
      if (op === "eigen_general") {
        try {
          const pair=generalEigenObject(result,Matrix);
          let sum=0;
          for(let i=0;i<pair.eigenvalues_real.length;i++)sum+=(i+1)*(pair.eigenvalues_real[i]+Math.abs(pair.eigenvalues_imag[i]));
          for(const x of pair.eigenvectors_real.values)sum+=x*x;
          for(const x of pair.eigenvectors_imag.values)sum+=x*x;
          if(!finite(sum))throw new Error('nonfinite complex eigenpair checksum');
          return sum;
        } finally {result?.vectorsReal?.dispose();result?.vectorsImag?.dispose();}
      }
      if (op === "eigen_symmetric") {
        try {
          const pair = eigenObject(result, Matrix);
          let sum = 0;
          for (let i = 0; i < pair.eigenvalues.length; i++) sum += (i+1)*pair.eigenvalues[i];
          for (const value of pair.eigenvectors.values) sum += value*value;
          if (!finite(sum)) throw new Error("nonfinite eigenpair checksum");
          return sum;
        } finally { result?.vectors?.dispose(); }
      }
      if (!finite(result)) throw new Error("WASM returned a nonfinite or nonnumeric scalar");
      return result;
    }
    function measure(iterations) {
      bounded();
      for (let i = 0; i < Math.max(5, Math.min(iterations, 100)); i++) consume();
      let checksum = 0;
      const start = performance.now();
      for (let i = 0; i < iterations; i++) checksum += consume();
      const elapsed_ms = performance.now() - start;
      if (!finite(elapsed_ms) || elapsed_ms < 0) throw new Error("browser clock produced an invalid elapsed time");
      if (!close(checksum, expected * iterations, 1e-7, 1e-8 * iterations)) throw new Error("live benchmark checksum differs from the independent reference");
      bounded();
      return { elapsed_ms, ns_per_op: elapsed_ms * 1e6 / iterations, checksum };
    }

    onProgress({ message: "Warming up and calibrating the WASM workload", completed: 0, total: samples });
    let probeIterations = 1, probe = measure(probeIterations);
    while (probe.elapsed_ms === 0 && probeIterations < MAX_ITERATIONS) {
      probeIterations = Math.min(MAX_ITERATIONS, probeIterations * 8);
      await yieldBatch();
      probe = measure(probeIterations);
    }
    if (probe.elapsed_ms === 0) throw new Error("browser timer resolution is too coarse for this workload");
    let iterations = Math.max(1, Math.min(MAX_ITERATIONS, Math.round(TARGET_MS * 1e6 / probe.ns_per_op)));
    await yieldBatch();
    probe = measure(iterations);
    if (probe.ns_per_op > 0) iterations = Math.max(1, Math.min(MAX_ITERATIONS, Math.round(TARGET_MS * 1e6 / probe.ns_per_op)));
    const measured = [];
    for (let index = 0; index < samples; index++) {
      await yieldBatch();
      const sample = measure(iterations);
      if (sample.elapsed_ms === 0) throw new Error("browser timer resolution is too coarse for this workload");
      measured.push(sample);
      onProgress({ message: `Measured sample ${index + 1}/${samples}`, completed: index + 1, total: samples });
    }
    const times = measured.map(sample => sample.ns_per_op);
    const median_ns = median(times);
    return { operation: op, size, seed, iterations, samples: measured, median_ns,
      mad_ns: median(times.map(value => Math.abs(value - median_ns))),
      min_ns: Math.min(...times), max_ns: Math.max(...times), checksum_verified: true };
  } finally {
    a?.dispose(); b?.dispose();
  }
}

async function initialize(bundle) {
  if (!bundle || typeof bundle.wrapper_source !== "string" || typeof bundle.module_source !== "string" || typeof bundle.wasm_base64 !== "string") {
    throw new Error("report does not contain an embedded WASM bundle");
  }
  const wrapperURL = URL.createObjectURL(new Blob([bundle.wrapper_source], { type: "text/javascript" }));
  const glueURL = URL.createObjectURL(new Blob([bundle.module_source], { type: "text/javascript" }));
  try {
    const { createMatrixAPI } = await import(wrapperURL);
    const bytes = Uint8Array.from(atob(bundle.wasm_base64), character => character.charCodeAt(0));
    return await createMatrixAPI({ moduleUrl: glueURL, wasmBinary: bytes, locateFile: () => "embedded.wasm" });
  } finally {
    URL.revokeObjectURL(wrapperURL);
    URL.revokeObjectURL(glueURL);
  }
}

if (typeof WorkerGlobalScope !== "undefined" && globalThis instanceof WorkerGlobalScope) {
  let started = false;
  globalThis.onmessage = async ({ data }) => {
    if (started) return;
    started = true;
    let checks;
    try {
      if (!data || (data.type !== "checks" && data.type !== "benchmark")) throw new Error("unsupported live request");
      const progress = detail => globalThis.postMessage({ type: "progress", ...detail });
      progress({ message: "Compiling embedded WebAssembly", completed: 0, total: 0 });
      const { Matrix } = await initialize(data.bundle);
      checks = await runChecks(Matrix, data.bundle.fixtures, progress);
      if (data.type === "checks") {
        globalThis.postMessage({ type: "done", kind: "checks", result: checks });
      } else {
        if (checks.passed !== checks.total) throw new Error("shared correctness checks failed; benchmark was not run");
        const result = await runBenchmark(Matrix, data.config, progress);
        globalThis.postMessage({ type: "done", kind: "benchmark", result, checks });
      }
    } catch (error) {
      globalThis.postMessage({ type: "error", message: errorMessage(error), ...(checks ? { checks } : {}) });
    }
  };
}
