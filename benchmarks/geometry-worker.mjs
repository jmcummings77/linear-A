/** Verified 3x3 geometry and genuine, separately instrumented WASM arithmetic. */
const OPERATIONS = new Set(["multiply", "add", "subtract", "scale", "transpose", "cross", "rotate"]);
const LIMIT = 10;
const TRACE_SOURCE = "ports/c/matrix.c";
const CALL_PATH = ["Matrix.multiplyWithTrace", "wm_multiply_trace", "m_multiply"];
const finite = value => typeof value === "number" && Number.isFinite(value);
const yieldTurn = () => new Promise(resolve => setTimeout(resolve, 0));

function close(actual, expected, absolute = 1e-10) {
  return finite(actual) && finite(expected)
    && Math.abs(actual - expected) <= Math.max(absolute, 1e-9 * Math.max(Math.abs(actual), Math.abs(expected)));
}

function entries(value, name) {
  if ((!Array.isArray(value) && !(value instanceof Float64Array)) || value.length !== 9) {
    throw new RangeError(`${name} must contain exactly nine row-major values`);
  }
  return Array.from(value, entry => {
    if (!finite(entry) || Math.abs(entry) > LIMIT) throw new RangeError(`${name} entries must be finite numbers in [-10, 10]`);
    return entry;
  });
}

function multiplyReference(left, right) {
  const result = Array(9).fill(0);
  for (let row = 0; row < 3; row++) {
    for (let col = 0; col < 3; col++) {
      for (let k = 0; k < 3; k++) result[row * 3 + col] += left[row * 3 + k] * right[k * 3 + col];
    }
  }
  return result;
}

function crossReference(a, b) {
  return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}

// Quaternion action on basis vectors is independent of the WASM Rodrigues implementation.
function rotationReference(axis, radians) {
  const largest = Math.max(...axis.map(Math.abs)), scaled = axis.map(x=>x/largest), norm = Math.hypot(...scaled);
  const sine = Math.sin(radians/2), q = [Math.cos(radians/2), ...scaled.map(x=>sine*x/norm)];
  const product = ([w,x,y,z], [r,i,j,k]) => [w*r-x*i-y*j-z*k, w*i+x*r+y*k-z*j, w*j-x*k+y*r+z*i, w*k+x*j-y*i+z*r];
  const conjugate = [q[0],-q[1],-q[2],-q[3]];
  const columns = [0,1,2].map(col=>product(product(q,[0,...[0,1,2].map(row=>Number(row===col))]),conjugate).slice(1));
  return Array.from({length:9},(_,i)=>columns[i%3][Math.floor(i/3)]);
}

function reference(operation, m, n, scalar, axis, radians) {
  if (operation === "rotate") return multiplyReference(rotationReference(axis,radians),m);
  if (operation === "cross") {
    const columns = [0,1,2].map(col=>crossReference([m[col],m[3+col],m[6+col]],axis));
    return Array.from({length:9},(_,i)=>columns[i%3][Math.floor(i/3)]);
  }
  if (operation === "multiply") return multiplyReference(m, n);
  if (operation === "add") return m.map((value, index) => value + n[index]);
  if (operation === "subtract") return m.map((value, index) => value - n[index]);
  if (operation === "scale") return m.map(value => value * scalar);
  return Array.from({ length: 9 }, (_, index) => m[(index % 3) * 3 + Math.floor(index / 3)]);
}

function determinantReference(a) {
  // The six permutation products are independent of the WASM pivoting algorithm.
  const terms = [a[0] * a[4] * a[8], a[1] * a[5] * a[6], a[2] * a[3] * a[7],
    -a[2] * a[4] * a[6], -a[1] * a[3] * a[8], -a[0] * a[5] * a[7]];
  return { value: terms.reduce((sum, value) => sum + value, 0),
    tolerance: Math.max(1e-10, 64 * Number.EPSILON * terms.reduce((sum, value) => sum + Math.abs(value), 0)) };
}

function applyReference(matrix, point) {
  return [0, 1, 2].map(row => matrix[row * 3] * point[0] + matrix[row * 3 + 1] * point[1] + matrix[row * 3 + 2] * point[2]);
}

function grid(dimensions) {
  const points = [];
  const coordinates = dimensions.map(count => Array.from({length:count}, (_, i) => count === 1 ? 0 : 2*i/(count-1)-1));
  for (const x of coordinates[0]) {
    for (const y of coordinates[1]) {
      for (const z of coordinates[2]) points.push([x, y, z]);
    }
  }
  return points;
}

/**
 * Matrix inputs and results use row-major storage; geometric vectors are columns.
 * Vectors are A*x, anchored at x (arrow head x+A*x), not displacement A*x-x.
 * Rotation playback samples R(axis,t*angle)*M; other operations blend M and R.
 */
export async function computeGeometry(Matrix, config, { TraceMatrix, onProgress = () => {} } = {}) {
  if (typeof Matrix !== "function") throw new TypeError("a WASM Matrix class is required");
  if (!config || typeof config !== "object" || Array.isArray(config)) throw new RangeError("geometry configuration must be an object");
  const operation = config.operation;
  if (!OPERATIONS.has(operation)) throw new RangeError("unsupported geometry operation");
  const m = entries(config.m, "M"), n = entries(config.n, "N");
  const dimensions = config.grid === undefined ? [5,5,5] : Array.isArray(config.grid) ? Array.from(config.grid) : null;
  if (!Array.isArray(dimensions) || dimensions.length !== 3
    || !dimensions.every(count => Number.isInteger(count) && count >= 1 && count <= 15)) {
    throw new RangeError("grid dimensions must contain three integers from 1 to 15");
  }
  const scalar = config.scalar === undefined ? 1.25 : config.scalar;
  if (!finite(scalar) || Math.abs(scalar) > LIMIT) throw new RangeError("scalar must be a finite number in [-10, 10]");
  let axis = null, radians = 0;
  if (operation === "rotate" || operation === "cross") {
    if (!Array.isArray(config.axis) || config.axis.length !== 3
      || !config.axis.every(x=>finite(x)&&Math.abs(x)<=LIMIT)) throw new RangeError("direction must contain three finite entries in [-10, 10]");
    axis = config.axis.slice();
    if (operation === "rotate") {
      if (!axis.some(x=>x!==0)) throw new RangeError("rotation axis must be nonzero");
      radians = config.radians;
      if (!finite(radians) || Math.abs(radians)>2*Math.PI) throw new RangeError("rotation angle must be finite and within one turn");
    }
  }
  if (typeof onProgress !== "function") throw new TypeError("onProgress must be a function");
  const owned = [];
  let verified = 0;
  const own = matrix => {
    if (!matrix || typeof matrix.dispose !== "function") throw new Error("WASM operation did not return an owned matrix");
    owned.push(matrix);
    return matrix;
  };
  const check = (actual, expected, label, absolute) => {
    if (!close(actual, expected, absolute)) throw new Error(`${label} differs from the independent reference`);
    verified++;
  };
  const checkArray = (actual, expected, label) => {
    if ((!Array.isArray(actual) && !(actual instanceof Float64Array)) || actual.length !== expected.length) {
      throw new Error(`${label} has an invalid value count or type`);
    }
    for (let index = 0; index < expected.length; index++) check(actual[index], expected[index], `${label}[${index}]`);
  };
  const values = (matrix, rows, cols, label) => {
    if (matrix.rows !== rows || matrix.cols !== cols) throw new Error(`${label} has an invalid shape`);
    const result = matrix.toArray();
    if ((!Array.isArray(result) && !(result instanceof Float64Array)) || result.length !== rows * cols) {
      throw new Error(`${label} has an invalid value count or type`);
    }
    return Array.from(result);
  };
  try {
    onProgress({ message: "Computing the 3 × 3 matrix operation in WebAssembly" });
    await yieldTurn();
    const a = own(new Matrix(3, 3, m)), b = own(new Matrix(3, 3, n));
    let vectorResult, direction;
    if (axis) direction = own(new Matrix(3,1,axis));
    if (operation === "cross") {
      const columns = [];
      for (let col=0;col<3;col++) {
        const source = own(new Matrix(3,1,[m[col],m[3+col],m[6+col]]));
        const crossed = own(source.cross(direction));
        columns.push(values(crossed,3,1,"crossed column"));
      }
      vectorResult = new Matrix(3,3,Array.from({length:9},(_,i)=>columns[i%3][Math.floor(i/3)]));
    } else if (operation === "rotate") {
      const rotation = own(Matrix.rotationAxisAngle(direction,radians));
      checkArray(values(rotation,3,3,"rotation"),rotationReference(axis,radians),"rotation matrix");
      vectorResult = rotation.multiply(a);
    }
    const result = own(vectorResult || (operation === "multiply" ? a.multiply(b)
      : operation === "add" ? a.add(b)
      : operation === "subtract" ? a.subtract(b)
      : operation === "scale" ? a.scale(scalar) : a.transpose()));
    const resultValues = values(result, 3, 3, "result");
    const expected = reference(operation, m, n, scalar, axis, radians);
    checkArray(resultValues, expected, "matrix result");
    checkArray(values(a, 3, 3, "M"), m, "unchanged M");
    checkArray(values(b, 3, 3, "N"), n, "unchanged N");
    const determinantM = a.determinant(), determinantResult = result.determinant();
    const expectedMDet = determinantReference(m), expectedResultDet = determinantReference(expected);
    check(determinantM, expectedMDet.value, "determinant of M", expectedMDet.tolerance);
    check(determinantResult, expectedResultDet.value, "determinant of result", expectedResultDet.tolerance);

    let eigenM = null;
    if (m.every((value, index) => value === m[(index % 3) * 3 + Math.floor(index / 3)])) {
      const eigen = a.eigenSymmetric();
      if (eigen?.vectors) own(eigen.vectors);
      if (!eigen || (!Array.isArray(eigen.values) && !(eigen.values instanceof Float64Array))
        || eigen.values.length !== 3 || !Array.from(eigen.values).every(finite)) throw new Error("invalid input eigensystem");
      const q = values(eigen.vectors, 3, 3, "input eigenvectors"), lambdas = Array.from(eigen.values);
      const inputScale = Math.max(...m.map(Math.abs));
      const normalized = inputScale ? m.map(value => value / inputScale) : m;
      const norm = Math.hypot(...normalized);
      for (let col = 0; col < 3; col++) {
        if (col && lambdas[col] < lambdas[col-1]) throw new Error("input eigenvalues are not sorted");
        for (let row = 0; row < 3; row++) {
          let av = 0;
          for (let k = 0; k < 3; k++) av += normalized[row*3+k] * q[k*3+col];
          check(av, (inputScale ? lambdas[col] / inputScale : lambdas[col]) * q[row*3+col], "input eigenpair residual", 1e-9*norm);
          let dot = 0;
          for (let k = 0; k < 3; k++) dot += q[k*3+row] * q[k*3+col];
          check(dot, row === col ? 1 : 0, "input eigenvector orthogonality", 1e-9);
        }
      }
      eigenM = { values: lambdas, imag_values: [0,0,0], vectors: q, imag_vectors: Array(9).fill(0) };
    } else {
      const eigen = a.eigenGeneral();
      if(eigen?.vectorsReal) own(eigen.vectorsReal);
      if(eigen?.vectorsImag) own(eigen.vectorsImag);
      const vr=Array.from(eigen.valuesReal),vi=Array.from(eigen.valuesImag);
      if(vr.length!==3||vi.length!==3||!vr.every(finite)||!vi.every(finite))throw new Error("invalid complex input eigensystem");
      const qr=values(eigen.vectorsReal,3,3,"input eigenvectors real"),qi=values(eigen.vectorsImag,3,3,"input eigenvectors imaginary");
      const inputScale=Math.max(...m.map(Math.abs))||1,normalized=m.map(x=>x/inputScale),norm=Math.hypot(...normalized);
      for(let col=0;col<3;col++) {
        let length=0,error=0;
        if(col&&(vr[col]<vr[col-1]||vr[col]===vr[col-1]&&vi[col]<vi[col-1]))throw new Error("input complex eigenvalues are not sorted");
        for(let row=0;row<3;row++) {
          const re=qr[row*3+col],im=qi[row*3+col];length=Math.hypot(length,re,im);
          let ar=0,ai=0;
          for(let k=0;k<3;k++) {ar+=normalized[row*3+k]*qr[k*3+col];ai+=normalized[row*3+k]*qi[k*3+col];}
          error=Math.hypot(error,ar-(vr[col]/inputScale*re-vi[col]/inputScale*im),ai-(vr[col]/inputScale*im+vi[col]/inputScale*re));
        }
        check(length,1,"complex input eigenvector norm",1e-9);
        if(error>1e-9*norm)throw new Error("complex input eigenpair residual exceeds tolerance");
      }
      // Residuals alone permit duplicating one valid pair. Check the complete
      // cubic spectrum independently through its characteristic coefficients.
      const roots=vr.map((re,i)=>[re/inputScale,vi[i]/inputScale]);
      const product=([a,b],[c,d])=>[a*c-b*d,a*d+b*c];
      const sum=items=>items.reduce(([r,i],[x,y])=>[r+x,i+y],[0,0]);
      const coefficients=[sum(roots),sum([product(roots[0],roots[1]),product(roots[0],roots[2]),product(roots[1],roots[2])]),product(product(roots[0],roots[1]),roots[2])];
      const expectedCoefficients=[normalized[0]+normalized[4]+normalized[8],
        normalized[0]*normalized[4]-normalized[1]*normalized[3]+normalized[0]*normalized[8]-normalized[2]*normalized[6]+normalized[4]*normalized[8]-normalized[5]*normalized[7],
        determinantReference(normalized).value];
      coefficients.forEach(([re,im],i)=>{check(re,expectedCoefficients[i],"input complex spectrum",1e-8);check(im,0,"input complex spectrum imaginary coefficient",1e-8);});
      eigenM={values:vr,imag_values:vi,vectors:qr,imag_vectors:qi};
    }

    const points = grid(dimensions);
    onProgress({ message: `Verifying ${points.length} field vectors from actual WASM matrix products` });
    const pointMatrix = own(new Matrix(3, points.length, [0, 1, 2].flatMap(axis => points.map(point => point[axis]))));
    const inputProduct = own(a.multiply(pointMatrix)), outputProduct = own(result.multiply(pointMatrix));
    const inputFlat = values(inputProduct, 3, points.length, "M times grid");
    const outputFlat = values(outputProduct, 3, points.length, "result times grid");
    const inputVectors = [], outputVectors = [];
    for (let index = 0; index < points.length; index++) {
      const input = [0, 1, 2].map(axis => inputFlat[axis * points.length + index]);
      const output = [0, 1, 2].map(axis => outputFlat[axis * points.length + index]);
      checkArray(input, applyReference(m, points[index]), `M*x at point ${index}`);
      checkArray(output, applyReference(expected, points[index]), `result*x at point ${index}`);
      inputVectors.push(input);
      outputVectors.push(output);
    }

    const rotationFrames = [];
    let rotationExtent = 0;
    if (operation === "rotate") {
      // Each playback matrix is computed by WASM; these are angle samples, not CPU frames.
      onProgress({message:"Computing and verifying the rotation path in WebAssembly"});
      for (let frame=0;frame<=900;frame++) {
        const angle = radians*(frame/900);
        let rotation, product;
        try {
          rotation = Matrix.rotationAxisAngle(direction,angle);
          product = rotation.multiply(a);
          const matrix = values(product,3,3,"rotation frame");
          checkArray(matrix,multiplyReference(rotationReference(axis,angle),m),"rotation path");
          rotationFrames.push(matrix);
          for (let row=0;row<3;row++) rotationExtent=Math.max(rotationExtent,.55*(Math.abs(matrix[row*3])+Math.abs(matrix[row*3+1])+Math.abs(matrix[row*3+2])));
        } finally { product?.dispose(); rotation?.dispose(); }
        if (frame%100===0) await yieldTurn();
      }
    }

    const steps = [];
    const traceAvailable = operation === "multiply" && typeof TraceMatrix === "function";
    if (traceAvailable) {
      onProgress({ message: "Validating 27 arithmetic events captured from the instrumented C multiply loop" });
      await yieldTurn();
      const traceA = own(new TraceMatrix(3, 3, m)), traceB = own(new TraceMatrix(3, 3, n));
      if (typeof traceA.multiplyWithTrace !== "function") throw new Error("trace WASM build has no arithmetic capture API");
      const captured = traceA.multiplyWithTrace(traceB);
      if (captured?.matrix) own(captured.matrix);
      if (!captured || captured.source !== TRACE_SOURCE || !Array.isArray(captured.steps) || captured.steps.length !== 27) {
        throw new Error("trace must contain 27 captured multiply events and a source location");
      }
      const tracedValues = values(captured.matrix, 3, 3, "traced result");
      checkArray(tracedValues, resultValues, "traced versus ordinary WASM result");
      const partial = Array(9).fill(0), expectedSums = Array(9).fill(0);
      for (let index = 0; index < captured.steps.length; index++) {
        const step = captured.steps[index];
        const row = Math.floor(index / 9), col = Math.floor(index / 3) % 3, k = index % 3;
        if (!step || step.row !== row || step.col !== col || step.k !== k
          || !Number.isSafeInteger(step.source_line) || step.source_line < 1) {
          throw new Error("trace event order or source line is invalid");
        }
        const left = m[row * 3 + k], right = n[k * 3 + col], product = left * right;
        expectedSums[row * 3 + col] += product;
        check(step.left, left, `trace ${index} left operand`);
        check(step.right, right, `trace ${index} right operand`);
        check(step.product, product, `trace ${index} product`);
        check(step.sum, expectedSums[row * 3 + col], `trace ${index} accumulator`);
        partial[row * 3 + col] = step.sum;
        steps.push({ row, col, k, left: step.left, right: step.right, product: step.product,
          sum: step.sum, source_line: step.source_line, matrix: partial.slice(), call_path: CALL_PATH.slice() });
      }
      checkArray(partial, tracedValues, "final captured accumulator matrix");
      checkArray(values(traceA, 3, 3, "trace M"), m, "unchanged trace M");
      checkArray(values(traceB, 3, 3, "trace N"), n, "unchanged trace N");
    }

    return { operation, m, n, scalar, axis, radians, rotation_frames:rotationFrames, rotation_extent:rotationExtent,
      eigen_m: eigenM, result: resultValues, determinant_m: determinantM,
      determinant_result: determinantResult, points, input_vectors: inputVectors, output_vectors: outputVectors,
      steps, trace_available: traceAvailable, trace_source: traceAvailable ? TRACE_SOURCE : null,
      trace_reason: traceAvailable ? null : operation === "multiply"
        ? "The separate arithmetic trace build is unavailable."
        : "Arithmetic capture is available for matrix multiplication.",
      checks: { passed: verified, total: verified },
      convention: "Column vectors: (MN)x = M(Nx). Field arrows at x have head x + A*x.",
      animation: operation === "rotate" ? "901 verified WASM angle samples of R(axis,t*angle)*M, not CPU trace frames."
        : "Linear blend of operators, not a physical trajectory." };
  } finally {
    for (const matrix of owned.reverse()) matrix.dispose();
  }
}

async function fromEmbedded(bundle, trace = false) {
  const source = bundle?.[trace ? "trace_module_source" : "module_source"];
  const binary = bundle?.[trace ? "trace_wasm_base64" : "wasm_base64"];
  if (typeof bundle?.wrapper_source !== "string" || typeof source !== "string" || typeof binary !== "string") {
    throw new Error(trace ? "embedded arithmetic trace module is incomplete" : "embedded geometry module is incomplete");
  }
  const wrapperURL = URL.createObjectURL(new Blob([bundle.wrapper_source], { type: "text/javascript" }));
  const moduleURL = URL.createObjectURL(new Blob([source], { type: "text/javascript" }));
  try {
    const { createMatrixAPI } = await import(wrapperURL);
    const bytes = Uint8Array.from(atob(binary), character => character.charCodeAt(0));
    return await createMatrixAPI({ moduleUrl: moduleURL, wasmBinary: bytes, locateFile: () => "embedded.wasm" });
  } finally {
    URL.revokeObjectURL(wrapperURL);
    URL.revokeObjectURL(moduleURL);
  }
}

if (typeof WorkerGlobalScope !== "undefined" && globalThis instanceof WorkerGlobalScope) {
  let started = false;
  const handleGeometry = async ({ data }) => {
    if (started) return;
    started = true;
    let sharedChecks;
    try {
      if (data?.type !== "geometry") throw new Error("unsupported geometry request");
      if (typeof data.bundle?.worker_source !== "string") throw new Error("shared correctness worker source is unavailable");
      const progress = event => globalThis.postMessage({ type: "progress", message: event.message });
      progress({ message: "Loading and checking the embedded WASM library" });
      const checksURL = URL.createObjectURL(new Blob([data.bundle.worker_source], { type: "text/javascript" }));
      let runChecks;
      try {
        ({ runChecks } = await import(checksURL));
      } finally {
        // The shared worker installs its own message handler when imported in a
        // worker. Restore this one before any further asynchronous initialization.
        globalThis.onmessage = handleGeometry;
        URL.revokeObjectURL(checksURL);
      }
      const { Matrix } = await fromEmbedded(data.bundle);
      sharedChecks = await runChecks(Matrix, data.bundle.fixtures, progress);
      if (sharedChecks.passed !== sharedChecks.total) throw new Error("shared correctness checks failed; geometry was not computed");
      let TraceMatrix;
      if (data.config?.operation === "multiply" && data.bundle.trace_module_source && data.bundle.trace_wasm_base64) {
        ({ Matrix: TraceMatrix } = await fromEmbedded(data.bundle, true));
      }
      const result = await computeGeometry(Matrix, data.config, { TraceMatrix, onProgress: progress });
      result.shared_checks = sharedChecks;
      globalThis.postMessage({ type: "geometry", result });
    } catch (error) {
      globalThis.postMessage({ type: "error", message: error instanceof Error ? error.message : String(error),
        ...(sharedChecks ? { checks: sharedChecks } : {}) });
    }
  };
  globalThis.onmessage = handleGeometry;
}
