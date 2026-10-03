/** Owned float64 matrices backed by the C implementation compiled to WebAssembly. */
const MAX_DIMENSION = 0xffffffff;
const MAX_ELEMENTS = Math.floor(MAX_DIMENSION / 8);

function dimension(value) {
  if (!Number.isSafeInteger(value) || value < 0 || value > MAX_DIMENSION) {
    throw new RangeError("matrix dimensions must be nonnegative wasm32 integers");
  }
  return value;
}

function finite(value) {
  if (typeof value !== "number" || !Number.isFinite(value)) throw new RangeError("expected a finite number");
  return value;
}

/**
 * Initialize once, then use the returned { Matrix, Factorization } classes. moduleUrl points to
 * the generated matrix.mjs; locateFile can locate its matrix.wasm in a browser.
 * wasmBinary optionally supplies Uint8Array bytes for an embedded/offline build.
 * Every Matrix and Factorization owns its allocation and must be disposed explicitly.
 */
export async function createMatrixAPI({
  moduleUrl = new URL("../../.build/wasm/matrix.mjs", import.meta.url),
  locateFile,
  wasmBinary,
} = {}) {
  const { default: createModule } = await import(String(moduleUrl));
  const runtimeOptions = {};
  if (locateFile !== undefined) runtimeOptions.locateFile = locateFile;
  if (wasmBinary !== undefined) runtimeOptions.wasmBinary = wasmBinary;
  const runtime = await createModule(runtimeOptions);
  const ownedHandle = Symbol("owned matrix handle");

  function checkedStatus() {
    if (runtime._wm_last_error() !== 0) {
      throw new RangeError(runtime.UTF8ToString(runtime._wm_error_message()));
    }
  }

  function checkedHandle(pointer) {
    checkedStatus();
    if (!pointer) throw new RangeError("WebAssembly matrix allocation failed");
    return pointer >>> 0;
  }

  let solveFactor, matrixPointer;
  const ownedFactor = Symbol("owned factor handle");

  class Matrix {
    static { matrixPointer = value => value.#alive(); }
    #pointer = 0;
    #rows = 0;
    #cols = 0;

    constructor(rows = 0, cols = 0, values, token) {
      if (token === ownedHandle) {
        this.#pointer = checkedHandle(rows);
        this.#rows = runtime._wm_rows(this.#pointer) >>> 0;
        this.#cols = runtime._wm_cols(this.#pointer) >>> 0;
        return;
      }
      this.#rows = dimension(rows);
      this.#cols = dimension(cols);
      const count = rows * cols;
      if (!Number.isSafeInteger(count) || count > MAX_ELEMENTS) throw new RangeError("matrix dimensions overflow wasm32 storage");
      let initial;
      if (values !== undefined) {
        if (!Array.isArray(values) && !(values instanceof Float64Array)) throw new TypeError("values must be an Array or Float64Array");
        if (values.length !== count) throw new RangeError("value count does not match matrix dimensions");
        initial = new Float64Array(count);
        // Read each value once, so array getters or concurrent shared-memory
        // writers cannot change an already-validated value before the copy.
        for (let i = 0; i < count; i++) initial[i] = finite(values[i]);
      }
      this.#pointer = checkedHandle(runtime._wm_create(rows, cols));
      try {
        if (initial !== undefined) runtime.HEAPF64.set(initial, (runtime._wm_data(this.#pointer) >>> 0) / 8);
      } catch (error) {
        this.dispose();
        throw error;
      }
    }

    static #fromHandle(pointer) { return new Matrix(pointer, 0, undefined, ownedHandle); }
    static identity(size) { return Matrix.#fromHandle(runtime._wm_identity(dimension(size))); }
    /** Right-handed active rotations of column vectors; angles are radians. */
    static rotation2D(radians) { return Matrix.#fromHandle(runtime._wm_rotation_2d(finite(radians))); }
    static rotationX(radians) { return Matrix.#fromHandle(runtime._wm_rotation_x(finite(radians))); }
    static rotationY(radians) { return Matrix.#fromHandle(runtime._wm_rotation_y(finite(radians))); }
    static rotationZ(radians) { return Matrix.#fromHandle(runtime._wm_rotation_z(finite(radians))); }
    static rotationAxisAngle(axis, radians) {
      if (!(axis instanceof Matrix)) throw new TypeError("axis must be a Matrix from the same WebAssembly instance");
      return Matrix.#fromHandle(runtime._wm_rotation_axis_angle(axis.#alive(), finite(radians)));
    }

    #alive() {
      if (!this.#pointer) throw new Error("matrix has been disposed");
      return this.#pointer;
    }

    #other(other) {
      if (!(other instanceof Matrix)) throw new TypeError("operand must be a Matrix from the same WebAssembly instance");
      return other.#alive();
    }

    #index(row, col) {
      this.#alive();
      if (!Number.isSafeInteger(row) || !Number.isSafeInteger(col) || row < 0 || col < 0 || row >= this.#rows || col >= this.#cols) {
        throw new RangeError("matrix index out of bounds");
      }
      return row * this.#cols + col;
    }

    get rows() { this.#alive(); return this.#rows; }
    get cols() { this.#alive(); return this.#cols; }
    get disposed() { return this.#pointer === 0; }

    get(row, col) {
      const index = this.#index(row, col);
      return runtime.HEAPF64[(runtime._wm_data(this.#pointer) >>> 0) / 8 + index];
    }

    set(row, col, value) {
      this.#index(row, col);
      runtime._wm_set(this.#pointer, row, col, finite(value));
      checkedStatus();
    }

    /** Return an independent JS copy; no heap views survive WASM memory growth. */
    toArray() {
      this.#alive();
      const start = (runtime._wm_data(this.#pointer) >>> 0) / 8;
      return runtime.HEAPF64.slice(start, start + this.#rows * this.#cols);
    }

    copy() { return Matrix.#fromHandle(runtime._wm_copy(this.#alive())); }
    add(other) { return Matrix.#fromHandle(runtime._wm_add(this.#alive(), this.#other(other))); }
    subtract(other) { return Matrix.#fromHandle(runtime._wm_subtract(this.#alive(), this.#other(other))); }
    scale(scalar) { return Matrix.#fromHandle(runtime._wm_scale(this.#alive(), finite(scalar))); }
    transpose() { return Matrix.#fromHandle(runtime._wm_transpose(this.#alive())); }
    multiply(other) { return Matrix.#fromHandle(runtime._wm_multiply(this.#alive(), this.#other(other))); }
    /** Cross two 3D row/column vectors; the result retains this matrix's shape. */
    cross(other) { return Matrix.#fromHandle(runtime._wm_cross(this.#alive(), this.#other(other))); }

    /** Economy SVD; values are copied, and u/vt must each be disposed. */
    svd(tolerance = 1e-12, maxSweeps = 100) {
      this.#alive();
      if (!Number.isFinite(tolerance) || tolerance <= 0 || tolerance >= 1 || !Number.isInteger(maxSweeps) || maxSweeps < 1 || maxSweeps > 10000) throw new RangeError("invalid SVD options");
      let u, values, vt;
      try {
        u = new Matrix(); values = new Matrix(); vt = new Matrix();
        runtime._wm_svd(this.#pointer,tolerance,maxSweeps,u.#pointer,values.#pointer,vt.#pointer);
        checkedStatus();
        for (const item of [u,values,vt]) {
          item.#rows=runtime._wm_rows(item.#pointer)>>>0; item.#cols=runtime._wm_cols(item.#pointer)>>>0;
        }
        const copied=values.toArray();
        return {u,values:copied,vt};
      } catch(error) { u?.dispose(); vt?.dispose(); throw error; }
      finally { values?.dispose(); }
    }

    /** Truncated SVD inverse; dispose the returned matrix. */
    solveRidge(rhs,lambda) {
      if(!Number.isFinite(lambda)||lambda<0)throw new RangeError('lambda must be finite and nonnegative');
      return Matrix.#fromHandle(runtime._wm_solve_ridge(this.#alive(),this.#other(rhs),lambda));
    }
    pseudoinverse(cutoff=Math.max(this.rows,this.cols)*Number.EPSILON) {
      this.#checkCutoff(cutoff);
      return Matrix.#fromHandle(runtime._wm_pseudoinverse(this.#alive(),cutoff));
    }
    solveMinimumNorm(rhs,cutoff=Math.max(this.rows,this.cols)*Number.EPSILON) {
      this.#checkCutoff(cutoff);
      return Matrix.#fromHandle(runtime._wm_solve_minimum_norm(this.#alive(),this.#other(rhs),cutoff));
    }
    #checkCutoff(cutoff) {
      if(typeof cutoff!=="number" || !Number.isFinite(cutoff)||cutoff<0||cutoff>1)throw new RangeError("invalid relative cutoff");
    }
    spectralDiagnostics(cutoff=Math.max(this.rows,this.cols)*Number.EPSILON) {
      this.#checkCutoff(cutoff);
      const result=Matrix.#fromHandle(runtime._wm_spectral_diagnostics(this.#alive(),cutoff));
      try {const [rank,reciprocalCondition,retainedReciprocalCondition]=result.toArray();return {rank,reciprocalCondition,retainedReciprocalCondition};}
      finally {result.dispose();}
    }
    factorLU() { return new Factorization(runtime._wm_factorize(this.#alive(),1),ownedFactor); }
    factorCholesky() { return new Factorization(runtime._wm_factorize(this.#alive(),2),ownedFactor); }
    factorQR() { return new Factorization(runtime._wm_factorize(this.#alive(),3),ownedFactor); }
    solve(rhs) { const factor=this.factorLU();try{return factor.solve(rhs);}finally{factor.dispose();} }
    leastSquares(rhs) { const factor=this.factorQR();try{return factor.solve(rhs);}finally{factor.dispose();} }
    static {
      solveFactor=(pointer,rhs)=>{
        if(!(rhs instanceof Matrix))throw new TypeError("right-hand side must be a Matrix from the same WebAssembly instance");
        return Matrix.#fromHandle(runtime._wm_factor_solve(pointer,rhs.#alive()));
      };
    }

    static {
      // Only the separate trace build exposes this operation. The events come
      // from the C inner loop; JS copies them before another call can reuse it.
      if (typeof runtime._wm_multiply_trace === "function") {
        Object.defineProperty(this.prototype, "multiplyWithTrace", { value(other) {
          const left = this.#alive();
          const right = this.#other(other);
          if ([this.#rows, this.#cols, other.#rows, other.#cols].some(value => value > 8)) {
            throw new RangeError("traced multiplication supports dimensions from 0 to 8");
          }
          const matrix = Matrix.#fromHandle(runtime._wm_multiply_trace(left, right));
          try {
            const count = runtime._wm_trace_count() >>> 0;
            const stride = runtime._wm_trace_stride() >>> 0;
            if (stride !== 8 || count > 512 || count !== this.#rows * other.#cols * this.#cols) {
              throw new Error("invalid WebAssembly multiplication trace");
            }
            const start = (runtime._wm_trace_data() >>> 0) / 8;
            const heap = runtime.HEAPF64;
            const steps = Array.from({ length: count }, (_, index) => {
              const offset = start + index * stride;
              return {
                row: heap[offset], col: heap[offset + 1], k: heap[offset + 2],
                left: heap[offset + 3], right: heap[offset + 4], product: heap[offset + 5],
                sum: heap[offset + 6], source_line: heap[offset + 7],
              };
            });
            return { matrix, steps, source: "ports/c/matrix.c" };
          } catch (error) {
            matrix.dispose();
            throw error;
          }
        } });
      }
    }

    row(index) {
      this.#alive();
      if (!Number.isSafeInteger(index) || index < 0 || index >= this.#rows) throw new RangeError("row index out of bounds");
      const result = Matrix.#fromHandle(runtime._wm_row(this.#pointer, index));
      try { return result.toArray(); } finally { result.dispose(); }
    }

    column(index) {
      this.#alive();
      if (!Number.isSafeInteger(index) || index < 0 || index >= this.#cols) throw new RangeError("column index out of bounds");
      const result = Matrix.#fromHandle(runtime._wm_column(this.#pointer, index));
      try { return result.toArray(); } finally { result.dispose(); }
    }

    trace() {
      const result = runtime._wm_trace(this.#alive());
      checkedStatus();
      return result;
    }

    /** Auto uses small/triangular shortcuts; Cholesky requires exact symmetry and positive definiteness. */
    determinant(algorithm = "auto") {
      const pointer = this.#alive();
      const selected = ["auto", "lu", "cholesky"].indexOf(algorithm);
      if (selected < 0) throw new RangeError("unknown determinant algorithm; expected auto, lu, or cholesky");
      const result = runtime._wm_determinant_algorithm(pointer, selected);
      checkedStatus();
      return result;
    }

    triangular() {
      const result = runtime._wm_triangular(this.#alive());
      checkedStatus();
      return { upper: Boolean(result & 1), lower: Boolean(result & 2) };
    }

    /** Real symmetric eigendecomposition; values ascend and vectors are columns.
     * The returned values are copied; dispose the independently owned vectors. */
    eigenSymmetric({ tolerance = 1e-12, maxSweeps = 50 } = {}) {
      const pointer = this.#alive();
      if (typeof tolerance !== "number" || !Number.isFinite(tolerance) || tolerance <= 0 || tolerance >= 1) {
        throw new RangeError("eigen tolerance must be finite and between zero and one");
      }
      if (!Number.isSafeInteger(maxSweeps) || maxSweeps < 1 || maxSweeps > 10000) {
        throw new RangeError("maxSweeps must be an integer from 1 to 10000");
      }
      const packed = Matrix.#fromHandle(runtime._wm_eigen_symmetric(pointer, tolerance, maxSweeps));
      try {
        const n = this.#rows, data = packed.toArray();
        const values = new Float64Array(n), vectors = new Float64Array(n * n);
        for (let row = 0; row < n; row++) {
          values[row] = data[row * (n + 1)];
          vectors.set(data.subarray(row * (n + 1) + 1, (row + 1) * (n + 1)), row * n);
        }
        return { values, vectors: new Matrix(n, n, vectors) };
      } finally { packed.dispose(); }
    }

    /** General complex eigenpairs of a real square matrix. Values are copied;
     * dispose both independently owned real and imaginary vector matrices. */
    eigenGeneral({ maxIterations = 1000 } = {}) {
      const pointer = this.#alive();
      if (!Number.isSafeInteger(maxIterations) || maxIterations < 1 || maxIterations > 100000) {
        throw new RangeError("maxIterations must be an integer from 1 to 100000");
      }
      const packed = Matrix.#fromHandle(runtime._wm_eigen_general(pointer, maxIterations));
      let vectorsReal;
      try {
        const n = this.#rows, data = packed.toArray(), stride = 2 * n + 2;
        const valuesReal = new Float64Array(n), valuesImag = new Float64Array(n);
        const real = new Float64Array(n * n), imag = new Float64Array(n * n);
        for (let row = 0; row < n; row++) {
          valuesReal[row] = data[row * stride]; valuesImag[row] = data[row * stride + 1];
          real.set(data.subarray(row * stride + 2, row * stride + 2 + n), row * n);
          imag.set(data.subarray(row * stride + 2 + n, (row + 1) * stride), row * n);
        }
        vectorsReal = new Matrix(n, n, real);
        const vectorsImag = new Matrix(n, n, imag);
        return { valuesReal, valuesImag, vectorsReal, vectorsImag };
      } catch (error) { vectorsReal?.dispose(); throw error; }
      finally { packed.dispose(); }
    }

    /** Benchmark consumption: sum first, middle, and last values inside WASM. */
    checksum() {
      const result = runtime._wm_checksum(this.#alive());
      checkedStatus();
      return result;
    }

    dispose() {
      if (this.#pointer) runtime._wm_destroy(this.#pointer);
      this.#pointer = 0;
    }
  }

  class Factorization {
    #pointer=0;
    constructor(pointer,token){
      if(token!==ownedFactor)throw new TypeError("Create factors through a Matrix factor method");
      this.#pointer=checkedHandle(pointer);
    }
    #alive(){if(!this.#pointer)throw new Error("factorization has been disposed");return this.#pointer;}
    get disposed(){return this.#pointer===0;}
    solve(rhs){return solveFactor(this.#alive(),rhs);}
    reciprocalCondition(){const value=runtime._wm_factor_rcond(this.#alive());checkedStatus();return value;}
    dispose(){if(this.#pointer)runtime._wm_factor_destroy(this.#pointer);this.#pointer=0;}
  }

  let csrPointer,iluPointer;
  class ILU0 {
    #pointer;#size;#nnz;
    static {iluPointer=f=>f.#alive();}
    constructor(a){this.#size=a.rows;this.#nnz=a.nnz;this.#pointer=checkedHandle(runtime._wm_ilu0_create(csrPointer(a)));}
    get size(){return this.#size;}get nnz(){return this.#nnz;}
    #alive(){if(!this.#pointer)throw new Error('ILU0 has been disposed');return this.#pointer;}
    dispose(){if(this.#pointer){runtime._wm_ilu0_destroy(this.#pointer);this.#pointer=0;}}
    apply(values){let b,out;try{b=new Matrix(values.length,1,values);out=new Matrix(runtime._wm_ilu0_apply(this.#alive(),matrixPointer(b)),0,undefined,ownedHandle);return out.toArray();}finally{b?.dispose();out?.dispose();}}
  }
  class CSRMatrix {
    static {csrPointer=a=>a.#alive();}
    #pointer=0; #rows; #cols; #offsets; #indices; #values;
    constructor(rows,cols,offsets,indices,values){
      this.#rows=dimension(rows);this.#cols=dimension(cols);
      this.#offsets=Array.from(offsets);this.#indices=Array.from(indices);this.#values=Array.from(values);
      let rp,ci,v;
      try{
        rp=new Matrix(this.#offsets.length,1,this.#offsets);ci=new Matrix(this.#indices.length,1,this.#indices);v=new Matrix(this.#values.length,1,this.#values);
        this.#pointer=checkedHandle(runtime._wm_csr_create(rows,cols,matrixPointer(rp),matrixPointer(ci),matrixPointer(v)));
      }finally{rp?.dispose();ci?.dispose();v?.dispose();}
    }
    get rows(){return this.#rows;}get cols(){return this.#cols;}get nnz(){return this.#values.length;}
    get rowOffsets(){return this.#offsets.slice();}get columnIndices(){return this.#indices.slice();}get values(){return this.#values.slice();}
    #alive(){if(!this.#pointer)throw new Error('CSR matrix has been disposed');return this.#pointer;}
    dispose(){if(this.#pointer){runtime._wm_csr_destroy(this.#pointer);this.#pointer=0;}}
    reverseCuthillMcKee(){let out;try{out=new Matrix(runtime._wm_csr_rcm(this.#alive()),0,undefined,ownedHandle);return Array.from(out.toArray());}finally{out?.dispose();}}
    static permuteVector(order,x,inverse=false){
      if(typeof inverse!=='boolean')throw new TypeError('inverse must be boolean');
      let p,b,out;try{p=new Matrix(order.length,1,order);b=new Matrix(x.length,1,x);out=new Matrix(runtime._wm_permute_vector(matrixPointer(p),matrixPointer(b),+inverse),0,undefined,ownedHandle);return Array.from(out.toArray());}finally{p?.dispose();b?.dispose();out?.dispose();}
    }
    permuteSymmetric(order){
      let p,out;try{p=new Matrix(order.length,1,order);out=new Matrix(runtime._wm_csr_permute(this.#alive(),matrixPointer(p)),0,undefined,ownedHandle);const data=out.toArray(),n=this.rows,m=this.nnz;return new CSRMatrix(n,n,data.slice(0,n+1),data.slice(n+1,n+1+m),data.slice(n+1+m));}finally{p?.dispose();out?.dispose();}
    }
    static fromDense(a){
      const rp=[0],ci=[],v=[];
      for(let i=0;i<a.rows;i++){for(let j=0;j<a.cols;j++){const value=a.get(i,j);if(value!==0){ci.push(j);v.push(value);}}rp.push(v.length);}
      return new CSRMatrix(a.rows,a.cols,rp,ci,v);
    }
    matvec(values){
      let x,out;
      try{x=new Matrix(values.length,1,values);out=new Matrix(runtime._wm_csr_matvec(this.#alive(),matrixPointer(x)),0,undefined,ownedHandle);return out.toArray();}
      finally{x?.dispose();out?.dispose();}
    }
    conjugateGradient(values,{rtol=1e-10,atol=0,maxIterations=1000,jacobi=false,capture=false}={}){
      if(!Number.isInteger(maxIterations)||maxIterations<0||maxIterations>100000||typeof jacobi!=='boolean'||typeof capture!=='boolean')throw new RangeError('invalid CG options');
      finite(rtol);finite(atol);
      let b,handle;
      try{
        b=new Matrix(values.length,1,values);
        handle=checkedHandle(runtime._wm_csr_cg(this.#alive(),matrixPointer(b),rtol,atol,maxIterations,+jacobi,+capture));
        const iterations=runtime._wm_cg_iterations(handle),reason=['converged','iteration_limit','breakdown','nonfinite'][runtime._wm_cg_reason(handle)];
        const copy=(field,count)=>{const offset=runtime._wm_cg_data(handle,field)>>>3;return Array.from(runtime.HEAPF64.subarray(offset,offset+count));};
        const x=copy(0,this.rows),residuals=copy(1,iterations+1),flat=capture?copy(2,(iterations+1)*this.rows):[];
        const iterates=capture?Array.from({length:iterations+1},(_,i)=>flat.slice(i*this.rows,(i+1)*this.rows)):[];
        return {x,residuals,iterates,iterations,reason,converged:reason==='converged'};
      }finally{if(handle)runtime._wm_cg_destroy(handle);b?.dispose();}
    }
    gmres(values,{restart=30,rtol=1e-10,atol=0,maxIterations=1000,jacobi=false,capture=false,preconditioner}={}){
      if(!Number.isInteger(restart)||restart<1||restart>1024||!Number.isInteger(maxIterations)||maxIterations<0||maxIterations>100000||typeof jacobi!=='boolean'||typeof capture!=='boolean')throw new RangeError('invalid GMRES options');
      finite(rtol);finite(atol);let b,handle;
      try{
        b=new Matrix(values.length,1,values);handle=checkedHandle(runtime._wm_csr_gmres_preconditioned(this.#alive(),matrixPointer(b),restart,rtol,atol,maxIterations,+jacobi,+capture,preconditioner?iluPointer(preconditioner):0));
        const iterations=runtime._wm_gmres_iterations(handle),reason=['converged','iteration_limit','breakdown','nonfinite','stagnation'][runtime._wm_gmres_reason(handle)];
        const copy=(field,count)=>{const offset=runtime._wm_gmres_data(handle,field)>>>3;return Array.from(runtime.HEAPF64.subarray(offset,offset+count));};
        const x=copy(0,this.rows),residuals=copy(1,iterations+1),estimatedResiduals=copy(2,iterations+1),flat=capture?copy(3,(iterations+1)*this.rows):[];
        const iterates=capture?Array.from({length:iterations+1},(_,i)=>flat.slice(i*this.rows,(i+1)*this.rows)):[];
        const offset=runtime._wm_gmres_restarts(handle)>>>2,restarts=Array.from(new Uint32Array(runtime.HEAPF64.buffer,offset*4,runtime._wm_gmres_restart_count(handle)));
        return {x,residuals,estimatedResiduals,iterates,restarts,iterations,reason,converged:reason==='converged'};
      }finally{if(handle)runtime._wm_gmres_destroy(handle);b?.dispose();}
    }

  }

  if (typeof Symbol.dispose === "symbol") {
    Object.defineProperty(Matrix.prototype, Symbol.dispose, { value: Matrix.prototype.dispose });
    Object.defineProperty(Factorization.prototype, Symbol.dispose, { value: Factorization.prototype.dispose });
  }
  return { Matrix, Factorization, CSRMatrix, ILU0 };
}
