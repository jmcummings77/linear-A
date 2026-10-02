# F# float64 matrices

The .NET 10 `LinearA.FSharp` library exposes `LinearA.Matrix`, an owned, flat
row-major `double array` with zero-based indexed access. Arithmetic methods return
independent matrices. `FromArray`, `Copy`, `ToArray`, `GetRow`, and `GetColumn`
copy their storage; only indexer assignments change an existing matrix.

```fsharp
open LinearA

let a = Matrix.FromArray(2, 3, [| 1.; 2.; 3.; 4.; 5.; 6. |])
let gram = a.Multiply(a.Transpose())
let scaled = gram.Scale(0.5)
printfn "%A" (scaled.ToArray()) // [| 7.; 16.; 16.; 38.5 |]
```

Construct a zero matrix with `Matrix(rows, cols)` and an identity matrix with
`Matrix.Identity(size)`. The API provides `Rows`, `Cols`, `matrix[row, column]`,
`Add`, `Subtract`, `Scale`, `Transpose`, `Multiply`, `Cross`, `Trace`, `Determinant`,
`EigenSymmetric`, `EigenGeneral`, `IsUpperTriangular`, `IsLowerTriangular`, and
`IsTriangular`.

Invalid dimensions, indices, and incompatible operands raise argument errors;
trace and determinant require square matrices. Empty dimensions are supported:
`Matrix(0, 0)` has trace zero and determinant one. Triangular checks use exact
zeros. `Determinant()` defaults to Auto, which shortcuts triangular matrices;
`Determinant("lu")` selects partial-pivot elimination and `Determinant("cholesky")`
requires exact symmetry and positive computed pivots. Both factorizations use a
copied work array and cubic arithmetic work, with scaled determinant products.
Determinants reject nonfinite inputs/results. Ordinary float64 rounding,
cancellation, underflow, and intermediate overflow still apply. Library arithmetic follows .NET floating-point
behavior; the comparison runner validates finite input and rejects nonfinite
output before JSON serialization. It does not estimate numerical conditioning.

`matrix.EigenSymmetric()` returns `(eigenvalues, eigenvectors)` in independent
storage. Eigenvalues are ascending; matching unit eigenvectors are the **columns**
of the returned matrix. For example, `Matrix.FromArray(2, 2, [| 2.; 1.; 1.; 2. |])`
has eigenvalues `[| 1.; 3. |]`. Only finite, exactly symmetric square matrices are
accepted. Empty input returns empty values and a 0×0 matrix.

The cyclic Jacobi solver accepts `tolerance` (default `1e-12`, finite and strictly
between zero and one) and `maxSweeps` (default `50`, positive). The tolerance is
relative to the input's Frobenius norm, so tiny eigenvalues alongside huge entries
may have poor relative accuracy. Nonconvergence and nonfinite eigenvalues raise
errors. Repeated eigenvalues can have any orthonormal basis within their eigenspace;
individual vector signs carry no mathematical significance. The symmetric API remains specialized for real orthonormal eigenvectors.

`matrix.EigenGeneral()` accepts any finite real square matrix and returns
`(System.Numerics.Complex array, System.Numerics.Complex[,])`: eigenvalues sorted
by real part then imaginary part, with corresponding unit right eigenvectors in
columns. A quarter-turn matrix `[| 0.; -1.; 1.; 0. |]` has values `-i` and `+i`.
The source is preserved; empty input returns empty arrays. Vectors for repeated
or defective eigenvalues need not be independent or orthogonal, so the result
does not promise an invertible eigenvector matrix.

The general solver uses power-of-two similarity balancing, Hessenberg reduction,
and real double-shift QR, adapted from the public-domain
[NIST JAMA algorithm](https://math.nist.gov/javanumerics/jama/). It uses machine
epsilon for deflation; `maxIterations` defaults to `1000` per unresolved root
and accepts `1..100000`. Nonconvergence and nonfinite results raise errors.
Scaling protects intermediate range and preserves exact triangular eigenvalues,
but this is a floating-point approximation without condition estimates; sensitive
eigenvalues and eigenvectors can lose accuracy. Complex column phases are arbitrary
mathematically; the solver chooses a deterministic phase for its output.

`a.Cross(b)` computes the right-handed 3D vector cross product. Both inputs must
be 3×1 or 1×3; they may use different orientations. The independent result keeps
`a`'s shape. Finite inputs are required and nonfinite products fail explicitly.

`Matrix.Rotation2D(radians)`, `RotationX`, `RotationY`, `RotationZ`, and
`Matrix.RotationAxisAngle(axis, radians)` construct independent rotation matrices.
These are active, right-handed rotations of column vectors: apply `R.Multiply(v)`.
For a row vector, apply `v.Multiply(R.Transpose())`. The axis accepts either
three-component vector shape, must be finite and nonzero, and is normalized without
mutating it. All angles are finite radians. Positive Z rotation takes X toward Y;
a positive 2D angle is counterclockwise. To rotate by A then B, compose `B.Multiply(A)`.

From the repository root:

```sh
dotnet fsi ports/fsharp/tests.fsx
dotnet build ports/fsharp/Runner.fsproj --configuration Release
printf '1 2 3 4' | dotnet ports/fsharp/bin/Release/net10.0/Runner.dll check determinant 2 2
dotnet ports/fsharp/bin/Release/net10.0/Runner.dll bench multiply 32 100 42
```

For use from another .NET project, reference `LinearA.FSharp.fsproj`. The runner
implements the [shared protocol](../../benchmarks/PROTOCOL.md), uses .NET's
monotonic `Stopwatch`, and includes allocation and naturally scheduled GC in the
timed loop. It uses handwritten loops and no BLAS or third-party numeric library.

## Linear systems

See [the shared solver guide](../SOLVING.md) for this implementation’s reusable
LU, Cholesky, and column-pivoted QR APIs, multiple right-hand sides, least squares,
condition diagnostics, ownership, and numerical limits.
