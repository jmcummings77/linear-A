# C# implementation

The C# implementation provides `Matrix<T>`, `DoubleMatrix`, and `IntegerMatrix`.
It shares the [float64 contract](../ports/README.md) and
[comparison protocol](../benchmarks/PROTOCOL.md) with the other implementations,
and also supports generic numeric types and legacy .NET APIs. Public indices
are zero-based. See the [project overview](../README.md) to choose a language.

Reference [linear-A.csproj](linear-A/linear-A.csproj) from a .NET application and
import the `linear_A` namespace. The source project has no external runtime
dependencies.

## Getting started

Install the [.NET 10 SDK](https://dotnet.microsoft.com/download/dotnet/10.0).
The SDK version is recorded in `global.json`. From the repository root:

```sh
dotnet restore --locked-mode
dotnet build --configuration Release --no-restore
dotnet test --configuration Release --no-build --no-restore
dotnet restore tests/NetStandardSmoke/NetStandardSmoke.csproj --locked-mode
dotnet run --project tests/NetStandardSmoke/NetStandardSmoke.csproj --configuration Release --no-restore
```

The library targets both .NET 10 and .NET Standard 2.1. `IntegerMatrix` and
`DoubleMatrix` remain available on both targets; `Matrix<T>` requires
the .NET 10 target. Tests run on .NET 10. A separate smoke runner explicitly loads
the .NET Standard 2.1 assembly and exercises its legacy API at runtime; CI runs
both checks on Linux, Windows, and macOS. This verifies that target on .NET 10,
not every runtime that implements .NET Standard. The library has no external
runtime dependencies.

## Generic matrices

[Generic math](https://learn.microsoft.com/dotnet/standard/generics/math) supplies
the numeric operators that were missing when this project started. `Matrix<T>`
uses `INumber<T>` to support types such as `int`, `double`, `decimal`, and
`BigInteger` without duplicating each matrix operation.

```csharp
using linear_A;

var a = new Matrix<double>(new double[,]
{
    { 1, 2, 3 },
    { 4, 5, 6 }
});

var gram = a * a.Transpose();
// [14, 32]
// [32, 77]

Console.WriteLine(gram.GetDeterminant()); // 54

var identity = new Matrix<double>(2, initializeAsIdentity: true);
var sum = gram + identity; // Creates a new matrix.
sum.Scale(0.5);            // Changes sum in place.
```

The array and copy constructors copy their input. `GetRow` and `GetColumn`
return copies as well. Operators, `DotProduct`, and `Transpose` create new
matrices; indexer assignments, `Scale`, `TryAddMatrix`, and `TrySubtractMatrix`
change the receiving matrix. Addition and subtraction require matching shapes;
multiplication requires the left column count to equal the right row count.

## Legacy matrices and conversions

`DoubleMatrix` and `IntegerMatrix` remain supported on both targets. Use
`DoubleMatrix.FromArray(double[,])` or `IntegerMatrix.FromArray(int[,])` to create
a matrix from a rectangular array, and `ToArray()` to export a rectangular copy.
Both boundaries copy their storage: later edits to either array or matrix do not
affect the other. The factories accept empty rectangular shapes such as 2×0;
the original dimension constructors still require positive dimensions.
Their storage and core operations share one internal implementation; determinant
elimination is also shared with `Matrix<T>`.
The old copy constructors still treat null as an empty 0×0 matrix; the new
factories reject null.

On .NET 10, `ToMatrix()` and `FromMatrix(...)` copy between the legacy type and
its matching generic type:

```csharp
var legacy = DoubleMatrix.FromArray(new double[,] { { 1, 2 }, { 3, 4 } });
Matrix<double> generic = legacy.ToMatrix();
DoubleMatrix copiedBack = DoubleMatrix.FromMatrix(generic);
double[,] snapshot = copiedBack.ToArray();
```

`IntegerMatrix` provides the same conversions with `Matrix<int>`. These are
independent copies, and they preserve empty rectangular shapes. Generic
conversion methods are absent from the .NET Standard target.

Both legacy classes implement typed `IEnumerable<double>` or `IEnumerable<int>`
in row-major order, so they work with LINQ. Their original nongeneric
`GetEnumerator()` remains available. Converting to a generic matrix does not
carry over legacy unchecked integer arithmetic: subsequent `Matrix<int>`
operations use checked arithmetic.

## Math and numeric behavior

- `Matrix<T>` uses checked arithmetic where the numeric type supports it.
  Bounded integer operations can throw `OverflowException`; use
  `Matrix<BigInteger>` for unbounded integer arithmetic. Floating-point
  calculations still round and can overflow or underflow. Intermediate values
  can overflow even when the final result would fit in ordinary matrix operations.
- `GetDeterminant()` selects a cubic-time algorithm for built-in numeric types:
  fraction-free Bareiss with `BigInteger` intermediates for integers and decimal,
  or scaled partial-pivot LU for `double`, `float`, and `Half`. Integer results
  convert back with overflow checking only at the end; decimal results are rounded
  once to the available decimal precision. Float and Half use double working precision.
- Select an algorithm with `GetDeterminant(DeterminantAlgorithm.Lu)`, `Bareiss`,
  `Cholesky`, or `Cofactor`. Cholesky supports binary floats and requires finite,
  exactly symmetric positive-definite input; it does roughly half the factorization
  arithmetic of LU. Cofactor retains the checked division-free implementation for
  small matrices and is the Auto fallback for custom `INumber<T>` types. Its cost
  grows factorially and bounded intermediate arithmetic can overflow.
- `IntegerMatrix` determinants use exact `BigInteger` intermediates. A result
  outside `int` range throws from `GetDeterminant`, while `TryGetDeterminant`
  returns false. Its original addition, subtraction, matrix multiplication, and
  scaling retain unchecked `int` behavior; the new cross product is checked.
- `DoubleMatrix` determinants use elimination with partial pivoting. This is
  faster than cofactor expansion, but rounding still matters, especially for
  nearly singular matrices. No conditioning estimate or approximate-zero
  tolerance is provided. Triangular checks compare entries to zero exactly.
- Legacy `IsInvertible` methods test nonsingularity over the real numbers;
  an integer matrix's inverse need not contain integers. The double version
  checks computed pivots, avoiding a false negative solely from determinant
  underflow. `Matrix<T>` deliberately has no general `IsInvertible` method:
  the meaning depends on the scalar type and numerical precision.
- The default constructor creates a 0×0 matrix, with trace zero and determinant
  one. `Matrix<T>` and the legacy `FromArray` factories also support rectangular
  empty shapes such as 2×0. The legacy dimension-based constructors require
  positive dimensions. Trace and determinant require a square matrix. Generic `TryGet` methods return false
  for incompatible shapes; arithmetic exceptions still propagate.

The tests include independently specified examples, rectangular matrices,
fractional inputs, determinant sign changes under row swaps, overflow cases,
and algebraic identities. They are checks on the implementation, not a claim
that every numerical edge case is handled.

## Eigenvalues and eigenvectors

`Matrix<T>`, `DoubleMatrix`, and `IntegerMatrix` expose
`GetSymmetricEigenDecomposition(tolerance: 1e-12, maxSweeps: 50)`:

```csharp
var symmetric = new Matrix<int>(new[,] { { 1, 1 }, { 1, 0 } });
var eigen = symmetric.GetSymmetricEigenDecomposition();
Console.WriteLine(eigen.EigenValues[1]); // Approximately 1.618033988749895.
double[] direction = eigen.EigenVectors.GetColumn(1);
```

The result contains ascending `double[] EigenValues` and a `DoubleMatrix`
`EigenVectors`. Its **columns** are the corresponding unit eigenvectors:
`A * Q = Q * diag(values)` and `Qᵀ * Q = I`. Results own their storage and the
input is unchanged. Repeated eigenvalues can have different valid orthonormal
bases; eigenvector signs alone do not change the solution.

Inputs must be finite, square, and exactly symmetric. Cyclic Jacobi rotations
support positive-definite, indefinite, singular, repeated-eigenvalue, and empty
matrices. Tolerance must be finite and strictly between zero and one; the sweep
limit must be positive. Failure to converge or a nonfinite eigenvalue throws.
Accuracy is relative to the matrix Frobenius norm, so a tiny eigenvalue beside
much larger entries may have a large relative error. Generic numeric entries
are approximated as doubles, which can lose integer or decimal precision.
For arbitrary finite real **square** matrices, use `GetEigenDecomposition(maxIterations: 1000)`:

```csharp
var rotation = new Matrix<double>(new[,] { { 0d, -1d }, { 1d, 0d } });
var general = rotation.GetEigenDecomposition();
System.Numerics.Complex lambda = general.EigenValues[0]; // -i
System.Numerics.Complex component = general.EigenVectors[0, 0];
```

The result owns a `Complex[] EigenValues` and `Complex[,] EigenVectors`. Values
are sorted by real part, then imaginary part; matching unit right eigenvectors
are columns, satisfying `A v = λ v`. Power-of-two similarity balancing precedes
Hessenberg reduction and real double-shift QR; the iteration bound is per
unconverged root (1–100000). The solver is adapted from the public-domain
[NIST/MathWorks JAMA](https://math.nist.gov/javanumerics/jama/) implementation.
Inputs are preserved. Empty, singular, repeated-spectrum and defective inputs
are supported. Failure to converge or finite-range failure is explicit.

General eigenvectors need not be orthogonal. A defective matrix cannot have a
complete independent eigenbasis, so returned columns may be dependent and do
not imply an invertible diagonalization. Balancing helps disparate scales;
ill-conditioned eigenvalues may still be sensitive to tiny input changes.
This adds complex **outputs**; the matrix entries remain real. Rectangular
matrices do not have an ordinary eigendecomposition.

Legacy `DoubleMatrix.GetEigenValues()` delegates to the symmetric solver, and
`GetEigenVectors()` returns a row-major flattened matrix of eigenvector columns.
Its `TryGetEigenValues` returns false for unsupported inputs or solver failure.
The old `IntegerMatrix` methods returning `int[]` are obsolete: integer matrices
can have noninteger spectra and normalized vectors. Their `Get` methods throw
`NotSupportedException`, and `TryGetEigenValues` returns false with an empty
array. Migrate these callers to `GetSymmetricEigenDecomposition()`.

## Cross products and rotations

`Matrix<T>.CrossProduct(other)` returns a new right-handed 3D vector cross
product. Both operands may independently be 3×1 columns or 1×3 rows; the result
retains the left operand's shape. Inputs are unchanged. Generic arithmetic is
checked, so bounded integer overflow throws, and floating-point inputs and
results must be finite.

```csharp
var x = new Matrix<int>(new[,] { { 1 }, { 0 }, { 0 } });
var y = new Matrix<int>(new[,] { { 0, 1, 0 } });
var z = x.CrossProduct(y); // 3×1 column: [0, 0, 1].
var rotation = MatrixRotation.CreateAxisAngle(z, Math.PI / 2);
```

`MatrixRotation.Create2D`, `CreateX`, `CreateY`, `CreateZ`, and `CreateAxisAngle`
return `DoubleMatrix` rotation matrices on both supported .NET targets.
`CreateAxisAngle` accepts a `DoubleMatrix` axis, with a `Matrix<T>` overload on
.NET 10. Angles are finite **radians**, with right-handed **active column-vector**
conventions: positive 2D angles rotate counterclockwise, and a positive quarter
turn about Z maps X to Y. Arbitrary axes must be finite, nonzero 3×1 or 1×3
vectors; stable normalization handles huge, tiny, and subnormal magnitudes.
The factories construct matrices; apply one with `R * columnVector`, or use
`rowVector * R.Transpose()` for a row vector. No Euler-angle or quaternion API
is exposed.

Legacy `DoubleMatrix` and `IntegerMatrix` offer `GetCrossProduct` for an
independent result. Their existing boolean `CrossProduct` methods now replace
the receiving vector on success and return true. Null, invalid shapes,
nonfinite input/results, or checked integer overflow return false without
changing the receiver. The result-returning forms throw on these failures.

## Development checks

Run the build and test commands above from the repository root, then check
formatting:

```sh
dotnet format --verify-no-changes --no-restore
dotnet format tests/NetStandardSmoke/NetStandardSmoke.csproj --verify-no-changes --no-restore
```

Use `dotnet format --no-restore` to apply the formatting rules. Compiler warnings
and enabled .NET analyzer warnings are treated as errors. When changing NuGet
references, run `dotnet restore` to refresh the committed `packages.lock.json`
files and repeat the checks. Update [global.json](../global.json) when changing
the SDK feature band.

The [C# comparison runner](../benchmarks/dotnet/Program.cs) uses `Matrix<double>`
for the shared float64 workloads and copies inputs before calling mutating
methods. The [benchmark guide](../benchmarks/README.md) describes verification,
algorithm comparisons, timings, and profiles. Contribution guidance for every
implementation is in [CONTRIBUTING.md](../CONTRIBUTING.md).
