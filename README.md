# linear-A

A small C# linear algebra library built for fun and for practicing TDD. This is a
learning project, with no intended production use case.

The library provides mutable, fixed-size matrices with addition, subtraction,
matrix multiplication, scalar multiplication, transpose, trace, determinants,
and triangular-matrix checks. Indices are zero-based.

## Getting started

Install the [.NET 10 SDK](https://dotnet.microsoft.com/download/dotnet/10.0).
The SDK version is recorded in `global.json`. From the repository root:

```sh
dotnet restore --locked-mode
dotnet build --configuration Release --no-restore
dotnet test --configuration Release --no-build --no-restore
```

The library targets both .NET 10 and .NET Standard 2.1. `IntegerMatrix` and
`DoubleMatrix` remain available on both targets; the new `Matrix<T>` requires
the .NET 10 target. Tests run on .NET 10. The library has no external runtime
dependencies.

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

## Math and numeric behavior

- `Matrix<T>` uses checked arithmetic where the numeric type supports it.
  Bounded integer operations can throw `OverflowException`; use
  `Matrix<BigInteger>` for unbounded integer arithmetic. Floating-point
  calculations still round and can overflow or underflow. Intermediate values
  can overflow even when the final result would fit; unsigned types can also
  fail when a cofactor requires a negative value.
- Generic determinants use division-free cofactor expansion, so integer
  division cannot truncate intermediate results. This is straightforward but
  has factorial cost; it is intended for small matrices.
- `IntegerMatrix` determinants use exact `BigInteger` intermediates. A result
  outside `int` range throws from `GetDeterminant`, while `TryGetDeterminant`
  returns false. Its other arithmetic retains the original unchecked `int`
  behavior.
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
  one. `Matrix<T>` also supports rectangular empty shapes such as 2×0. The
  legacy dimension-based constructors require positive dimensions. Trace and
  determinant require a square matrix. Generic `TryGet` methods return false
  for incompatible shapes; arithmetic exceptions still propagate.

The tests include independently specified examples, rectangular matrices,
fractional inputs, determinant sign changes under row swaps, overflow cases,
and algebraic identities. They are checks on the implementation, not a claim
that every numerical edge case is handled.

## Still on the wishlist

- Eigenvalues and eigenvectors. The old classes still contain unimplemented
  methods; their `Get` methods throw `NotImplementedException` for square
  matrices, and `TryGetEigenValues` also throws for square matrices.
- Cross products and rotations. The old `CrossProduct` method is a placeholder
  that always returns false.
- Benchmarks and faster generic determinant algorithms.
- Revisit the duplicated legacy classes once compatibility is no longer useful.

Contributions and suggested exercises are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md)
for the development workflow.
