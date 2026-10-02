# Contributing

This is a hobby project for learning linear algebra and C# through tests.
Small fixes, better examples, and suggestions for useful exercises are welcome.

Use the .NET SDK selected by `global.json`. Work on a branch and open a pull
request against `main` with a short explanation of the change and how it was
checked. Before submitting, run these commands from the repository root:

```sh
dotnet restore --locked-mode
dotnet build --configuration Release --no-restore
dotnet test --configuration Release --no-build --no-restore
dotnet format --verify-no-changes --no-restore
```

Run `dotnet format --no-restore` to apply the formatting rules. Compiler warnings
and enabled .NET analyzer warnings are treated as errors. GitHub Actions builds
both library targets and runs the tests on Linux, Windows, and macOS.

For math changes, include a small worked example with expected values calculated
independently. Add a regression test for a bug before changing the algorithm.
Use rectangular and nonsymmetric matrices where appropriate: square or identity
matrices can hide indexing mistakes. For floating-point assertions, choose a
tolerance suitable for the values under test and explain it when it is not
obvious. Keep integer overflow, floating-point rounding, and shape errors distinct.

`Matrix<T>` is the generic implementation on .NET 10. `IntegerMatrix` and
`DoubleMatrix` retain the original API and also build for .NET Standard 2.1.
Document any intentional difference in their behavior.

When changing dependencies, update the project references and run `dotnet restore`
to refresh the committed `packages.lock.json` files, then repeat the checks above.
Update `global.json` when moving to a newer SDK feature band. Dependabot checks
NuGet packages and GitHub Actions monthly.
