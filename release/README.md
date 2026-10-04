# Release kits and external consumers

All eleven implementations are peers. A release consists of a source kit for
each, conventional packages where applicable, a checked compatibility table,
versioned API documentation, checksums, and Git provenance. Packages are built
from extracted kits and tested in a separate temporary consumer directory.
Nothing imports the checkout. Python uses a fresh virtual environment, .NET a
fresh NuGet package directory, JavaScript a fresh npm project, and Julia an
explicit project/load path. Toolchain caches may be shared; package build tools
and FSharp.Core may require network access.

Run from the repository root on ARM64 with the documented language toolchains:

```sh
# Show the release version from manifest.json:
python3 release/build.py --version
python3 release/build.py
# Or exercise only the implementation you need:
python3 release/build.py --port python --output .build/python-release
# Prepare portable source kits without claiming installation verification:
python3 release/build.py --source-only --output .build/source-release
python3 -m unittest discover -s tests/release -v
```

The output directory must be empty to prevent mixing artifacts from different
revisions. `--port` is repeatable. The default output is `.build/releases/<version>`.
Throughout these instructions, `<version>` is the version in
[`manifest.json`](manifest.json), printed by `python3 release/build.py --version`.
`CC` and `CXX` select native compilers; other tools (`dotnet`, `cargo`, `npm`,
`julia`, `go`, `emcc`) are found on PATH. Emscripten needs its environment enabled.
The ARM64 kit requires an ARM64 host; other kits can be built separately.

## Artifacts and installation

First verify `SHA256SUMS` (`sha256sum -c SHA256SUMS`, or `shasum -a 256 -c
SHA256SUMS` on macOS). The index records version, commit, dirty status, selected
ports, installation-check status, byte lengths and hashes. It does not export
environment variables, tool paths, hostnames, or build logs. Source tarballs have
relative members, zero timestamps, and no user/group identifiers. The WASM
package uses the same executable-aware sanitization as the public reports.

| Port | Artifact | Consumer setup |
| --- | --- | --- |
| ARM64 assembly | Source `.tar.gz` | Compile `example/c.c`, `ports/c/matrix.c`, `ports/assembly/kernels.S`, `-I ports/c -DMATRIX_USE_ASM -std=c11 -lm` |
| C | Source `.tar.gz` | Compile `example/c.c` and `ports/c/matrix.c`, `-I ports/c -std=c11 -lm` |
| C++ | Source `.tar.gz` | Compile `example/cpp.cpp`, `-I ports/cpp -std=c++17`; retain sibling `ports/c` headers |
| C# | `LinearA.CSharp.<version>.nupkg` + source | Add a local NuGet source, reference `LinearA.CSharp` version `<version>` in a .NET 10 project |
| F# | `LinearA.FSharp.<version>.nupkg` + source | Same, using `LinearA.FSharp`; FSharp.Core is resolved by NuGet |
| Go | Source `.tar.gz` | Use a module `replace github.com/jmcummings77/linear-A/ports/go => /path/to/extracted/ports/go` |
| Julia | Source `.tar.gz` | In a fresh project, `using Pkg; Pkg.develop(path="/path/to/extracted/ports/julia")` |
| Python | `linear_a_python-<version>-*.whl` + source | `python3 -m venv .venv`; `.venv/bin/pip install /path/to/package.whl` |
| Rust | `linear-a-rust-<version>.crate` + source | Extract crate, use `linear-a-rust = { path = "/path/to/extracted/crate" }`; imported library is `linear_a` |
| TypeScript | `linear-a-typescript-<version>.tgz` + source | `npm install /path/to/package.tgz`; import `Matrix` from `linear-a-typescript` |
| WebAssembly | `linear-a-wasm-<version>.tgz` + source | `npm install /path/to/package.tgz`; import `createMatrixAPI` from `linear-a-wasm` |

Use `-ffp-contract=off` for native compilation to match the shared arithmetic
policy. Examples under `release/examples/` are the exact programs run by the
installation checks. Build/install orchestration is in `build.py`; artifacts
contain their corresponding example and API docs as well.

For a browser, serve the WASM package files over HTTP and import `index.mjs` by
URL. Keep `runtime/matrix.mjs` and `runtime/matrix.wasm` together. The runtime URL
is resolved relative to the wrapper, so moving the package does not break it.
You can supply `moduleUrl`, `locateFile`, or `wasmBinary` to customize hosting.
Dispose native allocations explicitly. TypeScript also requires a build step;
the npm archive contains compiled JavaScript and declarations.

These names and versions describe local artifacts; no package-registry
publication is implied. A root `v<version>` Git tag does not publish the Go submodule:
Go registry resolution needs a separate `ports/go/v<version>` tag. Julia is currently
consumed by path or repository subdirectory, not through General.

## Release maintenance

`manifest.json` owns the release version and shared capability declarations.
Keep the Python, Rust, Julia, and TypeScript versions equal to it. .NET and WASM
package versions are supplied from the manifest. Update `docs/api/<version>/`,
then regenerate its compatibility table using `python3 release/build.py
--write-table`. Earlier API documents remain versioned. The table is checked
before every build, while the shared port conformance suite checks capabilities.

CI builds and tests the complete artifact set on ARM64, then uploads it as a
workflow artifact. A manually dispatched run does the same without publishing a
release. Pushing a matching `v<version>` tag builds the kits and attaches them to
a **draft GitHub release** after all consumer checks pass. Review the draft's
assets and provenance before publishing. No package-registry tokens are used.
Do not tag a release until the shared port conformance workflow is green on the
same commit. Source-only builds are explicitly marked `not-run` and are never
used by the release workflow.
