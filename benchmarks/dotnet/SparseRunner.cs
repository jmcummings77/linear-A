using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using linear_A;
internal static class SparseRunner
{
    internal static void Run(string[] args)
    {
        if (args.Length != (args.Length > 1 && args[1] == "gmres" ? 12 : 11)) throw new ArgumentException("Invalid sparse protocol.");
        var op = args[1]; int Int(int i) => int.Parse(args[i], CultureInfo.InvariantCulture); double Number(int i) => double.Parse(args[i], CultureInfo.InvariantCulture);
        var rows = Int(2); var cols = Int(3); var nnz = Int(4); var iterations = Int(5); var rtol = Number(6); var atol = Number(7); var limit = Int(8); var jacobi = Int(9); var capture = Int(10);
        if (op is not ("spmv" or "dense" or "cg" or "gmres" or "ilu_setup" or "ilu_apply" or "rcm" or "permute" or "permutation_check" or "rcm_solve" or "ilu_solve") || rows < 0 || cols < 0 || nnz < 0 || iterations < 0 || (jacobi < 0 || jacobi > (op == "gmres" ? 3 : 1)) || capture is < 0 or > 1) throw new ArgumentException("Invalid sparse options.");
        var tokens = Console.In.ReadToEnd().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries); var count = op is "cg" or "gmres" ? rows : cols;
        if (tokens.Length != checked(rows + 1 + 2 * nnz + count)) throw new ArgumentException("Incorrect sparse input count.");
        var offset = 0; int[] Indices(int n) => Enumerable.Range(0, n).Select(_ => int.Parse(tokens[offset++], CultureInfo.InvariantCulture)).ToArray(); double[] Values(int n) => Enumerable.Range(0, n).Select(_ => double.Parse(tokens[offset++], CultureInfo.InvariantCulture)).ToArray();
        var rp = Indices(rows + 1); var ci = Indices(nnz); var v = Values(nnz); var b = Values(count); var a = new CSRMatrix(rows, cols, rp, ci, v);
        var factor = op == "ilu_apply" || (op == "gmres" && jacobi == 2) ? new ILU0(a) : null;
        Matrix<double>? dense = null, right = null;
        if (op == "dense") { dense = new Matrix<double>(rows, cols); for (var i = 0; i < rows; i++) for (var p = rp[i]; p < rp[i + 1]; p++) dense[i, ci[p]] = v[p]; right = new Matrix<double>(cols, 1); for (var i = 0; i < cols; i++) right[i, 0] = b[i]; }
        double[] Compute()
        {
            if (op is "rcm_solve" or "ilu_solve")
            {
                var p = op == "rcm_solve" ? a.ReverseCuthillMcKee() : Enumerable.Range(0, rows).ToArray(); var q = op == "rcm_solve" ? a.PermuteSymmetric(p) : a; var rhs = op == "rcm_solve" ? CSRMatrix.PermuteVector(p, b) : b;
                var orderedResult = q.Gmres(rhs, 20, rtol, atol, limit, false, false, new ILU0(q)); if (!orderedResult.Converged) throw new ArithmeticException("Ordering solve failed: " + orderedResult.Reason); var x = op == "rcm_solve" ? CSRMatrix.PermuteVector(p, orderedResult.X, true) : orderedResult.X; return new double[] { orderedResult.Iterations }.Concat(x).ToArray();
            }
            if (op == "rcm") return a.ReverseCuthillMcKee().Select(i => (double)i).ToArray();
            if (op is "permute" or "permutation_check")
            {
                if (b.Any(x => !double.IsFinite(x) || x < 0 || x >= rows || x != Math.Floor(x))) throw new ArgumentException("Invalid permutation.");
                var p = b.Select(x => (int)x).ToArray(); var q = a.PermuteSymmetric(p);
                if (op == "permute") return q.RowOffsets.Select(x => (double)x).Concat(q.ColumnIndices.Select(x => (double)x)).Concat(q.Values).ToArray();
                var x = Enumerable.Range(1, rows).Select(i => (double)i).ToArray(); var y = CSRMatrix.PermuteVector(p, x);
                return y.Concat(CSRMatrix.PermuteVector(p, y, true)).Concat(CSRMatrix.PermuteVector(p, q.Matvec(y), true)).ToArray();
            }
            if (op == "ilu_setup") { var f = new ILU0(a); return new double[] { f.Size, f.NNZ }; }
            if (op == "ilu_apply") return factor!.Apply(b);
            if (op == "spmv") return a.Matvec(b);
            if (op == "dense") { var value = dense! * right!; return Enumerable.Range(0, rows).Select(i => value[i, 0]).ToArray(); }
            if (op == "gmres")
            {
                var g = a.Gmres(b, Int(11), rtol, atol, limit, jacobi == 1, capture != 0, jacobi == 3 ? new ILU0(a) : factor);
                if (iterations > 0) { if (!g.Converged) throw new ArithmeticException("Benchmark GMRES did not converge: " + g.Reason); return g.X; }
                var packed = new List<double> { Array.IndexOf(new[] { "converged", "iteration_limit", "breakdown", "nonfinite", "stagnation" }, g.Reason), g.Iterations, g.Residuals.Length };
                packed.AddRange(g.X); packed.AddRange(g.Residuals); packed.AddRange(g.EstimatedResiduals); packed.Add(g.Restarts.Length); packed.AddRange(g.Restarts.Select(v => (double)v)); foreach (var frame in g.Iterates) packed.AddRange(frame); return packed.ToArray();
            }
            var r = a.ConjugateGradient(b, rtol, atol, limit, jacobi == 1, capture != 0);
            if (iterations > 0) { if (!r.Converged) throw new ArithmeticException("Benchmark CG did not converge: " + r.Reason); return r.X; }
            var result = new List<double> { Array.IndexOf(new[] { "converged", "iteration_limit", "breakdown", "nonfinite" }, r.Reason), r.Iterations, r.Residuals.Length };
            result.AddRange(r.X); result.AddRange(r.Residuals); foreach (var frame in r.Iterates) result.AddRange(frame); return result.ToArray();
        }
        if (iterations == 0) { var values = Compute(); Console.WriteLine(JsonSerializer.Serialize(new { rows = 1, cols = values.Length, values })); return; }
        for (var i = 0; i < 3; i++) Compute(); var checksum = 0.0; var timer = Stopwatch.StartNew();
        for (var i = 0; i < iterations; i++) foreach (var x in Compute()) checksum += x;
        var elapsed = timer.Elapsed.TotalNanoseconds;
        if (!double.IsFinite(checksum)) throw new ArithmeticException("Nonfinite benchmark checksum.");
        Console.WriteLine(JsonSerializer.Serialize(new { elapsed_ns = elapsed, iterations, checksum }));
    }
}
