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
        if (op is not ("spmv" or "dense" or "cg" or "gmres" or "ilu_setup" or "ilu_apply") || rows < 0 || cols < 0 || nnz < 0 || iterations < 0 || (jacobi < 0 || jacobi > (op == "gmres" ? 3 : 1)) || capture is < 0 or > 1) throw new ArgumentException("Invalid sparse options.");
        var tokens = Console.In.ReadToEnd().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries); var count = op is "cg" or "gmres" ? rows : cols;
        if (tokens.Length != checked(rows + 1 + 2 * nnz + count)) throw new ArgumentException("Incorrect sparse input count.");
        var offset = 0; int[] Indices(int n) => Enumerable.Range(0, n).Select(_ => int.Parse(tokens[offset++], CultureInfo.InvariantCulture)).ToArray(); double[] Values(int n) => Enumerable.Range(0, n).Select(_ => double.Parse(tokens[offset++], CultureInfo.InvariantCulture)).ToArray();
        var rp = Indices(rows + 1); var ci = Indices(nnz); var v = Values(nnz); var b = Values(count); var a = new CSRMatrix(rows, cols, rp, ci, v);
        var factor = op == "ilu_apply" || (op == "gmres" && jacobi == 2) ? new ILU0(a) : null;
        Matrix<double>? dense = null, right = null;
        if (op == "dense") { dense = new Matrix<double>(rows, cols); for (var i = 0; i < rows; i++) for (var p = rp[i]; p < rp[i + 1]; p++) dense[i, ci[p]] = v[p]; right = new Matrix<double>(cols, 1); for (var i = 0; i < cols; i++) right[i, 0] = b[i]; }
        double[] Compute()
        {
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
