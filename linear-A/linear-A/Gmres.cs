#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;

namespace linear_A;

/// <summary>Owned GMRES solution, residual histories, optional snapshots, and restart indices.</summary>
public sealed record GMRESResult(
    double[] X,
    bool Converged,
    int Iterations,
    string Reason,
    double[] Residuals,
    double[] EstimatedResiduals,
    double[][] Iterates,
    int[] Restarts);

public sealed partial class CSRMatrix
{
    /// <summary>Restarted GMRES from a zero initial guess, with optional right preconditioning.</summary>
    public GMRESResult Gmres(
        double[] b,
        int restart = 30,
        double relativeTolerance = 1e-10,
        double absoluteTolerance = 0,
        int maxIterations = 1000,
        bool jacobi = false,
        bool capture = false,
        ILU0? preconditioner = null)
    {
        // Validate inputs before preparing the right preconditioner.
        ArgumentNullException.ThrowIfNull(b);
        int n = Rows;
        if (Cols != n || b.Length != n || b.Any(value => !double.IsFinite(value)))
        {
            throw new ArgumentException("GMRES requires square matrix and finite matching vector.");
        }
        if (restart < 1 || restart > 1024 || maxIterations < 0 || maxIterations > 100000 ||
            !double.IsFinite(relativeTolerance) || relativeTolerance < 0 || relativeTolerance >= 1 ||
            !double.IsFinite(absoluteTolerance) || absoluteTolerance < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(restart), "Invalid GMRES options.");
        }
        if (preconditioner is not null && (jacobi || preconditioner.Size != n))
        {
            throw new ArgumentException("Incompatible preconditioner.");
        }

        var diagonal = Enumerable.Repeat(1.0, n).ToArray();
        if (jacobi)
        {
            for (int row = 0; row < n; row++)
            {
                bool found = false;
                for (int entry = offsets[row]; entry < offsets[row + 1]; entry++)
                {
                    if (indices[entry] == row)
                    {
                        diagonal[row] = values[entry];
                        found = true;
                        break;
                    }
                }
                if (!found || diagonal[row] == 0)
                {
                    throw new ArgumentException("Jacobi requires a nonzero diagonal.");
                }
            }
        }

        double[] ApplyRightPreconditioner(double[] vector)
        {
            return preconditioner is null
                ? vector.Select((value, index) => value / diagonal[index]).ToArray()
                : preconditioner.Apply(vector);
        }

        // Initialize the true residual and histories at x = 0. Only accepted
        // candidates advance the iteration count or become captured frames.
        var x = new double[n];
        var residual = (double[])b.Clone();
        var residualHistory = new List<double> { Norm(residual) };
        var estimatedResidualHistory = new List<double>(residualHistory);
        var iterates = new List<double[]>();
        var restarts = new List<int>();
        if (capture)
        {
            iterates.Add((double[])x.Clone());
        }

        GMRESResult Result(string reason) => new(
            (double[])x.Clone(),
            reason == "converged",
            residualHistory.Count - 1,
            reason,
            residualHistory.ToArray(),
            estimatedResidualHistory.ToArray(),
            iterates.ToArray(),
            restarts.ToArray());

        double threshold = Math.Max(absoluteTolerance, relativeTolerance * residualHistory[0]);
        int cycleSize = Math.Min(restart, Math.Min(n, maxIterations));
        if (!double.IsFinite(residualHistory[0]))
        {
            return Result("nonfinite");
        }
        if (residualHistory[0] <= threshold)
        {
            return Result("converged");
        }

        while (residualHistory.Count - 1 < maxIterations)
        {
            // Start a restart cycle from the last accepted true residual.
            if (residualHistory.Count > 1)
            {
                restarts.Add(residualHistory.Count - 1);
            }
            var cycleStart = (double[])x.Clone();
            double beta = Norm(residual);
            var basis = new List<double[]> { residual.Select(value => value / beta).ToArray() };
            var hessenberg = new double[cycleSize + 1, cycleSize];
            var cosines = new double[cycleSize];
            var sines = new double[cycleSize];
            var leastSquaresRhs = new double[cycleSize + 1];
            leastSquaresRhs[0] = beta;
            int steps = Math.Min(cycleSize, maxIterations - (residualHistory.Count - 1));

            for (int column = 0; column < steps; column++)
            {
                // Arnoldi expansion for A * M^-1. Right preconditioning keeps
                // the residual in the original system's coordinates.
                double[] work;
                try
                {
                    work = Matvec(ApplyRightPreconditioner(basis[column]));
                }
                catch (ArithmeticException)
                {
                    return Result("nonfinite");
                }
                catch (ArgumentException)
                {
                    // Matvec rejects nonfinite vectors produced by Jacobi division.
                    return Result("nonfinite");
                }
                double originalNorm = Norm(work);
                Reorthogonalize(basis, work, hessenberg, column);
                double remainderNorm = Norm(work);
                if (!double.IsFinite(originalNorm) || !double.IsFinite(remainderNorm))
                {
                    return Result("nonfinite");
                }
                for (int row = 0; row <= column; row++)
                {
                    if (!double.IsFinite(hessenberg[row, column]))
                    {
                        return Result("nonfinite");
                    }
                }

                // A remainder at roundoff scale cannot provide a reliable new
                // direction. This is "happy" breakdown only if the true residual
                // later confirms convergence. double.Epsilon is not machine epsilon.
                const double machineEpsilon = 2.220446049250313e-16;
                bool arnoldiBreakdown = remainderNorm <= 8 * machineEpsilon * originalNorm;
                hessenberg[column + 1, column] = arnoldiBreakdown ? 0 : remainderNorm;
                if (!arnoldiBreakdown)
                {
                    basis.Add(work.Select(value => value / remainderNorm).ToArray());
                }

                // Update the least-squares QR factorization. Hypot avoids the
                // unnecessary overflow/underflow of squaring both pivot entries.
                ApplyPreviousRotations(hessenberg, cosines, sines, column);
                double pivot = double.Hypot(hessenberg[column, column], hessenberg[column + 1, column]);
                if (!double.IsFinite(pivot))
                {
                    return Result("nonfinite");
                }
                if (pivot == 0)
                {
                    return Result("breakdown");
                }

                cosines[column] = hessenberg[column, column] / pivot;
                sines[column] = hessenberg[column + 1, column] / pivot;
                hessenberg[column, column] = pivot;
                hessenberg[column + 1, column] = 0;
                leastSquaresRhs[column + 1] = -sines[column] * leastSquaresRhs[column];
                leastSquaresRhs[column] = cosines[column] * leastSquaresRhs[column];

                // Solve the reduced triangular system, then form x_base + M^-1 V y.
                var coefficients = leastSquaresRhs[..(column + 1)];
                if (!BackSubstitute(hessenberg, coefficients, column))
                {
                    return Result("breakdown");
                }
                var candidate = CombineBasis(basis, coefficients, column);
                try
                {
                    candidate = ApplyRightPreconditioner(candidate);
                }
                catch (ArithmeticException)
                {
                    return Result("nonfinite");
                }
                for (int row = 0; row < n; row++)
                {
                    candidate[row] += cycleStart[row];
                }
                if (coefficients.Any(value => !double.IsFinite(value)) ||
                    candidate.Any(value => !double.IsFinite(value)) ||
                    !double.IsFinite(leastSquaresRhs[column + 1]))
                {
                    return Result("nonfinite");
                }

                // Recompute the true residual before accepting the candidate.
                // Roundoff can make the projected estimate overly optimistic;
                // it is recorded for diagnostics but never establishes convergence.
                double[] product;
                try
                {
                    product = Matvec(candidate);
                }
                catch (ArithmeticException)
                {
                    return Result("nonfinite");
                }
                catch (ArgumentException)
                {
                    return Result("nonfinite");
                }
                var candidateResidual = b.Select((value, index) => value - product[index]).ToArray();
                double residualNorm = Norm(candidateResidual);
                if (!double.IsFinite(residualNorm))
                {
                    return Result("nonfinite");
                }

                x = candidate;
                residual = candidateResidual;
                residualHistory.Add(residualNorm);
                estimatedResidualHistory.Add(Math.Abs(leastSquaresRhs[column + 1]));
                if (capture)
                {
                    iterates.Add((double[])x.Clone());
                }
                if (residualNorm <= threshold)
                {
                    return Result("converged");
                }
                if (arnoldiBreakdown)
                {
                    return Result("breakdown");
                }
            }

            // An exactly unchanged completed cycle takes precedence over the
            // iteration limit. Do not introduce another numerical tolerance here.
            if (x.SequenceEqual(cycleStart))
            {
                return Result("stagnation");
            }
        }
        return Result("iteration_limit");

        static double Norm(double[] vector)
        {
            // Accumulating with Hypot keeps representable norms usable even when
            // individual squares would overflow or underflow.
            double norm = 0;
            foreach (double value in vector)
            {
                norm = double.Hypot(norm, value);
            }
            return norm;
        }

        static void Reorthogonalize(List<double[]> basis, double[] work, double[,] hessenberg, int column)
        {
            // Two modified Gram-Schmidt passes reduce loss of orthogonality as
            // Krylov vectors become nearly dependent. Both passes contribute to H.
            for (int pass = 0; pass < 2; pass++)
            {
                for (int vector = 0; vector <= column; vector++)
                {
                    double projection = 0;
                    for (int row = 0; row < work.Length; row++)
                    {
                        projection += basis[vector][row] * work[row];
                    }
                    hessenberg[vector, column] += projection;
                    for (int row = 0; row < work.Length; row++)
                    {
                        work[row] -= projection * basis[vector][row];
                    }
                }
            }
        }

        static void ApplyPreviousRotations(double[,] hessenberg, double[] cosines, double[] sines, int column)
        {
            for (int row = 0; row < column; row++)
            {
                double top = cosines[row] * hessenberg[row, column] + sines[row] * hessenberg[row + 1, column];
                hessenberg[row + 1, column] = -sines[row] * hessenberg[row, column] + cosines[row] * hessenberg[row + 1, column];
                hessenberg[row, column] = top;
            }
        }

        static bool BackSubstitute(double[,] upper, double[] coefficients, int column)
        {
            for (int row = column; row >= 0; row--)
            {
                if (upper[row, row] == 0)
                {
                    return false;
                }
                double sum = 0;
                for (int next = row + 1; next <= column; next++)
                {
                    sum += upper[row, next] * coefficients[next];
                }
                coefficients[row] = (coefficients[row] - sum) / upper[row, row];
            }
            return true;
        }

        static double[] CombineBasis(List<double[]> basis, double[] coefficients, int column)
        {
            var correction = new double[basis[0].Length];
            for (int row = 0; row < correction.Length; row++)
            {
                double sum = 0;
                for (int vector = 0; vector <= column; vector++)
                {
                    sum += basis[vector][row] * coefficients[vector];
                }
                correction[row] = sum;
            }
            return correction;
        }
    }
}
#endif
