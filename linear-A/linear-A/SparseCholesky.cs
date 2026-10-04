#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;

namespace linear_A;

/// <summary>Reusable symbolic lower pattern. FillSteps uses -1 for original entries.</summary>
public sealed class SparseCholeskySymbolic
{
    private readonly int[] sourceRowOffsets;
    private readonly int[] sourceColumnIndices;
    private readonly int[] rowOffsets;
    private readonly int[] columnIndices;
    private readonly int[] fillSteps;

    public int Size { get; }
    public int NNZ => columnIndices.Length;
    public int FillCount => fillSteps.Count(step => step >= 0);
    public int[] RowOffsets => (int[])rowOffsets.Clone();
    public int[] ColumnIndices => (int[])columnIndices.Clone();
    public int[] FillSteps => (int[])fillSteps.Clone();

    public SparseCholeskySymbolic(CSRMatrix a) : this(a, false)
    {
    }

    internal SparseCholeskySymbolic(CSRMatrix a, bool incomplete)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.Rows != a.Cols)
        {
            throw new ArgumentException("Cholesky requires square matrix.");
        }
        Size = a.Rows;
        sourceRowOffsets = a.RowOffsets;
        sourceColumnIndices = a.ColumnIndices;

        // Symbolic analysis uses the symmetric union of the stored pattern.
        // Numeric symmetry is checked separately when values are factorized.
        var adjacency = CreateSymmetricAdjacency();
        if (!incomplete)
        {
            AddEliminationFill(adjacency);
        }

        // IC(0) skips elimination fill. Both variants store sorted lower rows
        // ending in a diagonal, even if the source has no stored diagonal.
        rowOffsets = new int[Size + 1];
        var columns = new List<int>();
        var entryFillSteps = new List<int>();
        for (int row = 0; row < Size; row++)
        {
            foreach (var neighbor in adjacency[row])
            {
                if (neighbor.Key < row)
                {
                    columns.Add(neighbor.Key);
                    entryFillSteps.Add(neighbor.Value);
                }
            }
            columns.Add(row);
            entryFillSteps.Add(-1);
            rowOffsets[row + 1] = columns.Count;
        }
        columnIndices = columns.ToArray();
        fillSteps = entryFillSteps.ToArray();
    }

    private SortedDictionary<int, int>[] CreateSymmetricAdjacency()
    {
        var adjacency = Enumerable.Range(0, Size)
            .Select(_ => new SortedDictionary<int, int>()).ToArray();
        for (int row = 0; row < Size; row++)
        {
            for (int entry = sourceRowOffsets[row]; entry < sourceRowOffsets[row + 1]; entry++)
            {
                int column = sourceColumnIndices[entry];
                if (row != column)
                {
                    adjacency[row][column] = -1;
                    adjacency[column][row] = -1;
                }
            }
        }
        return adjacency;
    }

    private static void AddEliminationFill(SortedDictionary<int, int>[] adjacency)
    {
        for (int pivot = 0; pivot < adjacency.Length; pivot++)
        {
            // Eliminating a vertex connects its remaining neighbors. Record
            // only the first pivot that creates each edge for FillSteps.
            var neighbors = adjacency[pivot].Keys.Where(column => column > pivot).ToArray();
            for (int right = 0; right < neighbors.Length; right++)
            {
                for (int left = 0; left < right; left++)
                {
                    int row = neighbors[right];
                    int column = neighbors[left];
                    if (!adjacency[row].ContainsKey(column))
                    {
                        adjacency[row][column] = pivot;
                        adjacency[column][row] = pivot;
                    }
                }
            }
        }
    }

    public SparseCholesky Factorize(CSRMatrix a)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.Rows != Size || a.Cols != Size ||
            !a.RowOffsets.SequenceEqual(sourceRowOffsets) ||
            !a.ColumnIndices.SequenceEqual(sourceColumnIndices))
        {
            throw new ArgumentException("Cholesky symbolic pattern mismatch.");
        }

        var sourceRows = ReadSymmetricValues(a.Values);
        var factorValues = new double[columnIndices.Length];
        for (int row = 0; row < Size; row++)
        {
            for (int entry = rowOffsets[row]; entry < rowOffsets[row + 1]; entry++)
            {
                int column = columnIndices[entry];
                double updatedValue = SubtractKnownProducts(
                    row, entry, sourceRows[row].GetValueOrDefault(column), factorValues);
                RequireFiniteFactor(updatedValue);

                if (row == column)
                {
                    // No pivoting, shifts, or tolerance repair is applied. IC(0)
                    // can fail this test even for SPD input because it drops fill.
                    if (updatedValue <= 0)
                    {
                        throw new ArithmeticException("Nonpositive Cholesky pivot.");
                    }
                    factorValues[entry] = Math.Sqrt(updatedValue);
                }
                else
                {
                    int pivotEntry = rowOffsets[column + 1] - 1;
                    factorValues[entry] = updatedValue / factorValues[pivotEntry];
                }
                RequireFiniteFactor(factorValues[entry]);
            }
        }
        return new SparseCholesky(new CSRMatrix(Size, Size, rowOffsets, columnIndices, factorValues));
    }

    private Dictionary<int, double>[] ReadSymmetricValues(double[] sourceValues)
    {
        var sourceRows = Enumerable.Range(0, Size)
            .Select(_ => new Dictionary<int, double>()).ToArray();
        for (int row = 0; row < Size; row++)
        {
            for (int entry = sourceRowOffsets[row]; entry < sourceRowOffsets[row + 1]; entry++)
            {
                sourceRows[row][sourceColumnIndices[entry]] = sourceValues[entry];
            }
        }

        // Require exact numeric symmetry, treating absent entries as zero.
        // A stored zero does not require a matching stored transpose entry.
        for (int row = 0; row < Size; row++)
        {
            foreach (var entry in sourceRows[row])
            {
                if (entry.Value != sourceRows[entry.Key].GetValueOrDefault(row))
                {
                    throw new ArgumentException("Cholesky requires symmetric values.");
                }
            }
        }
        return sourceRows;
    }

    private double SubtractKnownProducts(int row, int entry, double value, double[] factorValues)
    {
        int column = columnIndices[entry];
        int rowEntry = rowOffsets[row];
        int columnEntry = rowOffsets[column];
        int columnDiagonal = rowOffsets[column + 1] - 1;

        // Intersect the sorted prefixes to subtract L[row,k] * L[column,k].
        // Absent entries contribute zero, preserving the chosen fill policy.
        while (rowEntry < entry && columnEntry < columnDiagonal)
        {
            if (columnIndices[rowEntry] == columnIndices[columnEntry])
            {
                value -= factorValues[rowEntry] * factorValues[columnEntry];
                rowEntry++;
                columnEntry++;
            }
            else if (columnIndices[rowEntry] < columnIndices[columnEntry])
            {
                rowEntry++;
            }
            else
            {
                columnEntry++;
            }
        }
        return value;
    }

    private static void RequireFiniteFactor(double value)
    {
        // Keep finite-range failure distinct from a nonpositive finite pivot.
        if (!double.IsFinite(value))
        {
            throw new ArithmeticException("Nonfinite Cholesky factor.");
        }
    }
}

/// <summary>Owned lower factor supporting repeated right-hand sides.</summary>
public sealed class SparseCholesky
{
    private readonly CSRMatrix lower;

    internal SparseCholesky(CSRMatrix lower)
    {
        this.lower = lower;
    }

    public int Size => lower.Rows;
    public int NNZ => lower.NNZ;
    public CSRMatrix Lower => new(Size, Size, lower.RowOffsets, lower.ColumnIndices, lower.Values);

    public double[] Solve(double[] b)
    {
        ArgumentNullException.ThrowIfNull(b);
        if (b.Length != Size || b.Any(value => !double.IsFinite(value)))
        {
            throw new ArgumentException("Invalid Cholesky right-hand side.");
        }
        var x = (double[])b.Clone();
        var rowOffsets = lower.RowOffsets;
        var columnIndices = lower.ColumnIndices;
        var values = lower.Values;

        // Forward substitution: L y = b. Each sorted lower row ends in its diagonal.
        for (int row = 0; row < Size; row++)
        {
            int diagonal = rowOffsets[row + 1] - 1;
            for (int entry = rowOffsets[row]; entry < diagonal; entry++)
            {
                x[row] -= values[entry] * x[columnIndices[entry]];
            }
            x[row] /= values[diagonal];
            RequireFiniteSolution(x[row]);
        }

        // Back substitution: L^T x = y. Scatter each solved row through L's
        // lower entries, avoiding a separate transposed factor and its storage.
        for (int row = Size - 1; row >= 0; row--)
        {
            int diagonal = rowOffsets[row + 1] - 1;
            x[row] /= values[diagonal];
            RequireFiniteSolution(x[row]);
            for (int entry = rowOffsets[row]; entry < diagonal; entry++)
            {
                int column = columnIndices[entry];
                x[column] -= values[entry] * x[row];
                RequireFiniteSolution(x[column]);
            }
        }
        return x;
    }

    private static void RequireFiniteSolution(double value)
    {
        if (!double.IsFinite(value))
        {
            throw new ArithmeticException("Nonfinite Cholesky solve.");
        }
    }
}

/// <summary>Owned zero-fill incomplete Cholesky; no shifts or pivoting.</summary>
public sealed class IC0
{
    private readonly SparseCholesky factor;

    public IC0(CSRMatrix a)
    {
        factor = new SparseCholeskySymbolic(a, incomplete: true).Factorize(a);
    }

    public int Size => factor.Size;
    public int NNZ => factor.NNZ;
    public CSRMatrix Lower => factor.Lower;
    public double[] Apply(double[] b) => factor.Solve(b);
}
#endif
