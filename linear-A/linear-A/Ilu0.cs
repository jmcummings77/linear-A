#if NET10_0_OR_GREATER
using System;
using System.Linq;

namespace linear_A;

/// <summary>Owned ILU(0), without pivoting, fill or diagonal shifts.</summary>
public sealed class ILU0
{
    private readonly int[] rowOffsets;
    private readonly int[] columnIndices;
    private readonly int[] diagonalEntries;
    private readonly double[] values;

    public int Size { get; }
    public int NNZ => values.Length;

    public ILU0(CSRMatrix a)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.Rows != a.Cols)
        {
            throw new ArgumentException("ILU0 requires square matrix.");
        }

        // CSR accessors return copies. The factor owns the original sparsity
        // pattern and overwrites only its private copy of the numeric entries.
        Size = a.Rows;
        rowOffsets = a.RowOffsets;
        columnIndices = a.ColumnIndices;
        values = a.Values;
        diagonalEntries = new int[Size];

        LocateDiagonals();
        FactorizeInPlace();
    }

    private void LocateDiagonals()
    {
        for (int row = 0; row < Size; row++)
        {
            int diagonal = FindEntry(row, row);
            if (diagonal == rowOffsets[row + 1] || columnIndices[diagonal] != row)
            {
                throw new ArgumentException("ILU0 requires stored diagonal.");
            }
            diagonalEntries[row] = diagonal;
        }
    }

    private void FactorizeInPlace()
    {
        for (int row = 0; row < Size; row++)
        {
            // Store L multipliers below the diagonal and U on/above it.
            // Earlier rows already have nonzero pivots; L has an implicit unit diagonal.
            for (int entry = rowOffsets[row]; entry < diagonalEntries[row]; entry++)
            {
                int pivotRow = columnIndices[entry];
                values[entry] /= values[diagonalEntries[pivotRow]];
                RequireFinite(values[entry]);

                for (int upperEntry = diagonalEntries[pivotRow] + 1;
                    upperEntry < rowOffsets[pivotRow + 1]; upperEntry++)
                {
                    int column = columnIndices[upperEntry];
                    int targetEntry = FindEntry(row, column);
                    // ILU(0) drops updates outside the stored pattern. Adding
                    // fill here would change the preconditioner and its storage bound.
                    if (targetEntry < rowOffsets[row + 1] && columnIndices[targetEntry] == column)
                    {
                        values[targetEntry] -= values[entry] * values[upperEntry];
                        RequireFinite(values[targetEntry]);
                    }
                }
            }

            // No row swaps, shifts, or small-pivot tolerance are applied. A
            // nonzero tiny pivot is allowed, but any ensuing overflow is explicit.
            if (values[diagonalEntries[row]] == 0)
            {
                throw new ArithmeticException("Zero ILU0 pivot.");
            }
        }
    }

    private static void RequireFinite(double value)
    {
        if (!double.IsFinite(value))
        {
            throw new ArithmeticException("Nonfinite ILU0 arithmetic.");
        }
    }

    // Return the lower-bound position, which can be the row's end when absent.
    private int FindEntry(int row, int column)
    {
        int low = rowOffsets[row];
        int high = rowOffsets[row + 1];
        while (low < high)
        {
            int middle = low + (high - low) / 2;
            if (columnIndices[middle] < column)
            {
                low = middle + 1;
            }
            else
            {
                high = middle;
            }
        }
        return low;
    }

    public double[] Apply(double[] b)
    {
        ArgumentNullException.ThrowIfNull(b);
        if (b.Length != Size || b.Any(value => !double.IsFinite(value)))
        {
            throw new ArgumentException("Invalid ILU0 vector.");
        }
        var x = (double[])b.Clone();

        // Forward substitution: L y = b, with L's implicit unit diagonal.
        for (int row = 0; row < Size; row++)
        {
            for (int entry = rowOffsets[row]; entry < diagonalEntries[row]; entry++)
            {
                x[row] -= values[entry] * x[columnIndices[entry]];
            }
            RequireFinite(x[row]);
        }

        // Back substitution: U x = y. Report finite-range failure rather than
        // returning a nonfinite vector for use by the iterative solver.
        for (int row = Size - 1; row >= 0; row--)
        {
            for (int entry = diagonalEntries[row] + 1; entry < rowOffsets[row + 1]; entry++)
            {
                x[row] -= values[entry] * x[columnIndices[entry]];
            }
            x[row] /= values[diagonalEntries[row]];
            RequireFinite(x[row]);
        }
        return x;
    }
}
#endif
