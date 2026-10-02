using System;
#if NET10_0_OR_GREATER
using System.Numerics;
#endif

namespace linear_A;

/// <summary>Creates right-handed active rotations of column vectors. Angles are in radians.</summary>
public static class MatrixRotation
{
    internal static bool IsVector3(int rows, int columns) => (rows == 3 && columns == 1) || (rows == 1 && columns == 3);
    internal static bool IsFinite(double value) => !double.IsNaN(value) && !double.IsInfinity(value);

    /// <summary>Creates a counterclockwise 2D rotation.</summary>
    public static DoubleMatrix Create2D(double radians)
    {
        ValidateAngle(radians);
        var c = Math.Cos(radians);
        var s = Math.Sin(radians);
        return From(new[,] { { c, -s }, { s, c } });
    }

    /// <summary>Creates a right-handed rotation about the positive x axis.</summary>
    public static DoubleMatrix CreateX(double radians) => Rodrigues(1, 0, 0, radians);

    /// <summary>Creates a right-handed rotation about the positive y axis.</summary>
    public static DoubleMatrix CreateY(double radians) => Rodrigues(0, 1, 0, radians);

    /// <summary>Creates a right-handed rotation about the positive z axis.</summary>
    public static DoubleMatrix CreateZ(double radians) => Rodrigues(0, 0, 1, radians);

    /// <summary>Creates a right-handed rotation about a nonzero finite 3D axis.</summary>
    /// <remarks>The axis may be 3-by-1 or 1-by-3, is normalized stably, and is not changed.</remarks>
    public static DoubleMatrix CreateAxisAngle(DoubleMatrix axis, double radians)
    {
        if (axis == null) throw new ArgumentNullException(nameof(axis));
        if (!IsVector3(axis.RowCount, axis.ColumnCount)) throw new ArgumentException("Rotation axis must be 3-by-1 or 1-by-3.", nameof(axis));
        var a = axis.RowCount == 1 ? axis.GetRow(0) : axis.GetColumn(0);
        return NormalizeAndRotate(a, radians);
    }

#if NET10_0_OR_GREATER
    /// <summary>Approximates a generic numeric axis as doubles and creates its rotation.</summary>
    public static DoubleMatrix CreateAxisAngle<T>(Matrix<T> axis, double radians) where T : INumber<T>
    {
        ArgumentNullException.ThrowIfNull(axis);
        if (!IsVector3(axis.RowCount, axis.ColumnCount)) throw new ArgumentException("Rotation axis must be 3-by-1 or 1-by-3.", nameof(axis));
        var a = new double[3];
        for (var i = 0; i < 3; i++) a[i] = double.CreateChecked(axis[axis.RowCount == 1 ? 0 : i, axis.RowCount == 1 ? i : 0]);
        return NormalizeAndRotate(a, radians);
    }
#endif

    private static DoubleMatrix NormalizeAndRotate(double[] axis, double radians)
    {
        var scale = 0d;
        foreach (var value in axis)
        {
            if (!IsFinite(value)) throw new ArgumentException("Rotation axis must be finite.", nameof(axis));
            scale = Math.Max(scale, Math.Abs(value));
        }
        if (scale == 0) throw new ArgumentException("Rotation axis must be nonzero.", nameof(axis));
        var x = axis[0] / scale;
        var y = axis[1] / scale;
        var z = axis[2] / scale;
        var norm = Math.Sqrt(x * x + y * y + z * z);
        return Rodrigues(x / norm, y / norm, z / norm, radians);
    }

    private static void ValidateAngle(double radians)
    {
        if (!IsFinite(radians)) throw new ArgumentOutOfRangeException(nameof(radians), "Rotation angle must be finite.");
    }

    private static DoubleMatrix Rodrigues(double x, double y, double z, double radians)
    {
        ValidateAngle(radians);
        var c = Math.Cos(radians);
        var s = Math.Sin(radians);
        var halfSine = Math.Sin(radians / 2);
        var t = Math.Abs(radians) < 1 ? 2 * halfSine * halfSine : 1 - c;
        return From(new[,] {
            { c + x*x*t, x*y*t - z*s, x*z*t + y*s },
            { y*x*t + z*s, c + y*y*t, y*z*t - x*s },
            { z*x*t - y*s, z*y*t + x*s, c + z*z*t }
        });
    }

    private static DoubleMatrix From(double[,] values)
    {
        var result = new DoubleMatrix(values.GetLength(0), values.GetLength(1));
        for (var row = 0; row < result.RowCount; row++)
            for (var column = 0; column < result.ColumnCount; column++) result[row, column] = values[row, column];
        return result;
    }
}
