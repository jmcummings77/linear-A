namespace linear_A;

/// <summary>Selects a determinant algorithm for <see cref="Matrix{T}"/>.</summary>
public enum DeterminantAlgorithm
{
    /// <summary>Uses Bareiss for built-in integers and decimal, LU for binary floats, and cofactor expansion for other types.</summary>
    Auto,

    /// <summary>Uses checked, division-free cofactor expansion. Its factorial cost is suitable only for small matrices.</summary>
    Cofactor,

    /// <summary>Uses partial-pivoting elimination with double working precision. Supports double, float, and Half.</summary>
    Lu,

    /// <summary>Uses exact BigInteger elimination. Supports built-in integer types and decimal, with one final conversion.</summary>
    Bareiss,

    /// <summary>Uses Cholesky factorization for finite, exactly symmetric positive-definite double, float, or Half matrices.</summary>
    Cholesky,
}
