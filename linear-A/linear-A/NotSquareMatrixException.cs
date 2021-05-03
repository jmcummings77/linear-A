#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="NotSquareMatrixException.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   A custom exception thrown when an operation that is only valid on a square matrix is attempted on a non-square matrix
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A
{
    using System;

    /// <inheritdoc />
    /// <summary>
    ///     Exception thrown when an operation that is only valid on a square matrix is attempted on a non-square matrix
    /// </summary>
    public class NotSquareMatrixException : Exception
    {
    }
}