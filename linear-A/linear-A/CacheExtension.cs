#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="CacheExtensions.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   A class for caching the results of functions for dynamic programming.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A
{
    using System;
    using System.Collections.Concurrent;

    /// <summary>
    ///     A class for caching the results of functions for dynamic programming
    /// </summary>
    public static class CacheExtension
    {
        /// <summary>
        ///     A delegate for memoizing function results
        /// </summary>
        /// <param name="f"></param>
        /// <typeparam name="T"></typeparam>
        /// <typeparam name="TResult"></typeparam>
        /// <returns></returns>
        public static Func<T, TResult> Memoize<T, TResult>(this Func<T, TResult> f) => a => new ConcurrentDictionary<T, TResult>().GetOrAdd(a, f);
    }
}