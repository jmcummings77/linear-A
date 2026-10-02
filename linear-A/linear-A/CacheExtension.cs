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
        ///     Creates a delegate that caches results by input for the lifetime of the delegate.
        /// </summary>
        /// <remarks>Concurrent calls for the same uncached input may evaluate the function more than once.</remarks>
        /// <param name="f"></param>
        /// <typeparam name="T"></typeparam>
        /// <typeparam name="TResult"></typeparam>
        /// <returns></returns>
        public static Func<T, TResult> Memoize<T, TResult>(this Func<T, TResult> f) where T : notnull
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            var cache = new ConcurrentDictionary<T, TResult>();
            return value => cache.GetOrAdd(value, f);
        }
    }
}
