using System;
using NUnit.Framework;

namespace linear_A.Test;

public class CacheExtensionTests
{
    [Test]
    public void MemoizeRetainsResultsForEachInput()
    {
        var calls = 0;
        Func<int, int> square = value =>
        {
            calls++;
            return value * value;
        };
        var cached = square.Memoize();

        Assert.That(cached(3), Is.EqualTo(9));
        Assert.That(cached(3), Is.EqualTo(9));
        Assert.That(cached(4), Is.EqualTo(16));
        Assert.That(cached(4), Is.EqualTo(16));
        Assert.That(calls, Is.EqualTo(2));
    }

    [Test]
    public void MemoizeKeepsSeparateDelegatesIndependent()
    {
        Func<int, int> increment = value => value + 1;
        Func<int, int> decrement = value => value - 1;

        Assert.That(increment.Memoize()(2), Is.EqualTo(3));
        Assert.That(decrement.Memoize()(2), Is.EqualTo(1));
    }

    [Test]
    public void MemoizeRejectsNullFunction()
    {
        Func<int, int> function = null!;
        Assert.Throws<ArgumentNullException>(() => function.Memoize());
    }
}
