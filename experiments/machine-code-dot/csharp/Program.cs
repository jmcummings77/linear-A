using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json;
internal static class Program
{
    private static double sink;
    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static double Dot(double[] a, double[] b, int offset, int n)
    {
        var sum = 0.0;
        for (var i = offset; i < offset + n; i++) sum += a[i] * b[i];
        return sum;
    }
    public static void Main(string[] args)
    {
        var n = int.Parse(args[0], CultureInfo.InvariantCulture);
        var iterations = int.Parse(args[1], CultureInfo.InvariantCulture);
        var seed = int.Parse(args[2], CultureInfo.InvariantCulture);
        if (n < 0 || n > 1048576 || iterations < 1 || iterations > 100000000 || seed < 0 || seed > 1000000)
            throw new ArgumentException("out of range");
        var count = Math.Max(n, 1) * 32;
        var a = new double[count]; var b = new double[count];
        for (var i = 0; i < count; i++)
        {
            a[i] = ((i * 17 + seed * 13) % 101 - 50) / 16.0;
            b[i] = ((i * 29 + (seed + 1) * 7) % 103 - 51) / 16.0;
        }
        for (var i = 0; i < Math.Max(8, Math.Min(iterations, 128)); i++) sink = Dot(a, b, (i % 32) * n, n);
        var checksum = 0.0;
        var start = Stopwatch.GetTimestamp();
        for (var i = 0; i < iterations; i++) checksum += Dot(a, b, (i % 32) * n, n);
        var elapsed = Stopwatch.GetElapsedTime(start).TotalNanoseconds;
        if (!double.IsFinite(checksum)) throw new ArithmeticException("nonfinite result");
        Console.WriteLine(JsonSerializer.Serialize(new { elapsed_ns = elapsed, iterations, checksum }));
    }
}
