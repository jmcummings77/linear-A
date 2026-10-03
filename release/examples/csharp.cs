using System;
using linear_A;
var a = new Matrix<double>(new double[,] {{4,1},{2,3}});
var b = new Matrix<double>(new double[,] {{6},{8}});
var x = a.Solve(b);
var r = a * x;
for (var i=0;i<2;i++)
    if (!double.IsFinite(x[i,0]) || Math.Abs(x[i,0]-(i+1))>1e-12
        || !double.IsFinite(r[i,0]) || Math.Abs(r[i,0]-b[i,0])>1e-12) throw new Exception("incorrect solution");
Console.WriteLine("solution: 1, 2");

var decomposition = a.Svd();
if (decomposition.Values.Length != 2 || decomposition.Values[1] <= 0) throw new Exception("incorrect SVD");
