"""JSON check/benchmark CLI; see benchmarks/PROTOCOL.md."""

import json
import sys
from time import perf_counter_ns
from matrix import Matrix, SymmetricEigenDecomposition, GeneralEigenDecomposition, finite


DETERMINANTS = {"determinant", "determinant_lu", "determinant_cholesky", "determinant_spd_lu"}
OPERATIONS = {"add", "subtract", "scale", "transpose", "multiply", "trace", "triangular", "eigen_symmetric", "eigen_general",
              "solve", "solve_cholesky", "least_squares", "rcond", "cross", "rotation2d", "rotation3d"} | DETERMINANTS


def integer(text, minimum=0, maximum=None):
    value = int(text)
    if value < minimum or maximum is not None and value > maximum:
        raise ValueError("integer argument is outside the allowed range")
    return value


def operation(name, a, b=None, scalar=1.25):
    if name == "add":
        return a.add(b)
    if name == "subtract":
        return a.subtract(b)
    if name == "multiply":
        return a.multiply(b)
    if name == "solve":
        return a.solve(b)
    if name == "solve_cholesky":
        return a.factor_cholesky().solve(b)
    if name == "least_squares":
        return a.least_squares(b)
    if name == "rcond":
        return a.factor_lu().reciprocal_condition()
    if name == "cross":
        return a.cross(b)
    if name == "rotation2d":
        if (a.rows, a.cols) != (0, 0):
            raise ValueError("2D rotation construction expects a 0 by 0 input")
        return Matrix.rotation2d(scalar)
    if name == "rotation3d":
        return Matrix.rotation_axis_angle(a, scalar)
    if name == "scale":
        return a.scale(scalar)
    if name == "transpose":
        return a.transpose()
    if name == "trace":
        return a.trace()
    if name == "determinant":
        return a.determinant()
    if name in {"determinant_lu", "determinant_spd_lu"}:
        return a.determinant("lu")
    if name == "determinant_cholesky":
        return a.determinant("cholesky")
    if name == "triangular":
        return a.triangular()
    if name == "eigen_symmetric":
        return a.eigen_symmetric()
    if name == "eigen_general":
        return a.eigen_general()
    raise ValueError("unsupported operation")


def result_json(result):
    if isinstance(result, GeneralEigenDecomposition):
        return {"eigenvalues_real":result.values_real,"eigenvalues_imag":result.values_imag,
                "eigenvectors_real":result_json(result.vectors_real),"eigenvectors_imag":result_json(result.vectors_imag)}
    if isinstance(result, SymmetricEigenDecomposition):
        return {"eigenvalues": result.values, "eigenvectors": result_json(result.vectors)}
    if isinstance(result, Matrix):
        return {"rows": result.rows, "cols": result.cols, "values": result.values}
    if isinstance(result, tuple):
        return {"upper": result[0], "lower": result[1]}
    return {"value": finite(result)}


def check(args):
    if len(args) < 3 or args[0] not in OPERATIONS:
        raise ValueError("usage: check OP ROWS COLS [BROWS BCOLS | SCALAR]")
    name, rows, cols = args[0], integer(args[1]), integer(args[2])
    binary = name in {"add", "subtract", "multiply", "cross", "solve", "solve_cholesky", "least_squares"}
    scalar_operation = name in {"scale", "rotation2d", "rotation3d"}
    expected_args = 5 if binary else 4 if scalar_operation else 3
    if len(args) != expected_args:
        raise ValueError("incorrect number of operation arguments")
    brows, bcols = (integer(args[3]), integer(args[4])) if binary else (0, 0)
    scalar = finite(args[3]) if scalar_operation else 1.25
    values = [finite(token) for token in sys.stdin.read().split()]
    count = rows * cols
    if len(values) != count + brows * bcols:
        raise ValueError("input value count does not match matrix dimensions")
    a = Matrix(rows, cols, values[:count])
    b = Matrix(brows, bcols, values[count:]) if binary else None
    return result_json(operation(name, a, b, scalar))


def benchmark(args):
    if len(args) != 4 or args[0] not in OPERATIONS - {"triangular"}:
        raise ValueError("usage: bench OP SIZE ITERATIONS SEED")
    name, size, iterations, seed = args[0], integer(args[1], 1), integer(args[2], 1), integer(args[3], 0, 2147483646)
    required_size = {"cross":3, "rotation2d":2, "rotation3d":3}.get(name)
    if required_size is not None and size != required_size:
        raise ValueError("invalid size for vector or rotation benchmark")
    spd = name in {"determinant_cholesky", "determinant_spd_lu"}
    numerators = [(i * 17 + seed * 13) % 101 - 50 for i in range(size * size)]
    values = [(numerators[row * size + col] + numerators[col * size + row]) / 32.0
              if spd else numerators[row * size + col] / 16.0
              for row in range(size) for col in range(size)]
    a = Matrix(size, size, values)
    b = Matrix(size, size, [((i * 17 + (seed + 1) * 13) % 101 - 50) / 16.0 for i in range(size * size)])
    if name in DETERMINANTS:
        for i in range(size):
            a[i, i] += size * 4
    if name == "eigen_symmetric":
        a = Matrix(size, size, [2+(seed % 17)/16 if row == col else -1 if abs(row-col) == 1 else 0
                                for row in range(size) for col in range(size)])
    if name == "eigen_general":
        def entry(i,j):
            k=i//2
            if i==j: return 1+(seed%17)/16+k/8
            if i%2==0 and j==i+1: return -(.5+k/16)
            if i%2==1 and j==i-1: return .5+k/16
            return ((i*3+j*5+seed)%11-5)/32 if i<j else 0
        a=Matrix(size,size,[entry(i,j) for i in range(size) for j in range(size)])
    scalar = .5 if name in {"rotation2d", "rotation3d"} else 1.25
    if name == "cross":
        a, b = Matrix(3, 1, a.values[:3]), Matrix(3, 1, b.values[:3])
        b[2,0] = -b[2,0]
    elif name == "rotation2d":
        a = Matrix()
    elif name == "rotation3d":
        a = Matrix(3, 1, [1,2,3])
    for _ in range(max(5, min(iterations, 100))):
        operation(name, a, b, scalar)
    checksum = 0.0
    start = perf_counter_ns()
    for _ in range(iterations):
        result = operation(name, a, b, scalar)
        if isinstance(result, GeneralEigenDecomposition):
            checksum += sum((i+1)*(r+abs(im)) for i,(r,im) in enumerate(zip(result.values_real,result.values_imag)))
            checksum += sum(x*x for x in result.vectors_real.values+result.vectors_imag.values)
        elif isinstance(result, SymmetricEigenDecomposition):
            checksum += sum((index+1)*value for index, value in enumerate(result.values))
            checksum += sum(value*value for value in result.vectors.values)
        else:
            checksum += result.checksum() if isinstance(result, Matrix) else result
        del result
    elapsed = perf_counter_ns() - start
    return {"elapsed_ns": elapsed, "iterations": iterations, "checksum": finite(checksum)}


def main(args=None):
    args = sys.argv[1:] if args is None else args
    if not args or args[0] not in {"check", "bench"}:
        raise ValueError("usage: runner.py check|bench ...")
    result = check(args[1:]) if args[0] == "check" else benchmark(args[1:])
    print(json.dumps(result, allow_nan=False, separators=(",", ":")))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, TypeError, IndexError, OverflowError, MemoryError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
