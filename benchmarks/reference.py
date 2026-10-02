"""Independent small-matrix oracles and shared workloads (standard library only)."""
from fractions import Fraction
from itertools import permutations
import math
import random

OPERATIONS = ("add", "subtract", "scale", "transpose", "multiply", "trace", "determinant")


def matrix(rows, cols, values):
    assert len(values) == rows * cols
    return {"rows": rows, "cols": cols, "values": list(values)}


def determinant(a):
    """Leibniz formula: independent of the elimination algorithms in the ports."""
    n = a["rows"]
    assert n == a["cols"] and n <= 5
    total = Fraction(0)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = Fraction(-1 if inversions % 2 else 1)
        for i, j in enumerate(p):
            term *= Fraction(a["values"][i * n + j])
        total += term
    return float(total)


def expected(op, a, b=None, scalar=1.25):
    r, c, values = a["rows"], a["cols"], a["values"]
    if op == "add":
        return matrix(r, c, [x + y for x, y in zip(values, b["values"])])
    if op == "subtract":
        return matrix(r, c, [x - y for x, y in zip(values, b["values"])])
    if op == "scale":
        return matrix(r, c, [x * scalar for x in values])
    if op == "transpose":
        return matrix(c, r, [values[i * c + j] for j in range(c) for i in range(r)])
    if op == "multiply":
        return matrix(r, b["cols"], [sum(Fraction(values[i*c+k]) * Fraction(b["values"][k*b["cols"]+j])
               for k in range(c)) for i in range(r) for j in range(b["cols"])])
    if op == "trace":
        return {"value": sum(values[i*c+i] for i in range(r))}
    if op == "determinant":
        return {"value": determinant(a)}
    if op == "triangular":
        return {"upper": r == c and all(values[i*c+j] == 0 for i in range(r) for j in range(i)),
                "lower": r == c and all(values[i*c+j] == 0 for i in range(r) for j in range(i+1,c))}
    raise ValueError(op)


def fixtures():
    cases = []
    def add(name, op, a, b=None, scalar=1.25, invalid=False):
        case = {"name": name, "op": op, "a": a, "b": b, "scalar": scalar, "invalid": invalid}
        if not invalid:
            case["expected"] = expected(op, a, b, scalar)
        cases.append(case)
    a = matrix(2, 3, [0.5, -2, 3, 4.25, 5, -6])
    b = matrix(2, 3, [1, 2.5, -1, -4, 0, 6.25])
    for op in ("add", "subtract"):
        add("fractional rectangular " + op, op, a, b)
        add("invalid shape " + op, op, a, matrix(3, 2, b["values"]), invalid=True)
    add("fractional scale", "scale", a, scalar=-0.25)
    add("wide transpose", "transpose", a)
    add("tall transpose", "transpose", matrix(3, 2, a["values"]))
    add("rectangular multiply", "multiply", a, matrix(3, 2, [4, -2, 0.5, 5, 3, 1.25]))
    add("invalid multiply", "multiply", a, b, invalid=True)
    square = matrix(3, 3, [6, 1, 1, 4, -2, 5, 2, 8, 7])
    for op in ("trace", "determinant"):
        add("dense " + op, op, square)
        add("nonsquare " + op, op, a, invalid=True)
        add("empty " + op, op, matrix(0, 0, []))
        add("scalar " + op, op, matrix(1, 1, [-0.5]))
    add("determinant requires pivot", "determinant", matrix(3, 3, [0,2,1,3,1,4,2,5,6]))
    add("singular determinant", "determinant", matrix(3,3,[1,2,3,2,4,6,4,5,6]))
    add("dense 4x4 determinant", "determinant", matrix(4,4,[3,2,0,1,4,0,1,2,3,0,2,1,9,2,3,1]))
    add("fractional determinant", "determinant", matrix(2,2,[0.5,1.25,-2,0.75]))
    rng = random.Random(601)
    for index in range(8):
        size = 2 + index % 3
        add("seeded determinant %d" % index, "determinant", matrix(size,size,[rng.randrange(-8,9)/4 for _ in range(size*size)]))
    for name, mat in [("upper",matrix(3,3,[1,2,3,0,4,5,0,0,6])),
                      ("lower",matrix(3,3,[1,0,0,2,4,0,3,5,6])),
                      ("dense",square),("rectangle",a),("empty",matrix(0,0,[])),
                      ("tiny off diagonal",matrix(2,2,[1,5e-324,5e-324,1]))]:
        add(name + " triangular", "triangular", mat)
    add("empty transpose", "transpose", matrix(2,0,[]))
    add("empty inner product", "multiply", matrix(2,0,[]), matrix(0,3,[]))
    add("empty outer product", "multiply", matrix(0,2,[]), matrix(2,3,[1,2,3,4,5,6]))
    add("empty addition", "add", matrix(0,3,[]), matrix(0,3,[]))
    return cases


def generated(size, seed, diagonal=False):
    values = [((index*17+seed*13)%101-50)/16.0 for index in range(size*size)]
    if diagonal:
        for i in range(size):
            values[i*size+i] += size*4
    return matrix(size,size,values)


def generated_spd(size, seed):
    """Exact dyadic, symmetric, strictly diagonally dominant positive definite input."""
    a = generated(size, seed)
    values = [(a["values"][i*size+j] + a["values"][j*size+i])/2
              + (4*size if i == j else 0) for i in range(size) for j in range(size)]
    return matrix(size, size, values)


def benchmark_determinant(a, denominator=16):
    """Exact fraction-free elimination on dyadic integer numerators."""
    n = a["rows"]
    if not n:
        return 1.0
    scaled = [Fraction(value) * denominator for value in a["values"]]
    if any(value.denominator != 1 for value in scaled):
        raise ValueError("Input is not exact at the requested denominator")
    rows = [[int(scaled[i*n+j]) for j in range(n)] for i in range(n)]
    previous, sign = 1, 1
    for k in range(n-1):
        pivot = next((r for r in range(k,n) if rows[r][k]), None)
        if pivot is None:
            return 0.0
        if pivot != k:
            rows[k],rows[pivot] = rows[pivot],rows[k]
            sign = -sign
        value = rows[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = rows[i][j]*value-rows[i][k]*rows[k][j]
                assert numerator % previous == 0
                rows[i][j] = numerator//previous
            rows[i][k] = 0
        previous = value
    return float(Fraction(sign*rows[-1][-1],denominator**n))


def benchmark_checksum(op, size, seed):
    if op in ("cross", "rotation2d", "rotation3d"):
        from vector_reference import vector_checksum
        return vector_checksum(op, size, seed)
    if op == "eigen_general":
        from general_eigen_reference import general_eigen_checksum
        return general_eigen_checksum(size, seed)
    if op == "eigen_symmetric":
        from eigen_reference import eigen_checksum
        return eigen_checksum(size, seed)
    if op in ("determinant_spd", "determinant_spd_lu", "determinant_cholesky"):
        return benchmark_determinant(generated_spd(size, seed), 32)
    if op in ("determinant_lu", "determinant_cofactor", "determinant_small"):
        op = "determinant"
    a, b = generated(size,seed,op=="determinant"),generated(size,seed+1)
    if op == "determinant":
        return benchmark_determinant(a)
    if op == "trace":
        return sum(a["values"][i*size+i] for i in range(size))
    total = 0.0
    for index in (0,size*size//2,size*size-1):
        i,j = divmod(index,size)
        if op == "add": value = a["values"][index]+b["values"][index]
        elif op == "subtract": value = a["values"][index]-b["values"][index]
        elif op == "scale": value = a["values"][index]*1.25
        elif op == "transpose": value = a["values"][j*size+i]
        elif op == "multiply": value = sum(a["values"][i*size+k]*b["values"][k*size+j] for k in range(size))
        else: raise ValueError(op)
        total += value
    return total


def finite_number(value):
    """A finite JSON number, excluding bool (which Python treats as an int)."""
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def assert_result(actual, wanted):
    if not isinstance(actual, dict):
        raise AssertionError("Result must be a JSON object")
    if set(actual) != set(wanted):
        raise AssertionError("Wrong result keys: %r; expected %r" % (list(actual),list(wanted)))
    for key in wanted:
        x,y = actual[key],wanted[key]
        if key in ("rows","cols"):
            if type(x) is not int or x < 0:
                raise AssertionError("%s must be a nonnegative integer" % key)
            if x != y: raise AssertionError("%s: %r != %r" % (key,x,y))
        elif key in ("upper","lower"):
            if type(x) is not bool:
                raise AssertionError("%s must be a boolean" % key)
            if x != y: raise AssertionError("%s: %r != %r" % (key,x,y))
        else:
            if key == "values" and not isinstance(x, list):
                raise AssertionError("values must be an array")
            xs,ys = (x,y) if key=="values" else ([x],[y])
            if len(xs)!=len(ys): raise AssertionError("Wrong number of values")
            for index,(left,right) in enumerate(zip(xs,ys)):
                if not finite_number(left) or not math.isclose(left,float(right),rel_tol=1e-9,abs_tol=1e-10):
                    raise AssertionError("%s[%d]: %r != %r" % (key,index,left,right))
