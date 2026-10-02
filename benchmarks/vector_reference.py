"""Independent vector-cross and quaternion rotation references for the shared runners."""
from fractions import Fraction
import math

from reference import matrix


def _vector(a):
    if (a["rows"], a["cols"]) not in ((3, 1), (1, 3)) or len(a["values"]) != 3:
        raise ValueError("Expected a three-component row or column vector")
    if not all(math.isfinite(value) for value in a["values"]):
        raise ValueError("Expected finite vector entries")
    return a["values"]


def _quaternion_product(a, b):
    w, x, y, z = a
    r, i, j, k = b
    return [w*r-x*i-y*j-z*k, w*i+x*r+y*k-z*j,
            w*j-x*k+y*r+z*i, w*k+x*j-y*i+z*r]


def _quaternion_rotation(axis, angle):
    if not math.isfinite(angle):
        raise ValueError("Expected a finite angle")
    largest = max(map(abs, axis))
    if not largest:
        raise ValueError("Rotation axis must be nonzero")
    scaled = [value/largest for value in axis]
    norm = math.hypot(*scaled)
    sine = math.sin(angle/2)
    q = [math.cos(angle/2)] + [sine*value/norm for value in scaled]
    conjugate = [q[0], -q[1], -q[2], -q[3]]
    columns = []
    for index in range(3):
        basis = [0.0] + [float(component == index) for component in range(3)]
        columns.append(_quaternion_product(_quaternion_product(q, basis), conjugate)[1:])
    return [columns[col][row] for row in range(3) for col in range(3)]


def vector_expected(op, a, b=None, scalar=.5):
    if op == "cross":
        left, right = _vector(a), _vector(b)
        values = [float(Fraction(left[(i+1) % 3])*Fraction(right[(i+2) % 3])
                        - Fraction(left[(i+2) % 3])*Fraction(right[(i+1) % 3])) for i in range(3)]
        if not all(map(math.isfinite, values)):
            raise ValueError("Cross product is outside the finite float64 range")
        return matrix(a["rows"], a["cols"], values)
    if op == "rotation2d":
        if (a["rows"], a["cols"]) != (0, 0):
            raise ValueError("2D rotation construction expects an empty input")
        rotation = _quaternion_rotation([0, 0, 1], scalar)
        return matrix(2, 2, [rotation[0], rotation[1], rotation[3], rotation[4]])
    if op == "rotation3d":
        return matrix(3, 3, _quaternion_rotation(_vector(a), scalar))
    raise ValueError("Unknown vector operation: " + str(op))


def vector_inputs(op, size, seed):
    """Inputs used by every timed vector/rotation workload; scalar angles are radians."""
    expected_size = {"cross": 3, "rotation2d": 2, "rotation3d": 3}.get(op)
    if type(size) is not int or expected_size is None or size != expected_size:
        raise ValueError("Invalid vector benchmark operation or size")
    if type(seed) is not int or not 0 <= seed <= 2147483646:
        raise ValueError("Invalid benchmark seed")
    if op == "cross":
        left = [((index*17+seed*13) % 101-50)/16 for index in range(3)]
        right = [((index*17+(seed+1)*13) % 101-50)/16 for index in range(3)]
        right[2] = -right[2]
        return matrix(3, 1, left), matrix(3, 1, right), .5
    if op == "rotation2d":
        return matrix(0, 0, []), None, .5
    return matrix(3, 1, [1, 2, 3]), None, .5


def vector_checksum(op, size, seed):
    a, b, scalar = vector_inputs(op, size, seed)
    values = vector_expected(op, a, b, scalar)["values"]
    return values[0] + values[len(values)//2] + values[-1]


def vector_fixtures():
    cases = []

    def add(name, op, a, b=None, scalar=.5, invalid=False):
        case = {"name": name, "op": op, "a": a, "b": b, "scalar": scalar, "invalid": invalid}
        if not invalid:
            case["expected"] = vector_expected(op, a, b, scalar)
        cases.append(case)

    column = lambda values: matrix(3, 1, values)
    row = lambda values: matrix(1, 3, values)
    add("cross right handed basis", "cross", column([1,0,0]), column([0,1,0]))
    add("cross reversed row basis", "cross", row([0,1,0]), row([1,0,0]))
    add("cross mixed row and column", "cross", row([1,2,3]), column([-2,4,1]))
    add("cross mixed column and row", "cross", column([1,2,3]), row([-2,4,1]))
    add("cross fractional vectors", "cross", column([.5,-1.25,2.5]), column([3,.25,-2]))
    add("cross parallel vectors", "cross", column([1,2,3]), column([2,4,6]))
    add("cross zero vector", "cross", row([0,0,0]), column([-2,4,1]))
    add("cross mixed magnitudes", "cross", column([1e300,0,0]), column([0,1e-300,0]))
    add("cross large finite result", "cross", column([1e100,0,0]), column([0,1e100,0]))
    add("cross empty vector rejected", "cross", matrix(0,0,[]), column([1,0,0]), invalid=True)
    add("cross short vector rejected", "cross", matrix(2,1,[1,2]), column([1,0,0]), invalid=True)
    add("cross matrix rejected", "cross", matrix(2,2,[1,2,3,4]), column([1,0,0]), invalid=True)
    add("cross invalid right vector", "cross", column([1,0,0]), matrix(1,2,[1,2]), invalid=True)
    add("cross nonfinite result rejected", "cross", column([1e308,0,0]), column([0,1e308,0]), invalid=True)

    empty = matrix(0, 0, [])
    for name, angle in (("identity", 0), ("quarter turn", math.pi/2), ("clockwise quarter turn", -math.pi/2),
                        ("fractional angle", .5), ("tiny angle", 1e-10)):
        add("2D rotation " + name, "rotation2d", empty, scalar=angle)
    add("2D rotation rejects input matrix", "rotation2d", matrix(1,1,[1]), scalar=.5, invalid=True)
    add("2D rotation rejects empty rectangle", "rotation2d", matrix(0,2,[]), scalar=.5, invalid=True)

    for name, axis, angle in (("X quarter turn", [1,0,0], math.pi/2),
                              ("Y quarter turn", [0,1,0], math.pi/2),
                              ("Z clockwise quarter turn", [0,0,1], -math.pi/2),
                              ("arbitrary axis", [1,2,3], .5),
                              ("signed axis", [-1,2,-3], .5),
                              ("zero angle", [1,2,3], 0),
                              ("tiny angle", [1,2,3], 1e-10),
                              ("huge axis", [1e300,2e300,3e300], .5),
                              ("tiny axis", [1e-300,2e-300,3e-300], .5),
                              ("overflowing unscaled norm", [1.7e308,1.7e308,1.7e308], .5),
                              ("subnormal axis", [5e-324,0,0], math.pi/2)):
        add("3D rotation " + name, "rotation3d", column(axis), scalar=angle)
    add("3D rotation row axis", "rotation3d", row([1,2,3]), scalar=.5)
    for name, axis in (("zero axis", column([0,0,0])), ("empty axis", empty),
                       ("short axis", matrix(2,1,[1,2])), ("matrix axis", matrix(2,2,[1,2,3,4]))):
        add("3D rotation rejects " + name, "rotation3d", axis, invalid=True)
    return cases
