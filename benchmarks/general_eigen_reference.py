"""Independent analytic spectra and complex right-eigenpair residual checks."""
from fractions import Fraction
import math
from reference import finite_number, matrix


def generated_general_eigen(size, seed):
    """Nonnormal upper block-triangular input with real 2x2 rotation blocks."""
    def entry(i, j):
        k = i // 2
        a, b = 1 + (seed % 17)/16 + k/8, .5 + k/16
        if i == j:
            return a
        if i % 2 == 0 and j == i+1:
            return -b
        if i % 2 == 1 and j == i-1:
            return b
        return ((i*3+j*5+seed) % 11-5)/32 if i < j else 0
    return matrix(size, size, [entry(i, j) for i in range(size) for j in range(size)])


def general_eigen_spectrum(size, seed):
    return [[1+(seed % 17)/16+(i//2)/8,
             0 if i == size-1 and size % 2 else (.5+(i//2)/16)*(-1 if i % 2 == 0 else 1)]
            for i in range(size)]


def general_eigen_checksum(size, seed):
    return math.fsum((i+1)*(re+abs(im)) for i, (re,im) in enumerate(general_eigen_spectrum(size,seed)))+size


def general_eigen_fixtures():
    cases = []
    def add(name, values, size, spectrum, **extra):
        cases.append(dict(name=name, op="eigen_general", a=matrix(size,size,values), b=None,
                          invalid=False, expected_complex_eigenvalues=spectrum, **extra))
    add("empty general eigensystem", [], 0, [])
    add("scalar general eigensystem", [-7.25], 1, [[-7.25,0]])
    add("zero general eigensystem", [0]*9, 3, [[0,0]]*3)
    add("quarter turn complex pair", [0,-1,1,0], 2, [[0,-1],[0,1]])
    add("nonsymmetric real pair", [1,2,3,4], 2, [[(5-math.sqrt(33))/2,0],[(5+math.sqrt(33))/2,0]])
    add("defective Jordan block", [2,1,0,2], 2, [[2,0],[2,0]])
    add("nilpotent Jordan block", [0,1,0,0,0,1,0,0,0], 3, [[0,0]]*3)
    add("repeated diagonal general eigensystem", [2,0,0,0,-1,0,0,0,2], 3, [[-1,0],[2,0],[2,0]])
    add("mixed real and complex spectrum", [1,-2,3,2,1,4,0,0,-3], 3, [[-3,0],[1,-2],[1,2]])
    add("lower triangular real spectrum", [3,0,0,2,-1,0,1,4,2], 3, [[-1,0],[2,0],[3,0]])
    for name, scale in (("tiny",1e-200),("huge",1e200)):
        add(name+" general complex pair", [scale,-2*scale,2*scale,scale], 2, [[scale,-2*scale],[scale,2*scale]])
    add("mixed extreme diagonal general spectrum", [1e300,0,0,1e-300], 2, [[1e-300,0],[1e300,0]])
    add("balanced extreme complex pair", [0,1e300,-1e-300,0], 2, [[0,-1],[0,1]], spectrum_scale=1, componentwise=True)
    add("subnormal general scalar", [5e-324], 1, [[5e-324,0]])
    # Exact rational Householder similarity fills the matrix without an eigensolver.
    source = generated_general_eigen(4,3)
    q = [[Fraction(int(i==j))-Fraction(1,2) for j in range(4)] for i in range(4)]
    dense = [float(sum(q[i][k]*Fraction(source['values'][k*4+l])*q[l][j]
                       for k in range(4) for l in range(4))) for i in range(4) for j in range(4)]
    add("dense real similarity with complex spectrum",dense,4,general_eigen_spectrum(4,3))
    for size,seed in ((3,0),(8,9),(16,17)):
        add("analytic general block spectrum %d" % size,generated_general_eigen(size,seed)['values'],size,
            general_eigen_spectrum(size,seed))
    cases.append(dict(name="nonsquare general eigensystem",op="eigen_general",a=matrix(2,3,[1,2,3,4,5,6]),b=None,invalid=True))
    return cases


def assert_general_eigen_result(actual, a, expected=None, tolerance=1e-9, spectrum_scale=None, componentwise=False):
    n=a['rows']
    keys={'eigenvalues_real','eigenvalues_imag','eigenvectors_real','eigenvectors_imag'}
    if a['cols'] != n or not isinstance(actual,dict) or set(actual)!=keys:
        raise AssertionError('Expected a square complex eigensystem with split real/imaginary parts')
    arrays=[]
    for key in ('eigenvalues_real','eigenvalues_imag'):
        value=actual[key]
        if not isinstance(value,list) or len(value)!=n or not all(map(finite_number,value)):
            raise AssertionError('Eigenvalues must be finite arrays of length n')
        arrays.append(value)
    for key in ('eigenvectors_real','eigenvectors_imag'):
        value=actual[key]
        if (not isinstance(value,dict) or set(value)!={'rows','cols','values'}
                or type(value['rows']) is not int or type(value['cols']) is not int
                or value['rows']!=n or value['cols']!=n or not isinstance(value['values'],list)
                or len(value['values'])!=n*n or not all(map(finite_number,value['values']))):
            raise AssertionError('Eigenvectors must be finite n-by-n matrices')
        arrays.append(value['values'])
    vr,vi,qr,qi=arrays
    spectrum=[complex(r,i) for r,i in zip(vr,vi)]
    if any((vr[i],vi[i])>(vr[i+1],vi[i+1]) for i in range(n-1)):
        raise AssertionError('General eigenvalues must be lexicographically sorted')
    scale=max([abs(x) for x in a['values']]+[abs(x) for x in vr+vi]+[0]) or 1
    norm=math.sqrt(math.fsum((x/scale)**2 for x in a['values']))
    residual=0
    for col in range(n):
        vector=[complex(qr[i*n+col],qi[i*n+col]) for i in range(n)]
        if abs(math.fsum(abs(v)**2 for v in vector)-1)>tolerance*max(n,1):
            raise AssertionError('Complex eigenvectors must be unit columns')
        errors=[abs(sum((a['values'][i*n+k]/scale)*vector[k] for k in range(n))
                    - complex(vr[col]/scale,vi[col]/scale)*vector[i]) for i in range(n)]
        if componentwise:
            for row in range(n):
                terms=[a['values'][row*n+k]*vector[k] for k in range(n)]
                target=spectrum[col]*vector[row]
                error=abs(sum(terms)-target)
                bound=tolerance*(sum(abs(term) for term in terms)+abs(target))
                if error>bound:
                    raise AssertionError('Componentwise eigenpair residual exceeds tolerance')
        column_residual=math.hypot(*errors)
        if column_residual>tolerance*max(norm,abs(spectrum[col]/scale)):
            raise AssertionError('Complex eigenpair residual exceeds tolerance')
        residual=max(residual,column_residual/(norm or 1))
    if expected is not None:
        if len(expected)!=n:
            raise AssertionError('Expected spectrum must contain n eigenvalues')
        # Match multiplicities independently of rounding within repeated real parts.
        remaining=[complex(r,i) for r,i in expected]
        comparison_scale=spectrum_scale if spectrum_scale is not None else scale
        for value in spectrum:
            closest=min(range(len(remaining)),key=lambda j: abs(value/comparison_scale-remaining[j]/comparison_scale))
            error=abs(value/comparison_scale-remaining.pop(closest)/comparison_scale)
            if error>tolerance*(max(norm,1) if spectrum_scale is None else 1):
                raise AssertionError('General eigenvalue differs from the independent spectrum')
    return {'relative_residual':residual}
