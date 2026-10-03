#!/usr/bin/env python3
"""Seeded numerical properties, differential checks, replay and bounded reduction."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import time

import run as harness
from reference import matrix, expected, assert_result
from solve_reference import exact_solve, exact_least_squares
from eigen_reference import assert_eigen_result
from publication import PublicSanitizer

ROOT = harness.ROOT
SANITIZER = PublicSanitizer(root=ROOT)
VERSION = 1


def plain(value):
    if isinstance(value, dict): return {k: plain(v) for k, v in value.items()}
    if isinstance(value, list): return [plain(v) for v in value]
    if isinstance(value, str): return SANITIZER.text(value)
    return value


def generated(seed):
    """Small dyadic families: bounded products, exact cancellation, controlled rank."""
    rng = random.Random(seed)
    def dense(r, c): return matrix(r, c, [rng.randrange(-16, 17)/4 for _ in range(r*c)])
    def case(family, op, a, b=None):
        return dict(generator_version=VERSION, seed=seed, family=family, op=op, a=a, b=b)
    m, k, n = [rng.randrange(0, 6) for _ in range(3)]
    a, b = dense(m, k), dense(k, n)
    yield case('rectangular and empty products', 'multiply', a, b)
    yield case('transpose involution', 'transpose_twice', a)
    size = 1+rng.randrange(4)
    a = dense(size, size)
    yield case('dense determinant and transpose identity', 'determinant_transpose', a)
    singular = copy.deepcopy(a)
    if size > 1: singular['values'][size:2*size] = singular['values'][:size]
    else: singular['values'][0] = 0
    yield case('duplicate row determinant', 'determinant', singular)
    b = matrix(size, size, [-v for v in a['values']])
    yield case('exact cancellation', 'add', a, b)
    yield case('product transpose identity', 'product_transpose', a, dense(size, size))
    # Strict diagonal dominance keeps these generated systems away from rank thresholds.
    dominant = copy.deepcopy(a)
    for i in range(size):
        dominant['values'][i*size+i] = 1 + sum(abs(v) for v in a['values'][i*size:(i+1)*size])
    rhs = dense(size, 2)
    yield case('diagonally dominant multiple RHS', 'solve', dominant, rhs)
    # A^T A + I is SPD, constructed exactly from bounded dyadic numbers.
    spd = expected('multiply', expected('transpose', a), a)
    spd['values'] = [float(v) for v in spd['values']]
    for i in range(size): spd['values'][i*size+i] += 1
    yield case('Gram plus identity', 'solve_cholesky', spd, rhs)
    yield case('symmetric eigen residual and orthogonality', 'eigen_symmetric', spd)
    tall = matrix(size+2, size, dominant['values'] + dense(2, size)['values'])
    yield case('full rank noisy least squares', 'least_squares', tall, dense(size+2, 2))


def validate(case):
    if case.get('generator_version') != VERSION: raise ValueError('unsupported generator version')
    if case.get('op') not in ('multiply','add','transpose_twice','determinant','determinant_transpose',
                               'product_transpose','solve','solve_cholesky','least_squares','eigen_symmetric'):
        raise ValueError('unsupported property')
    for name in ('a', 'b'):
        a = case.get(name)
        if a is None:
            if name == 'a': raise ValueError('missing A')
            continue
        if not isinstance(a, dict): raise ValueError('matrix must be an object')
        if any(type(a.get(d)) is not int or not 0 <= a[d] <= 6 for d in ('rows','cols')):
            raise ValueError('matrix dimensions must be integers in 0..6')
        values = a.get('values')
        if not isinstance(values, list) or len(values) != a['rows']*a['cols']:
            raise ValueError('wrong value count')
        if any(type(v) not in (int,float) or not math.isfinite(v) or abs(v)>1024 for v in values):
            raise ValueError('expected bounded finite values')
    a, b, op = case['a'], case.get('b'), case['op']
    if op in ('multiply','product_transpose','add','solve','solve_cholesky','least_squares') and b is None:
        raise ValueError('missing B')
    if op in ('multiply','product_transpose') and a['cols'] != b['rows']: raise ValueError('product shape')
    if op == 'add' and (a['rows'],a['cols']) != (b['rows'],b['cols']): raise ValueError('sum shape')
    if op in ('determinant','determinant_transpose','solve','solve_cholesky','eigen_symmetric'):
        if a['rows'] != a['cols'] or a['rows'] > 5: raise ValueError('square matrix required, size <= 5')
    if op in ('solve','solve_cholesky','least_squares') and a['rows'] != b['rows']: raise ValueError('RHS shape')
    if op in ('solve_cholesky','eigen_symmetric'):
        n=a['rows']
        if any(a['values'][i*n+j] != a['values'][j*n+i] for i in range(n) for j in range(n)):
            raise ValueError('symmetric family required')
    if op == 'solve_cholesky':
        # Exact Sylvester criterion keeps reducer candidates in the SPD domain.
        for n in range(1,a['rows']+1):
            if expected('determinant', crop(a,n,n))['value'] <= 0: raise ValueError('SPD required')


def oracle(case):
    validate(case)
    op, a, b = case['op'], case['a'], case.get('b')
    if op == 'transpose_twice': return a
    if op in ('determinant','determinant_transpose'): return expected('determinant',a)
    if op == 'product_transpose': return expected('transpose', expected('multiply',a,b))
    if op in ('solve','solve_cholesky'): return exact_solve(a,b)
    if op == 'least_squares': return exact_least_squares(a,b)
    if op == 'eigen_symmetric': return None
    return expected(op,a,b)


class Failure(Exception):
    def __init__(self, kind, detail):
        self.kind, self.detail = kind, SANITIZER.text(str(detail))[-2000:]
        super().__init__(self.detail)


def request(implementation, op, a, b=None):
    args=['check',op,str(a['rows']),str(a['cols'])]; values=a['values'][:]
    if b is not None: args += [str(b['rows']),str(b['cols'])]; values += b['values']
    try:
        result=harness.execute(implementation['runner']+args,stdin=' '.join(format(float(v),'.17g') for v in values),
                               env=implementation['env'],timeout=30)
    except subprocess.TimeoutExpired: raise Failure('timeout','runner exceeded 30 seconds')
    except OSError as error: raise Failure('launch',error)
    if result.returncode: raise Failure('exit',result.stderr or 'runner exited without diagnostics')
    try: return json.loads(result.stdout)
    except (ValueError,TypeError): raise Failure('protocol','runner did not emit JSON')


def evaluate(implementation, case):
    wanted=oracle(case)  # Invalid reduction candidates must not turn into runner failures.
    a,b,op=case['a'],case.get('b'),case['op']
    def call(operation,left,right=None): return request(implementation,operation,left,right)
    try:
        if op == 'transpose_twice':
            transposed=call('transpose',a); assert_result(transposed,expected('transpose',a))
            actual=call('transpose',transposed)
        elif op == 'product_transpose':
            product=call('multiply',a,b); assert_result(product,expected('multiply',a,b))
            actual=call('transpose',product)
            other=call('multiply',expected('transpose',b),expected('transpose',a))
            assert_result(other,wanted); assert_result(actual,other)
        elif op == 'determinant_transpose':
            actual=call('determinant',a)
            other=call('determinant',expected('transpose',a)); assert_result(other,wanted)
            assert_result(actual,other)
        else: actual=call(op,a,b)
        if op == 'eigen_symmetric': assert_eigen_result(actual,a)
        else: assert_result(actual,wanted)
        return actual
    except (AssertionError,KeyError,TypeError,ValueError,OverflowError) as error:
        raise Failure('numerical',error)


def crop(a, rows, cols):
    return matrix(rows,cols,[a['values'][i*a['cols']+j] for i in range(rows) for j in range(cols)])


def candidates(case):
    """Monotone dimensions, then values. Domain-invalid candidates are skipped."""
    a,b,op=case['a'],case.get('b'),case['op']
    if op in ('multiply','product_transpose'):
        for m,k,n in [(a['rows']-1,a['cols'],b['cols']), (a['rows'],a['cols']-1,b['cols']), (a['rows'],a['cols'],b['cols']-1)]:
            if min(m,k,n)<0: continue
            c=copy.deepcopy(case); c['a']=crop(a,m,k);c['b']=crop(b,k,n);yield c
    elif op in ('determinant','determinant_transpose','solve','solve_cholesky','eigen_symmetric'):
        n=a['rows']-1
        if n>=0:
            c=copy.deepcopy(case);c['a']=crop(a,n,n)
            if b is not None:c['b']=crop(b,n,b['cols'])
            yield c
    elif op in ('transpose_twice','add','least_squares'):
        for r,s in [(a['rows']-1,a['cols']),(a['rows'],a['cols']-1)]:
            if min(r,s)<0:continue
            c=copy.deepcopy(case);c['a']=crop(a,r,s)
            if b is not None:c['b']=crop(b,r,s if op=='add' else b['cols'])
            yield c
    for key in ('a','b'):
        if case.get(key) is None:continue
        for index,value in enumerate(case[key]['values']):
            choices=[0]
            if abs(value)>1: choices.append(math.copysign(1,value))
            for simpler in choices:
                if abs(simpler)>=abs(value):continue
                c=copy.deepcopy(case);c[key]['values'][index]=simpler;yield c


def reduce_case(case, predicate, budget=100):
    current=copy.deepcopy(case); attempts=0; deadline=time.monotonic()+30
    while attempts<budget and time.monotonic()<deadline:
        improved=False
        for candidate in candidates(current):
            if attempts>=budget or time.monotonic()>=deadline:break
            attempts+=1
            try:oracle(candidate)
            except (ValueError,AssertionError,ZeroDivisionError):continue
            if predicate(candidate):current=candidate;improved=True;break
        if not improved:break
    return current,attempts


def run_cases(implementations,cases,output,budget=100):
    output.mkdir(parents=True,exist_ok=True); outcomes=[];failures=[]
    for case in cases:
        oracle(case);results=[]
        for implementation in implementations:
            entry=dict(implementation=implementation['id'],seed=case.get('seed'),family=case['family'],op=case['op'],passed=True)
            try:
                result=evaluate(implementation,case)
                # Compare scalar/matrix outputs only; eigenvector bases need not match.
                if results and case['op']!='eigen_symmetric':
                    try:assert_result(result,results[0][1])
                    except AssertionError as error:raise Failure('differential',error)
                results.append((implementation,result))
            except Failure as failure:
                entry.update(passed=False,kind=failure.kind,detail=failure.detail)
                def persists(candidate):
                    try:
                        actual=evaluate(implementation,candidate)
                        if failure.kind=='differential':
                            other=evaluate(results[0][0],candidate)
                            try:assert_result(actual,other)
                            except AssertionError:return True
                    except Failure as next_failure:return next_failure.kind==failure.kind
                    return False
                reduced,attempts=reduce_case(case,persists,budget)
                artifact=dict(schema_version=1,implementation=implementation['id'],kind=failure.kind,
                              detail=failure.detail,original=case,reduced=reduced,reduction_attempts=attempts,
                              reduction_budget=budget,globally_minimal=False,
                              reproduced_after_reduction=persists(reduced))
                digest=hashlib.sha256(json.dumps(plain(artifact),sort_keys=True).encode()).hexdigest()[:16]
                path=output/('failure-'+digest+'.json');path.write_text(json.dumps(plain(artifact),indent=2)+'\n')
                entry['artifact']=path.name;failures.append(artifact)
            outcomes.append(entry)
        print('Checked seed %s: %s'%(case.get('seed'),case['family']),flush=True)
    return outcomes,failures


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--implementations',nargs='+',choices=harness.NAMES,default=list(harness.NAMES))
    parser.add_argument('--no-build',action='store_true')
    parser.add_argument('--seed',type=int,default=20261002)
    parser.add_argument('--seeds',type=int,default=2)
    parser.add_argument('--reduce-budget',type=int,default=100)
    parser.add_argument('--replay',type=Path,help='replay a saved original/reduced failure or standalone case')
    parser.add_argument('--output',type=Path,default=ROOT/'.build/bugfinding')
    args=parser.parse_args()
    if not 1<=args.seeds<=100 or not 0<=args.reduce_budget<=1000:parser.error('seeds 1..100; reduction budget 0..1000')
    cases=[]
    if args.replay:
        data=json.loads(args.replay.read_text());cases=[data.get('reduced',data)]
    else:
        for seed in range(args.seed,args.seed+args.seeds):cases.extend(generated(seed))
        for path in sorted((ROOT/'tests/bugfinding/regressions').glob('*.json')):
            data=json.loads(path.read_text());cases.append(data.get('reduced',data))
    implementations=[harness.build_one(name,args.no_build) for name in args.implementations]
    checks,failures=run_cases(implementations,cases,args.output,args.reduce_budget)
    summary=dict(schema_version=1,generator_version=VERSION,created_at=datetime.now(timezone.utc).isoformat(),
                 revision=harness.require(['git','rev-parse','HEAD']).strip(),
                 dirty=bool(harness.require(['git','status','--porcelain']).strip()),
                 implementations=[dict(id=i['id'],toolchain=i['toolchain']) for i in implementations],
                 cases=len(cases),checks=checks,failures=len(failures))
    (args.output/'summary.json').write_text(json.dumps(plain(summary),indent=2)+'\n')
    print('%d/%d checks passed'%(sum(c['passed'] for c in checks),len(checks)))
    return 1 if failures else 0


if __name__=='__main__':raise SystemExit(main())
