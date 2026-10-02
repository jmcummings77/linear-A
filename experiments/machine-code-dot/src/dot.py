from array import array
import json
import math
import sys
from time import perf_counter_ns

def dot(a, b, offset, n):
    total = 0.0
    for i in range(offset, offset + n):
        total += a[i] * b[i]
    return total

def main():
    n, iterations, seed = map(int, sys.argv[1:])
    assert 0 <= n <= 1048576 and 1 <= iterations <= 100000000 and 0 <= seed <= 1000000
    count = max(n, 1) * 32
    a = array('d', (((i * 17 + seed * 13) % 101 - 50) / 16 for i in range(count)))
    b = array('d', (((i * 29 + (seed + 1) * 7) % 103 - 51) / 16 for i in range(count)))
    for i in range(max(8, min(iterations, 128))):
        dot(a, b, (i % 32) * n, n)
    checksum = 0.0
    start = perf_counter_ns()
    for i in range(iterations):
        checksum += dot(a, b, (i % 32) * n, n)
    elapsed = perf_counter_ns() - start
    assert math.isfinite(checksum)
    print(json.dumps(dict(elapsed_ns=elapsed, iterations=iterations, checksum=checksum)))

if __name__ == '__main__':
    main()
