"""Classical graph fixtures and an analytic prism control, not a new expander construction.

The Petersen graph and triangular prism are cubic nonbipartite Ramanujan
graphs. Larger prisms C_m square K_2 have known spectra but are not a
Ramanujan family. Shifted graph Laplacians are SPD by their edge energy.
"""
import math


SHIFT = 0.25


def petersen_graph():
    """Outer pentagon, inner pentagram, and matching spokes; vertices 0..9."""
    edges = set()
    for i in range(5):
        edges.add(tuple(sorted((i, (i + 1) % 5))))
        edges.add(tuple(sorted((5 + i, 5 + (i + 2) % 5))))
        edges.add((i, 5 + i))
    return 10, sorted(edges)


def prism_graph(cycle_size):
    """C_m square K_2, with two m-cycles joined by a matching."""
    if type(cycle_size) is not int or cycle_size < 3:
        raise ValueError('a prism requires an integer cycle size of at least three')
    edges = set()
    for i in range(cycle_size):
        for layer in (0, 1):
            offset = layer * cycle_size
            edges.add(tuple(sorted((offset + i, offset + (i + 1) % cycle_size))))
        edges.add((i, cycle_size + i))
    return 2 * cycle_size, sorted(edges)


def prism_spectrum(cycle_size):
    """Adjacency eigenvalues from cycle Fourier modes times the K_2 eigenvectors."""
    prism_graph(cycle_size)  # Validate the cycle size, including excluding booleans.
    return sorted(2 * math.cos(2 * math.pi * k / cycle_size) + sign
                  for k in range(cycle_size) for sign in (-1, 1))


def _validated_edges(vertex_count, edges, shift):
    if type(vertex_count) is not int or vertex_count < 1:
        raise ValueError('graph vertex count must be a positive integer')
    if not math.isfinite(shift) or shift <= 0:
        raise ValueError('Laplacian shift must be finite and positive')
    edges = list(edges)
    if any(len(edge) != 2 or any(type(i) is not int for i in edge)
           or not 0 <= edge[0] < edge[1] < vertex_count for edge in edges):
        raise ValueError('edges must be undirected pairs with 0 <= u < v < n')
    if len(set(edges)) != len(edges):
        raise ValueError('graph edges must be unique')
    return edges


def shifted_laplacian(vertex_count, edges, shift=SHIFT):
    """Canonical CSR for shift I + D - H, without a dense intermediate."""
    edges = _validated_edges(vertex_count, edges, shift)
    neighbors = [[] for _ in range(vertex_count)]
    for u, v in edges:
        neighbors[u].append(v)
        neighbors[v].append(u)
    offsets, indices, values = [0], [], []
    for i, adjacent in enumerate(neighbors):
        entries = {j: -1.0 for j in adjacent}
        entries[i] = shift + len(adjacent)
        for j in sorted(entries):
            indices.append(j)
            values.append(entries[j])
        offsets.append(len(indices))
    return dict(rows=vertex_count, cols=vertex_count, offsets=offsets,
                indices=indices, values=values)


def graph_product(vertex_count, edges, vector, shift=SHIFT):
    """Independent edge-difference action; does not read or multiply a CSR matrix."""
    edges = _validated_edges(vertex_count, edges, shift)
    if len(vector) != vertex_count or not all(math.isfinite(x) for x in vector):
        raise ValueError('graph vector must have matching size and finite entries')
    result = [shift * x for x in vector]
    for u, v in edges:
        difference = vector[u] - vector[v]
        result[u] += difference
        result[v] -= difference
    return result


def graph_energy(vertex_count, edges, vector, shift=SHIFT):
    """x^T A x = shift sum(x_i^2) + sum_edges (x_u - x_v)^2 > 0 for x != 0."""
    edges = _validated_edges(vertex_count, edges, shift)
    if len(vector) != vertex_count or not all(math.isfinite(x) for x in vector):
        raise ValueError('graph vector must have matching size and finite entries')
    return shift * math.fsum(x*x for x in vector) + math.fsum(
        (vector[u] - vector[v])**2 for u, v in edges)


def manufactured_solution(vertex_count):
    """Deterministic dyadic variation, avoiding the constant eigenvector alone."""
    return [1 + ((17*i + 3) % 31 - 15) / 32 for i in range(vertex_count)]


def fixtures():
    cases = []
    for name, (n, edges) in [('Petersen', petersen_graph()),
                            ('triangular prism', prism_graph(3))]:
        a = shifted_laplacian(n, edges)
        truth = manufactured_solution(n)
        rhs = graph_product(n, edges, truth)
        for op, jacobi, suffix in [('spmv', 0, 'SpMV'), ('cg', 0, 'CG'),
                                  ('cg', 1, 'Jacobi CG')]:
            cases.append(dict(name='graph ' + name + ' ' + suffix, a=a,
                              b=truth if op == 'spmv' else rhs, op=op,
                              expected=rhs if op == 'spmv' else truth,
                              reason='converged', invalid=False,
                              options=dict(rtol=1e-10, atol=0, limit=100,
                                           jacobi=jacobi, capture=1)))
    return cases
