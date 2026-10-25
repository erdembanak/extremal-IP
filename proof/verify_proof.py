"""Finite checks supporting paper/triangle-free-extremal-formula.tex; not a proof of the theorem.

Uses only the Python standard library. Checks small graph facts, exact counting
identities, and the finite table of five-cycle capacity patterns.
"""

from itertools import combinations, product


def edge(a, b):
    return tuple(sorted((a, b)))


def independent_sets(vertices, edges, size):
    return [
        subset
        for subset in combinations(vertices, size)
        if all(edge(a, b) not in edges for a, b in combinations(subset, 2))
    ]


def triangle_free(vertices, edges):
    return not any(
        all(edge(a, b) in edges for a, b in combinations(triple, 2))
        for triple in combinations(vertices, 3)
    )


def check_small_graphs():
    # H: K_{2,2} between {a,b} and {u,v}, and paths a-p-q-b, u-x-y-v.
    a, b, u, v, p, q, x, y = range(8)
    vertices = tuple(range(8))
    H = {
        edge(s, t)
        for s, t in (
            (a, u), (a, v), (b, u), (b, v),
            (a, p), (p, q), (q, b), (u, x), (x, y), (y, v),
        )
    }
    assert len(H) == 10 and triangle_free(vertices, H)
    assert independent_sets(vertices, H, 3)
    assert not independent_sets(vertices, H, 4)

    allowed = {
        edge(s, t)
        for s, t in combinations(vertices, 2)
        if edge(s, t) not in H
        and triangle_free(vertices, H | {edge(s, t)})
    }
    assert allowed == {edge(s, t) for s in (p, q) for t in (x, y)}
    for size in range(5):
        for extra in combinations(allowed, size):
            feasible = triangle_free(vertices, H | set(extra))
            is_matching = len({z for e in extra for z in e}) == 2 * size
            assert feasible == is_matching

    M = {edge(j, (j + s) % 8) for j in range(8) for s in (1, 4)}
    assert len(M) == 12 and triangle_free(vertices, M)
    neighborhoods = {
        tuple(j for j in vertices if edge(i, j) in M) for i in vertices
    }
    assert len(neighborhoods) == 8
    assert set(independent_sets(vertices, M, 3)) == neighborhoods
    assert not independent_sets(vertices, M, 4)

    # One perfect matching augmentation is M in this cyclic order.
    cycle = (a, p, q, b, u, x, y, v)
    mapped_M = {edge(cycle[i], cycle[j]) for i, j in M}
    assert mapped_M == H | {edge(p, x), edge(q, y)}

    # Verify the neighborhood-sum identity formally, using coefficient vectors.
    for i in vertices:
        coefficients = [0] * 8
        for j in ((i - 1) % 8, i, (i + 1) % 8):
            for k in vertices:
                coefficients[k] += int(edge(j, k) in M)
        assert coefficients == [1 + int(k == i) for k in vertices]

    # The contact configuration used in the five-bag extension contains H.
    # Labels 0..4 are the cycle; 5 is exceptional; 6,7 belong to bags 0,1.
    contact = {edge(j, (j + 1) % 5) for j in range(5)} | {
        edge(5, 6), edge(5, 7), edge(6, 4), edge(6, 1),
        edge(7, 0), edge(7, 2),
    }
    mapping = dict(zip((a, b, u, v, p, q, x, y), (0, 2, 7, 1, 4, 3, 5, 6)))
    assert {edge(mapping[i], mapping[j]) for i, j in H} <= contact
    return "small graph facts: passed"


def check_capacity_patterns():
    # All ordered nonnegative five-tuples summing to three, not just orbit reps.
    patterns = [p for p in product(range(4), repeat=5) if sum(p) == 3]
    assert len(patterns) == 35
    independent_pairs = [
        pair for pair in combinations(range(5), 2)
        if (pair[0] - pair[1]) % 5 not in (1, 4)
    ]
    independent_index_sets = [(i,) for i in range(5)] + independent_pairs
    for t in range(3, 101):
        for p in patterns:
            sizes = [2 * t + 1 - p[(i + 2) % 5] - p[(i + 3) % 5]
                     for i in range(5)]
            assert sum(sizes) == 10 * t - 1
            assert all(sizes[(i - 1) % 5] + sizes[(i + 1) % 5]
                       == 4 * t - 1 + p[i] for i in range(5))
            weights = [sizes[i] * p[i] for i in range(5)]
            lower_bound = max(
                2 * sum(weights[i] for i in I) - sum(weights)
                for I in independent_index_sets
            )
            assert lower_bound >= 2 * t - 1
    return "35 capacity patterns, 3 <= t <= 100: passed"


def check_one_exception_algebra():
    checked = 0
    for t in range(2, 21):
        h = 2 * t
        for x0 in range(1, h + 1):
            for u in range(h + 2):
                for v in range(h + 2):
                    if x0 + u + v > 2 * h + 1:
                        continue
                    a = h - x0
                    if u and v:
                        bound = (h + 1 + 2 * a - u - v
                                 + u * (v - 1) + v * (u - 1))
                        assert bound == h - 1 + 2 * a + 2 * (u - 1) * (v - 1)
                    else:
                        bound = h + 1 + 2 * a
                    assert bound >= h - 1
                    checked += 1
    return f"one-exception algebra, {checked} integer tuples: passed"


def check_odd_excess_imbalance():
    # Check the new general inequality beyond q=1,3. Enumerate every load
    # vector of the given total and every bipartition of its positive support.
    def compositions(total, length):
        if length == 1:
            yield (total,)
        else:
            for first in range(total + 1):
                for rest in compositions(total - first, length - 1):
                    yield (first,) + rest

    checked = 0
    for q in (1, 3, 5, 7, 9):
        k = (q + 1) // 2
        for p in compositions(q, 5):
            support = [i for i in range(5) if p[i]]
            partitions = []
            for colors in product((0, 1), repeat=len(support)):
                color = dict(zip(support, colors))
                if any(color[i] == color[(i + 1) % 5]
                       for i in support if (i + 1) % 5 in color):
                    continue
                I = [i for i in support if color[i] == 0]
                if 2 * sum(p[i] for i in I) > q:
                    partitions.append(I)
            # A proper induced subgraph of C5 is bipartite; full support is not.
            assert bool(partitions) == (len(support) < 5)
            for h in (q, q + 1, q + 7):
                w = [p[i] * (h - p[(i + 2) % 5] - p[(i + 3) % 5])
                     for i in range(5)]
                for I in partitions:
                    P = sum(p[i] for i in I)
                    signed_weight = 2 * sum(w[i] for i in I) - sum(w)
                    scalar_bound = h * (2 * P - q) - P * (P - 1)
                    assert signed_weight >= scalar_bound >= h - k * (k - 1)
                    checked += 1
    return f"odd-excess bipartition inequality, {checked} instances: passed"


def check_knapsack():
    checked = 0
    for d in range(4, 101):
        bonus = d // 2
        # Include all z satisfying the hypotheses used in the exchange proof.
        for z in range(d + 1, d + bonus):
            dp = [0] * (5 * z + 1)
            for m in range(1, len(dp)):
                dp[m] = dp[m - 1]
                for j in range(d, min(m, z) + 1):
                    utility = bonus if j == z else j - d + 1
                    dp[m] = max(dp[m], dp[m - j] + utility)
                k, r = divmod(m, z)
                formula = k * bonus + (r - d + 1 if r >= d else 0)
                assert dp[m] == formula, (d, z, m, dp[m], formula)
                checked += 1
    return f"knapsack dynamic programming, {checked} parameter triples: passed"


def check_construction_family():
    # Independently verify the extended B(d,d+t) construction, including
    # all one-vertex deletions, over the range of the written Hall argument.
    graphs = deletions = 0
    for d in range(3, 25):
        for t in range((d - 2) // 2 + 1):
            i = d + t
            L, R = set(range(i)), set(range(i, 2 * i))
            HL, HR = set(range(d - 1, i)), set(range(i + d - 1, 2 * i))
            v = 2 * i
            E = {edge(a, b) for a in L for b in R
                 if not (a in HL and b in HR)
                 and not (a < d - 1 and b < i + d - 1
                          and (b - i - a) % (d - 1) < t)}
            E |= {edge(v, x) for x in HL | HR}
            N = {x: set() for x in range(v + 1)}
            for a, b in E:
                N[a].add(b)
                N[b].add(a)
            assert [len(N[x]) for x in range(v)] == [d] * v
            assert len(N[v]) == 2 * t + 2 <= d
            assert len(E) == d * (d + t) + t + 1
            assert all(not (N[a] & N[b]) for a, b in E)

            def perfect_matching(left, right):
                mate = {}

                def augment(x, seen):
                    for y in sorted(N[x] & right):
                        if y not in seen:
                            seen.add(y)
                            if y not in mate or augment(mate[y], seen):
                                mate[y] = x
                                return True
                    return False

                return all(augment(x, set()) for x in sorted(left))

            assert perfect_matching(L, R)
            for deleted in range(v):
                if deleted in L:
                    left, right = L - {deleted}, R - {min(HR)}
                else:
                    left, right = L - {min(HL)}, R - {deleted}
                assert perfect_matching(left, right)
            graphs += 1
            deletions += v + 1
    return f"construction family, {graphs} graphs and {deletions} deletions: passed"


if __name__ == "__main__":
    for check in (check_small_graphs, check_capacity_patterns,
                  check_one_exception_algebra, check_odd_excess_imbalance,
                  check_knapsack, check_construction_family):
        print(check(), flush=True)
    print("These finite checks supplement the written proof; they do not replace it.")
