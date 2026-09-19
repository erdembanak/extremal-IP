"""CP-SAT port of extremal-IP's fixed-order triangle-free graph model.

Enumeration blocks edge assignments, including all matching-witness assignments
for that graph. Completeness is reported only after residual INFEASIBLE.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from collections import deque
from functools import lru_cache
import hashlib
from importlib.metadata import version
from itertools import combinations
import json
import os
from pathlib import Path
import time

import networkx as nx
import pynauty
from ortools.sat.python import cp_model


@dataclass(frozen=True)
class Spec:
    vertices: int = 17
    degree: int = 7
    edges: int | None = 58
    alpha_max: int | None = None
    factor_critical: bool = False
    degree_order: bool = True
    degrees: tuple[int, ...] | None = None
    row_lex: bool = False

    def validate(self):
        n, d = self.vertices, self.degree
        if n < 1 or not 0 <= d < n:
            raise ValueError('Require vertices >= 1 and 0 <= degree < vertices')
        if self.edges is not None and not 0 <= self.edges <= n * d // 2:
            raise ValueError('Edge count must be between 0 and floor(n*d/2)')
        if self.alpha_max is not None and not 1 <= self.alpha_max <= n:
            raise ValueError('Require 1 <= alpha-max <= vertices')
        if self.factor_critical and n % 2 == 0:
            raise ValueError('Factor-critical graphs must have odd order')
        if self.row_lex and n > 32:
            raise ValueError('Row-lex encoding supports at most 32 vertices')
        if self.degrees is not None:
            if len(self.degrees) != n or any(not 0 <= v <= d for v in self.degrees):
                raise ValueError('Invalid degree sequence')
            if self.edges is None or sum(self.degrees) != 2 * self.edges:
                raise ValueError('Degree sequence must sum to twice the fixed edge count')
            if tuple(sorted(self.degrees, reverse=True)) != self.degrees:
                raise ValueError('Degree sequence must be nonincreasing')


def degree_sequences(n, d, edges):
    """All nonincreasing degree multisets; some may not be graphical."""
    def partitions(total, cap, slots):
        if total == 0:
            yield ()
        elif slots:
            for first in range(min(total, cap), 0, -1):
                for rest in partitions(total - first, first, slots - 1):
                    yield (first,) + rest

    deficit = n * d - 2 * edges
    if deficit < 0:
        return []
    return sorted({tuple(sorted([d - x for x in part] + [d] * (n - len(part)),
                               reverse=True))
                   for part in partitions(deficit, d, n)})


def build_model(spec):
    spec.validate()
    n = spec.vertices
    model = cp_model.CpModel()
    x = {edge: model.new_bool_var(f'e_{edge[0]}_{edge[1]}')
         for edge in combinations(range(n), 2)}
    for i, j, k in combinations(range(n), 3):
        model.add(x[i, j] + x[i, k] + x[j, k] <= 2)
    degrees = [sum(var for edge, var in x.items() if v in edge) for v in range(n)]
    for v, degree in enumerate(degrees):
        model.add(degree <= spec.degree)
        if spec.degrees is not None:
            model.add(degree == spec.degrees[v])
    if spec.degree_order:
        for v in range(n - 1):
            model.add(degrees[v] >= degrees[v + 1])
    if spec.row_lex:
        # A lex-min upper-triangular adjacency vector under degree-preserving
        # permutations satisfies every pair-transposition lex comparison.
        # For i<j, the first changed edge is determined by the first vertex k
        # (in increasing order, excluding i,j) where x_ik != x_jk.
        for i, j in combinations(range(n), 2):
            if spec.degrees is not None and spec.degrees[i] != spec.degrees[j]:
                continue
            others = [k for k in range(n) if k not in (i, j)]
            lhs = sum((1 << (len(others) - t - 1)) * x[min(i, k), max(i, k)]
                      for t, k in enumerate(others))
            rhs = sum((1 << (len(others) - t - 1)) * x[min(j, k), max(j, k)]
                      for t, k in enumerate(others))
            constraint = model.add(lhs <= rhs)
            if spec.degrees is None:
                equal_degree = model.new_bool_var(f'equal_degree_{i}_{j}')
                model.add(degrees[i] == degrees[j]).only_enforce_if(equal_degree)
                model.add(degrees[i] != degrees[j]).only_enforce_if(equal_degree.Not())
                constraint.only_enforce_if(equal_degree)
    if spec.edges is None:
        model.maximize(sum(x.values()))
    else:
        model.add(sum(x.values()) == spec.edges)
    if spec.alpha_max is not None:
        for subset in combinations(range(n), spec.alpha_max + 1):
            model.add_bool_or([x[edge] for edge in combinations(subset, 2)])
    if spec.factor_critical:
        # Existential perfect matching after each vertex deletion. Blocking only
        # x prevents repeatedly enumerating alternative matching witnesses.
        for removed in range(n):
            matching = {edge: model.new_bool_var(f'm_{removed}_{edge[0]}_{edge[1]}')
                        for edge in x if removed not in edge}
            for edge, var in matching.items():
                model.add(var <= x[edge])
            for v in range(n):
                if v != removed:
                    model.add_exactly_one([var for edge, var in matching.items() if v in edge])
    return model, x


def block_graph(model, x, edges):
    # Safe because each enumeration run fixes the cardinality of E.
    # With zero edges this adds false, excluding the unique empty graph.
    model.add(sum(x[edge] for edge in edges) <= len(edges) - 1)


def properties(n, edges):
    """Independent graph verification: bitset independence, NetworkX blossom."""
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    graph.add_edges_from(edges)
    if len(edges) != graph.number_of_edges() or any(u == v for u, v in edges):
        raise ValueError('Duplicate edge or loop')
    adj = [sum(1 << u for u in graph[v]) for v in range(n)]

    @lru_cache(None)
    def alpha(mask):
        if not mask:
            return 0
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        return max(alpha(rest), 1 + alpha(rest & ~adj[v]))

    def matching_size(g):
        return len(nx.max_weight_matching(g, maxcardinality=True))

    return {
        'triangle_free': all(not (adj[u] & adj[v]) for u, v in edges),
        'edges': graph.number_of_edges(),
        'degrees': [graph.degree(v) for v in range(n)],
        'alpha': alpha((1 << n) - 1),
        'matching': matching_size(graph),
        'factor_critical': n % 2 == 1 and all(
            matching_size(graph.subgraph([u for u in range(n) if u != v])) == (n - 1) // 2
            for v in range(n)),
    }


def check_graph(spec, edges):
    if any(not 0 <= u < v < spec.vertices for u, v in edges):
        raise ValueError('Edges must be sorted pairs within the vertex set')
    props = properties(spec.vertices, edges)
    valid = props['triangle_free'] and max(props['degrees']) <= spec.degree
    valid &= spec.edges is None or props['edges'] == spec.edges
    valid &= spec.alpha_max is None or props['alpha'] <= spec.alpha_max
    valid &= not spec.factor_critical or props['factor_critical']
    valid &= not spec.degree_order or props['degrees'] == sorted(props['degrees'], reverse=True)
    valid &= spec.degrees is None or tuple(props['degrees']) == spec.degrees
    if not valid:
        raise ValueError(f'Graph fails independent checks: {props}')
    return props


def certificate(n, edges):
    adj = {v: [] for v in range(n)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return pynauty.certificate(pynauty.Graph(n, adjacency_dict=adj)).hex()


def relabelled_copies(n, edges, limit):
    """Up to limit distinct isomorphic copies, preserving each vertex's degree.

    Adjacent transpositions within equal-degree classes generate their full
    permutation group. A bounded BFS gives valid extra cuts, not a claim of
    exhausting that orbit. Exact edge sets suppress automorphism duplicates.
    """
    original = tuple(sorted(edges))
    yield original
    if limit == 0:
        return
    degrees = [sum(v in edge for edge in original) for v in range(n)]
    groups = {}
    for v, degree in enumerate(degrees):
        groups.setdefault(degree, []).append(v)
    swaps = [(a, b) for group in groups.values() for a, b in zip(group, group[1:])]
    seen, queue = {original}, deque([original])
    while queue:
        current = queue.popleft()
        for a, b in swaps:
            def perm(v):
                return b if v == a else a if v == b else v
            image = tuple(sorted(tuple(sorted((perm(u), perm(v)))) for u, v in current))
            if image not in seen:
                seen.add(image)
                queue.append(image)
                yield image
                if len(seen) >= limit + 1:
                    return


def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def make_solver(seconds, workers, seed):
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    return solver


def run(spec, output, seconds=60, workers=1, seed=0, max_solutions=0, resume=False,
        relabel_blocks=0):
    """Return summary. Time/solution limits are per invocation, not proof bounds."""
    spec.validate()
    if seconds <= 0 or workers < 1 or max_solutions < 0 or relabel_blocks < 0:
        raise ValueError('Positive time/workers and nonnegative solution limit required')
    if spec.edges is None and resume:
        raise ValueError('Resume is available only for fixed-edge enumeration')
    start = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        'schema': 1, 'spec': asdict(spec),
        'relabel_blocks': relabel_blocks,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'packages': {name: version(name) for name in ('ortools', 'networkx', 'pynauty')},
    }
    # JSON normalizes tuple-valued degree sequences into lists.
    manifest = json.loads(json.dumps(manifest))
    records = output / 'solutions.jsonl'
    if resume:
        if json.loads((output / 'manifest.json').read_text()) != manifest:
            raise ValueError('Resume requires the same model, code, and package versions')
    else:
        if any(output.iterdir()):
            raise ValueError('Output directory is not empty; use --resume or a fresh directory')
        write_json(output / 'manifest.json', manifest)
        records.touch()
    # A failed/interrupted invocation must never leave a stale COMPLETE marker.
    write_json(output / 'summary.json', {'complete': False, 'stop_reason': 'RUNNING'})
    model, x = build_model(spec)
    seen, assignments, blocked = set(), set(), set()

    def exclude(edges):
        for image in relabelled_copies(spec.vertices, edges, relabel_blocks):
            if image not in blocked:
                block_graph(model, x, image)
                blocked.add(image)
    count = 0
    if resume:
        with records.open() as src:
            for line in src:
                row = json.loads(line)  # Fail closed on a truncated/corrupt record.
                edges = tuple(tuple(edge) for edge in row['edges'])
                props = check_graph(spec, edges)
                cert = certificate(spec.vertices, edges)
                if (props != row['properties'] or cert != row['certificate'] or
                        edges in assignments or row['new_isomorphism'] != (cert not in seen)):
                    raise ValueError('Invalid checkpoint record')
                seen.add(cert)
                assignments.add(edges)
                count += 1
                exclude(edges)

    complete, reason, status_name = False, 'TIME_LIMIT', 'NOT_RUN'
    found_now, last_bound = 0, None
    with records.open('a') as sink:
        while True:
            remaining = seconds - (time.monotonic() - start)
            if remaining <= 0:
                break
            if max_solutions and found_now >= max_solutions:
                reason = 'SOLUTION_LIMIT'
                break
            solver = make_solver(remaining, workers, seed)
            status = solver.solve(model)
            status_name = solver.status_name(status)
            if status == cp_model.INFEASIBLE:
                complete, reason = True, 'INFEASIBLE'
                break
            if status == cp_model.MODEL_INVALID:
                raise RuntimeError(solver.response_stats())
            if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                reason = status_name
                break
            edges = tuple(edge for edge, var in x.items() if solver.value(var))
            if edges in blocked:
                raise RuntimeError('A blocked assignment was returned again')
            props = check_graph(spec, edges)
            cert = certificate(spec.vertices, edges)
            is_new = cert not in seen
            row = {'edges': edges, 'properties': props, 'certificate': cert,
                   'new_isomorphism': is_new}
            sink.write(json.dumps(row) + '\n')
            sink.flush()
            os.fsync(sink.fileno())
            assignments.add(edges)
            seen.add(cert)
            count += 1
            found_now += 1
            print(f'graph {count}: edges={len(edges)}, unique={len(seen)}, '
                  f'alpha={props["alpha"]}, factor_critical={props["factor_critical"]}', flush=True)
            if spec.edges is None:
                last_bound = solver.best_objective_bound
                complete, reason = status == cp_model.OPTIMAL, status_name
                break
            exclude(edges)
    summary = {
        'mode': 'optimize' if spec.edges is None else 'enumerate',
        'complete': complete, 'stop_reason': reason, 'last_solver_status': status_name,
        'labelled_solutions': count, 'nonisomorphic_graphs': len(seen),
        'blocked_assignments': len(blocked),
        'new_solutions_this_run': found_now, 'elapsed_seconds': time.monotonic() - start,
        'workers': workers, 'seed': seed,
    }
    if spec.edges is None:
        summary['objective_upper_bound'] = last_bound
        summary['complete_meaning'] = 'Optimization settled; NOT an enumeration'
    else:
        summary['complete_meaning'] = 'Every isomorphism class satisfying manifest spec covered'
    write_json(output / 'summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vertices', type=int, default=17)
    parser.add_argument('--degree', type=int, default=7)
    parser.add_argument('--edges', type=int, default=58)
    parser.add_argument('--alpha-max', type=int)
    parser.add_argument('--factor-critical', action='store_true')
    parser.add_argument('--no-degree-order', action='store_true')
    parser.add_argument('--row-lex', action='store_true',
                        help='Safe pair-transposition symmetry constraints within degree classes')
    parser.add_argument('--degree-pattern', type=int, help='Zero-based index from --list-patterns')
    parser.add_argument('--list-patterns', action='store_true')
    parser.add_argument('--optimize', action='store_true', help='Maximize edges instead of enumerating')
    parser.add_argument('--time-limit', type=float, default=60, help='Seconds for this invocation')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--max-solutions', type=int, default=0, help='0 means no solution-count limit')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--relabel-blocks', type=int, default=0,
                        help='Also exclude up to this many isomorphic copies per graph')
    args = parser.parse_args()
    try:
        spec = Spec(args.vertices, args.degree, None if args.optimize else args.edges,
                    args.alpha_max, args.factor_critical, not args.no_degree_order,
                    row_lex=args.row_lex)
        spec.validate()
        if args.optimize and (args.list_patterns or args.degree_pattern is not None):
            raise ValueError('Degree patterns require a fixed edge count')
        if args.list_patterns or args.degree_pattern is not None:
            patterns = degree_sequences(args.vertices, args.degree, args.edges)
            if args.list_patterns:
                for index, pattern in enumerate(patterns):
                    print(index, ','.join(map(str, pattern)))
                return
            if not 0 <= args.degree_pattern < len(patterns):
                raise ValueError('Degree-pattern index out of range')
            spec = Spec(**{**asdict(spec), 'degrees': patterns[args.degree_pattern]})
        if args.output is None:
            raise ValueError('--output is required')
        summary = run(spec, args.output, args.time_limit, args.workers, args.seed,
                      args.max_solutions, args.resume, args.relabel_blocks)
        print(json.dumps(summary, indent=2))
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(2, f'error: {exc}\n')


if __name__ == '__main__':
    main()
