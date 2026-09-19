"""Complement-variable audit of the unfiltered 17-vertex catalogue.

No independence-number or factor-criticality constraint is imposed.

Uses the same CP-SAT solver, but no matching witnesses, no isomorphic-copy
exclusions, and the opposite adjacency lex ordering. Also enumerates any
non-factor-critical graphs and reports them separately.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time

import networkx as nx
from ortools.sat.python import cp_model

from enumerate_graphs import properties, write_json
from catalogue_isomorphism import isomorphic


ROOT = Path(__file__).parent / 'runs' / 'unrestricted'


def sequences(n, total, cap):
    if n == 0:
        if total == 0:
            yield ()
        return
    for first in range(min(total, cap), -1, -1):
        rest = total - first
        if rest <= (n - 1) * first:
            for tail in sequences(n - 1, rest, first):
                yield (first,) + tail


def nxgraph(edges):
    graph = nx.Graph()
    graph.add_nodes_from(range(17))
    graph.add_edges_from(edges)
    return graph


def audit_case(e, pattern, degrees, timeout=300):
    started = time.monotonic()
    model = cp_model.CpModel()
    pairs = list(combinations(range(17), 2))
    # y=1 means ABSENCE of an edge in the original triangle-free graph.
    y = {pair: model.new_bool_var(f'absent_{pair[0]}_{pair[1]}') for pair in pairs}
    model.add(sum(y.values()) == 136 - e)
    for i, j, k in combinations(range(17), 3):
        model.add_bool_or([y[i, j], y[i, k], y[j, k]])
    for v in range(17):
        model.add(sum(var for pair, var in y.items() if v in pair) == 16 - degrees[v])
    # Lex-min NONEDGE vector, equivalently lex-max EDGE vector. This is
    # deliberately opposite to the production symmetry orientation.
    for a, b in combinations(range(17), 2):
        if degrees[a] != degrees[b]:
            continue
        others = [v for v in range(17) if v not in (a, b)]
        terms = []
        for index, v in enumerate(others):
            weight = 2 ** (14 - index)
            terms.append(weight * (y[tuple(sorted((a, v)))] - y[tuple(sorted((b, v)))]))
        model.add(sum(terms) <= 0)
    catalogue = json.loads((ROOT / f'catalogue-{e}.json').read_text())
    expected = [nxgraph(row['edges']) for row in catalogue]
    found_classes, noncritical, assignments = set(), [], 0
    complete, final_status = False, 'NOT_RUN'
    while time.monotonic() - started < timeout:
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 2
        solver.parameters.random_seed = 19
        solver.parameters.max_time_in_seconds = timeout - (time.monotonic() - started)
        status = solver.solve(model)
        final_status = solver.status_name(status)
        if status == cp_model.INFEASIBLE:
            complete = True
            break
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            break
        absent = {pair: solver.value(var) for pair, var in y.items()}
        edges = [pair for pair in pairs if not absent[pair]]
        props = properties(17, edges)
        if (not props['triangle_free'] or props['alpha'] > 8 or props['edges'] != e
                or tuple(props['degrees']) != degrees):
            raise RuntimeError('Audit encoding produced an invalid graph')
        assignments += 1
        graph = nxgraph(edges)
        matches = [i for i, target in enumerate(expected) if isomorphic(graph, target)]
        if len(matches) != 1:
            write_json(ROOT / f'audit-unexpected-{e}-{pattern}.json', {'edges': edges})
            raise RuntimeError('Audit found a missing or multiply represented graph class')
        found_classes.add(matches[0])
        if not props['factor_critical']:
            noncritical.append(edges)
        model.add_bool_or([y[pair].Not() if absent[pair] else y[pair] for pair in pairs])
    result = {'edges': e, 'degrees': degrees, 'pattern': pattern,
              'complete': complete, 'last_solver_status': final_status,
              'assignments': assignments, 'catalogue_indices': sorted(found_classes),
              'non_factor_critical_graphs': noncritical,
              'elapsed_seconds': time.monotonic() - started}
    write_json(ROOT / f'audit-e{e}-p{pattern}.json', result)
    print(f'audit e{e}-p{pattern}: {final_status}; {assignments} assignments; '
          f'classes {sorted(found_classes)}; {len(noncritical)} noncritical', flush=True)
    return result


def main():
    cases = [(e, i, deg) for e in (58, 57)
             for i, deg in enumerate(sorted(sequences(17, 2 * e, 7)))]
    if len(cases) != 10:
        raise RuntimeError('Unexpected degree-case count')
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs = [pool.submit(audit_case, *case) for case in cases]
        for job in as_completed(jobs):
            results.append(job.result())
    agreement = all(result['complete'] for result in results)
    for e in (58, 57):
        found = {i for result in results if result['edges'] == e for i in result['catalogue_indices']}
        expected = set(range(len(json.loads((ROOT / f'catalogue-{e}.json').read_text()))))
        agreement &= found == expected
    summary = {'complete_agreement': agreement,
               'method': 'Complement encoding, opposite lex order, no matching constraints, '
                         'no isomorphic-copy cuts; same CP-SAT solver; NetworkX isomorphism comparison',
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'isomorphism_source_sha256': hashlib.sha256(
                   Path(__file__).with_name('catalogue_isomorphism.py').read_bytes()).hexdigest(),
               'cases': sorted(results, key=lambda row: (-row['edges'], row['pattern']))}
    write_json(ROOT / 'audit.json', summary)
    print(json.dumps({'complete_agreement': agreement}), flush=True)
    if not agreement:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
