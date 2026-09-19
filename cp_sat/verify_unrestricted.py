"""Check unfiltered graph files independently of CP-SAT and nauty."""
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
from catalogue_isomorphism import isomorphic

ROOT = Path(__file__).resolve().parent


def main():
    data = ROOT / 'data/unrestricted'
    index = json.loads((data / 'index.json').read_text())
    assert index['complete'] and index['independence_constraint'] is None
    assert index['factor_critical_constraint'] is False
    for line in (data / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ')
        assert hashlib.sha256((data / name).read_bytes()).hexdigest() == digest
    graphs, entries, distribution = {}, {}, Counter()
    for row in index['graphs']:
        obj = json.loads((data / row['file']).read_text())
        assert obj['num_vertices'] == 17
        edges = obj['edges']
        assert all(0 <= u < v < 17 for u, v in edges)
        assert len(edges) == len({tuple(e) for e in edges}) == row['edges']
        g = nx.Graph()
        g.add_nodes_from(range(17))
        g.add_edges_from(edges)
        assert max(dict(g.degree()).values()) <= 7
        assert not any(nx.triangles(g).values())
        alpha = max(map(len, nx.find_cliques(nx.complement(g))))
        assert alpha == row['alpha'] == obj['properties']['alpha']
        assert len(nx.max_weight_matching(g, maxcardinality=True)) == 8
        critical = True
        for v in g:
            h = g.copy()
            h.remove_node(v)
            critical &= len(nx.max_weight_matching(h, maxcardinality=True)) == 8
        assert critical == row['factor_critical'] == obj['properties']['factor_critical']
        decoded = nx.from_graph6_bytes(obj['graph6'].encode())
        assert set(decoded.edges()) == set(g.edges())
        graphs[row['name']], entries[row['name']] = g, row
        distribution[(row['edges'], alpha, critical)] += 1
    assert distribution == {(58, 7, True): 1, (58, 8, True): 1,
                            (57, 7, True): 3, (57, 8, True): 18}
    for a, b in combinations(graphs.values(), 2):
        assert not isomorphic(a, b)
    parents = {name: g for name, g in graphs.items() if g.number_of_edges() == 58}
    children = {name: g for name, g in graphs.items() if g.number_of_edges() == 57}
    actual = {name: Counter() for name in children}
    for name, g in parents.items():
        for edge in g.edges():
            h = g.copy()
            h.remove_edge(*edge)
            matches = [child for child, target in children.items() if isomorphic(h, target)]
            assert len(matches) == 1
            actual[matches[0]][name] += 1
    for name, counts in actual.items():
        assert dict(counts) == entries[name]['single_edge_deletion_parents']
    assert sum(bool(counts) for counts in actual.values()) == 6
    for size in (57, 58):
        saved = list(nx.read_graph6(data / f'graphs-{size}.g6'))
        expected = [g for g in graphs.values() if g.number_of_edges() == size]
        assert len(saved) == len(expected)
        assert all(set(a.edges()) == set(b.edges()) for a, b in zip(saved, expected))
    print('PASS: 2 classes at 58 edges, 21 at 57; alpha and matching properties,')
    print('      distinct isomorphism classes, graph6, checksums, all edge-deletion links.')
    print('Exhaustiveness relies on the separate recorded solver runs.')


if __name__ == '__main__':
    main()
