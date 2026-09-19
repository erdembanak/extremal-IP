"""Check published graphs without importing the enumeration model or nauty."""
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parent


def graph(n, edges):
    assert all(0 <= u < v < n for u, v in edges)
    assert len(edges) == len({tuple(e) for e in edges})
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    return g


def main():
    data = ROOT / 'data'
    for line in (data / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ')
        assert hashlib.sha256((data / name).read_bytes()).hexdigest() == digest
    records = json.loads((data / 'index.json').read_text())['graphs']
    assert [r['name'] for r in records] == ['F58', 'F57_44', 'F57_43', 'F57_33']
    graphs, cones = [], []
    for row in records:
        obj = json.loads((data / row['graph']).read_text())
        assert obj['num_vertices'] == 17
        g = graph(17, obj['edges'])
        assert g.number_of_edges() == row['edges']
        assert max(dict(g.degree()).values()) <= 7
        assert sum(nx.triangles(g).values()) == 0
        assert max(map(len, nx.find_cliques(nx.complement(g)))) == 7
        assert len(nx.max_weight_matching(g, maxcardinality=True)) == 8
        for vertex in g:
            h = g.copy()
            h.remove_node(vertex)
            assert len(nx.max_weight_matching(h, maxcardinality=True)) == 8
        decoded = nx.from_graph6_bytes(obj['graph6'].encode())
        assert set(decoded.edges()) == set(g.edges())
        cobj = json.loads((data / row['cone']).read_text())
        assert cobj['num_vertices'] == 18 and cobj['apex'] == 17
        assert cobj['partition_is_certificate'] is False
        cone = graph(18, cobj['edges_part1'] + cobj['edges_part2'])
        expected = nx.complement(g)
        expected.add_edges_from((v, 17) for v in g)
        assert set(cone.edges()) == set(expected.edges())
        decoded = nx.from_graph6_bytes(cobj['graph6'].encode())
        assert set(decoded.edges()) == set(cone.edges())
        graphs.append(g)
        cones.append(cone)
    for a, b in combinations(graphs, 2):
        assert not nx.is_isomorphic(a, b)
    bags = [range(4), range(4, 8), range(8, 11), range(11, 14), range(14, 17)]
    blowup = {tuple(sorted((u, v))) for i in range(5)
              for u in bags[i] for v in bags[(i + 1) % 5]}
    assert set(graphs[0].edges()) == blowup
    counts = Counter()
    for edge in graphs[0].edges():
        h = graphs[0].copy()
        h.remove_edge(*edge)
        matches = [i for i in range(1, 4) if nx.is_isomorphic(h, graphs[i])]
        assert len(matches) == 1
        counts[matches[0]] += 1
    assert counts == {1: 16, 2: 24, 3: 18}
    for cone in cones[1:]:
        assert set(cones[0].edges()) < set(cone.edges())
    for size in (57, 58):
        saved = json.loads((ROOT / f'runs/full/catalogue-{size}.json').read_text())
        expected = [g for g in graphs if g.number_of_edges() == size]
        assert len(saved) == len(expected)
        for row in saved:
            assert sum(nx.is_isomorphic(graph(17, row['edges']), g) for g in expected) == 1
    status = json.loads((ROOT / 'runs/full/status.json').read_text())
    audit = json.loads((ROOT / 'runs/full/audit.json').read_text())
    assert status['complete'] and status['completed_cases'] == 10
    assert status['unique_graphs'] == {'57': 3, '58': 1}
    assert audit['complete_agreement']
    print('PASS: four graphs, all properties, cone inputs, graph6, checksums,')
    print('      deletion multiplicities 16/24/18, and recorded catalogue agreement.')
    print('Exhaustiveness still relies on the recorded CP-SAT infeasibility results.')


if __name__ == '__main__':
    main()
