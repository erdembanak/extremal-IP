"""Export the four explicit representatives and compare them to the catalogue."""
from collections import Counter
import json
from pathlib import Path

import networkx as nx

from enumerate_graphs import certificate, properties, write_json


def main():
    root = Path(__file__).parent / 'runs' / 'full'
    bags = [range(0, 4), range(4, 8), range(8, 11), range(11, 14), range(14, 17)]
    base = {tuple(sorted((u, v))) for i in range(5)
            for u in bags[i] for v in bags[(i + 1) % 5]}
    records = []
    for name, removed in [('F58', None), ('F57_44', (0, 4)),
                          ('F57_43', (4, 8)), ('F57_33', (8, 11))]:
        edges = base if removed is None else base - {removed}
        props = properties(17, sorted(edges))
        cert = certificate(17, edges)
        saved = json.loads((root / f'catalogue-{len(edges)}.json').read_text())
        if sum(row['certificate'] == cert for row in saved) != 1:
            raise RuntimeError('Explicit representative disagrees with the catalogue')
        graph = nx.Graph()
        graph.add_nodes_from(range(17))
        graph.add_edges_from(edges)
        cone = nx.complement(graph)
        cone.add_node(17)
        cone.add_edges_from((v, 17) for v in range(17))
        records.append({'name': name, 'bags': [list(b) for b in bags],
                        'deleted_edge': removed, 'edges': sorted(edges),
                        'properties': props, 'certificate': cert,
                        'graph6': nx.to_graph6_bytes(graph, header=False).decode().strip(),
                        'cone_edges': sorted(tuple(sorted(e)) for e in cone.edges()),
                        'cone_graph6': nx.to_graph6_bytes(cone, header=False).decode().strip()})
        print(name, dict(Counter(props['degrees'])), props['alpha'], props['factor_critical'])
    for e in (57, 58):
        saved = {r['certificate'] for r in json.loads((root / f'catalogue-{e}.json').read_text())}
        explicit = {r['certificate'] for r in records if len(r['edges']) == e}
        if saved != explicit:
            raise RuntimeError('Explicit descriptions do not cover the entire catalogue')
    deletion_classes = {certificate(17, base - {edge}) for edge in base}
    if deletion_classes != {r['certificate'] for r in records[1:]}:
        raise RuntimeError('Unexpected single-edge-deletion class')
    write_json(root / 'explicit-representatives.json', records)
    (root / 'candidates.g6').write_text('\n'.join(r['graph6'] for r in records) + '\n')
    (root / 'cones.g6').write_text('\n'.join(r['cone_graph6'] for r in records) + '\n')


if __name__ == '__main__':
    main()
