"""Export completed unfiltered catalogues as individual JSON and graph6 files."""
from collections import Counter
import hashlib
import json
from pathlib import Path

import networkx as nx

from enumerate_graphs import certificate, properties, write_json

ROOT = Path(__file__).resolve().parent


def main():
    run = ROOT / 'runs/unrestricted'
    status = json.loads((run / 'status.json').read_text())
    if not status['complete'] or status['completed_cases'] != 10:
        raise RuntimeError('Refusing to export an incomplete catalogue as complete')
    output = ROOT / 'data/unrestricted'
    output.mkdir(parents=True, exist_ok=True)
    index = {'schema_version': 1, 'num_vertices': 17, 'maximum_degree': 7,
             'triangle_free': True, 'independence_constraint': None,
             'factor_critical_constraint': False, 'complete': True,
             'graphs': [], 'layers': {}}
    deletions = {}
    for size in (58, 57):
        rows = sorted(json.loads((run / f'catalogue-{size}.json').read_text()),
                      key=lambda row: row['certificate'])
        assert len({r['certificate'] for r in rows}) == len(rows)
        lines, distribution = [], Counter()
        for number, row in enumerate(rows, 1):
            name = f'U{size}_{number:03d}'
            edges = row['edges']
            props = properties(17, edges)
            cert = certificate(17, edges)
            assert props == row['properties'] and cert == row['certificate']
            assert props['triangle_free'] and props['edges'] == size
            assert max(props['degrees']) <= 7 and props['alpha'] <= 8
            graph = nx.Graph()
            graph.add_nodes_from(range(17))
            graph.add_edges_from(edges)
            graph6 = nx.to_graph6_bytes(graph, header=False).decode().strip()
            obj = {'name': name, 'num_vertices': 17, 'edges': edges,
                   'properties': props, 'graph6': graph6,
                   'certificate': cert}
            write_json(output / f'{name}.json', obj)
            entry = {'name': name, 'file': f'{name}.json', 'edges': size,
                     'alpha': props['alpha'], 'factor_critical': props['factor_critical']}
            if size == 58:
                for edge in edges:
                    child = certificate(17, [e for e in edges if e != edge])
                    deletions.setdefault(child, {}).setdefault(name, 0)
                    deletions[child][name] += 1
            else:
                entry['single_edge_deletion_parents'] = deletions.get(cert, {})
            index['graphs'].append(entry)
            distribution[(props['alpha'], props['factor_critical'])] += 1
            lines.append(graph6)
        (output / f'graphs-{size}.g6').write_text('\n'.join(lines) + '\n')
        restricted = json.loads((ROOT / f'runs/full/catalogue-{size}.json').read_text())
        assert {r['certificate'] for r in rows if r['properties']['alpha'] <= 7} == {
            r['certificate'] for r in restricted}
        index['layers'][str(size)] = {
            'count': len(rows), 'distribution': [
                {'alpha': a, 'factor_critical': fc, 'count': count}
                for (a, fc), count in sorted(distribution.items())]}
    write_json(output / 'index.json', index)
    (output / 'SHA256SUMS').write_text(''.join(
        hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.name + '\n'
        for p in sorted(output.iterdir()) if p.suffix in ('.json', '.g6')))
    print(json.dumps(index['layers'], indent=2))


if __name__ == '__main__':
    main()
