"""Export individual candidates and downstream-compatible cone inputs."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    records = json.loads((ROOT / 'runs/full/explicit-representatives.json').read_text())
    output = ROOT / 'data'
    output.mkdir(exist_ok=True)
    index = {'schema_version': 1, 'vertex_labels': 'zero-based',
             'scope': 'n=17, triangle-free, maximum degree <=7, alpha <=7, edges in {57,58}',
             'graphs': []}
    for row in records:
        name = row['name']
        graph = {'name': name, 'num_vertices': 17, 'edges': row['edges'],
                 'bags': row['bags'], 'deleted_edge': row['deleted_edge'],
                 'properties': row['properties'], 'graph6': row['graph6']}
        cone = {'name': 'cone_' + name, 'num_vertices': 18, 'apex': 17,
                'source_graph': name, 'edges_part1': row['cone_edges'],
                'edges_part2': [], 'graph6': row['cone_graph6'],
                'partition_is_certificate': False,
                'note': 'All edges are stored in part1 as solver input, not as a planar partition.'}
        for filename, value in [(name + '.json', graph), ('cone_' + name + '.json', cone)]:
            (output / filename).write_text(json.dumps(value, indent=2) + '\n')
        index['graphs'].append({'name': name, 'edges': len(row['edges']),
                                'graph': name + '.json', 'cone': 'cone_' + name + '.json'})
    (output / 'index.json').write_text(json.dumps(index, indent=2) + '\n')
    (output / 'SHA256SUMS').write_text(''.join(
        hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.name + '\n'
        for p in sorted(output.glob('*.json'))))


if __name__ == '__main__':
    main()
