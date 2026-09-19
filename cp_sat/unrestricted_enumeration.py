"""Enumerate all 17-vertex triangle-free graphs of maximum degree 7 and size 57/58.

No independence or factor-criticality constraints are imposed."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import subprocess
import sys

from enumerate_graphs import degree_sequences, write_json


def aggregate(root):
    cases, catalogues = [], {57: {}, 58: {}}
    for edges in (58, 57):
        for index, degrees in enumerate(degree_sequences(17, 7, edges)):
            path = root / f'e{edges}-p{index}'
            status_path = path / 'summary.json'
            status = json.loads(status_path.read_text()) if status_path.exists() else {
                'complete': False, 'stop_reason': 'NOT_STARTED'}
            count, unique = 0, set()
            if (path / 'solutions.jsonl').exists():
                for line in (path / 'solutions.jsonl').read_text().splitlines():
                    try:
                        row = json.loads(line)
                    except json.JSONDecodeError:
                        # A concurrent writer may not yet have finished its last
                        # row. The worker's resume path performs strict validation.
                        continue
                    count += 1
                    unique.add(row['certificate'])
                    catalogues[edges].setdefault(row['certificate'], row)
            cases.append({'edges': edges, 'pattern': index, 'degrees': degrees,
                          'complete': status['complete'], 'stop_reason': status['stop_reason'],
                          'saved_assignments': count, 'unique_graphs': len(unique)})
    for edges, catalogue in catalogues.items():
        write_json(root / f'catalogue-{edges}.json', list(catalogue.values()))
    result = {'updated_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'n=17, triangle-free, maximum degree <=7; no alpha or factor-criticality restriction',
              'complete': all(case['complete'] for case in cases),
              'completed_cases': sum(case['complete'] for case in cases),
              'unique_graphs': {str(e): len(c) for e, c in catalogues.items()},
              'cases': cases}
    write_json(root / 'status.json', result)
    return result


def run_case(root, edges, pattern, seconds, workers, round_no):
    output = root / f'e{edges}-p{pattern}'
    cmd = [sys.executable, str(Path(__file__).with_name('enumerate_graphs.py')),
           '--edges', str(edges),
           '--degree-pattern', str(pattern), '--row-lex', '--relabel-blocks', '128',
           '--workers', str(workers), '--time-limit', str(seconds),
           '--seed', str(round_no), '--output', str(output)]
    if (output / 'manifest.json').exists():
        cmd.append('--resume')
    log = root / f'e{edges}-p{pattern}.log'
    with log.open('a') as sink:
        sink.write(f'\nROUND {round_no} at {datetime.now(timezone.utc).isoformat()}\n')
        sink.flush()
        proc = subprocess.run(cmd, stdout=sink, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f'e{edges}-p{pattern} exited {proc.returncode}; inspect {log}')
    return edges, pattern


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('cp_sat/runs/unrestricted'))
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--round-seconds', type=float, default=300)
    parser.add_argument('--max-rounds', type=int, default=0, help='0: continue until complete')
    parser.add_argument('--status-only', action='store_true')
    args = parser.parse_args()
    if args.jobs < 1 or args.workers < 1 or args.round_seconds <= 0 or args.max_rounds < 0:
        parser.error('Invalid limits')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    # Never treat an old, filtered catalogue as completion of this family.
    for manifest in root.glob('*/manifest.json'):
        spec = json.loads(manifest.read_text())['spec']
        if (spec['vertices'] != 17 or spec['degree'] != 7 or
                spec['alpha_max'] is not None or spec['factor_critical']):
            parser.error(f'Incompatible restricted checkpoint: {manifest}')
    status = aggregate(root)
    if args.status_only:
        print(json.dumps(status, indent=2))
        return
    round_no = 0
    while not status['complete'] and (not args.max_rounds or round_no < args.max_rounds):
        round_no += 1
        pending = [(c['edges'], c['pattern']) for c in status['cases'] if not c['complete']]
        budget = min(args.round_seconds * 2 ** (round_no - 1), 3600)
        print(f'Round {round_no}: {len(pending)} cases, {budget}s per case, '
              f'{args.jobs} jobs x {args.workers} workers', flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_case, root, e, p, budget, args.workers, round_no)
                       for e, p in pending]
            for future in as_completed(futures):
                e, p = future.result()
                status = aggregate(root)
                case = next(c for c in status['cases'] if c['edges'] == e and c['pattern'] == p)
                print(f'e{e}-p{p}: {case["stop_reason"]}, '
                      f'{case["unique_graphs"]} unique; '
                      f'{status["completed_cases"]}/10 cases complete', flush=True)
        status = aggregate(root)
    print(json.dumps(status, indent=2), flush=True)


if __name__ == '__main__':
    main()
