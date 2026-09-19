"""Compare CP-SAT enumeration with exhaustive small graphs, not with itself."""
import contextlib
import io
from itertools import combinations
import json
from pathlib import Path
import tempfile
import unittest

import networkx as nx
from ortools.sat.python import cp_model

from enumerate_graphs import Spec, build_model, certificate, degree_sequences, properties, relabelled_copies, run


def brute_force(spec):
    result = set()
    for edges in combinations(combinations(range(spec.vertices), 2), spec.edges):
        graph = nx.Graph()
        graph.add_nodes_from(range(spec.vertices))
        graph.add_edges_from(edges)
        if max(dict(graph.degree()).values()) > spec.degree:
            continue
        if any(nx.triangles(graph).values()):
            continue
        if spec.alpha_max is not None:
            if any(graph.subgraph(s).number_of_edges() == 0
                   for s in combinations(range(spec.vertices), spec.alpha_max + 1)):
                continue
        if spec.factor_critical:
            if not all(nx.is_perfect_matching(graph.subgraph(set(graph) - {v}),
                                             nx.max_weight_matching(
                                                 graph.subgraph(set(graph) - {v}),
                                                 maxcardinality=True)) for v in graph):
                continue
        result.add(edges)
    return result


class EnumerationTests(unittest.TestCase):
    def execute(self, spec, path, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return run(spec, path, seconds=30, **kwargs)

    def rows(self, path):
        return [json.loads(line) for line in (path / 'solutions.jsonl').read_text().splitlines()]

    def test_labelled_enumeration_matches_brute_force(self):
        spec = Spec(5, 3, 4, 3, degree_order=False)
        expected = brute_force(spec)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = self.execute(spec, path)
            actual = {tuple(tuple(e) for e in row['edges']) for row in self.rows(path)}
        self.assertTrue(summary['complete'])
        self.assertEqual(actual, expected)

    def test_symmetry_and_isomorphic_exclusions_preserve_classes(self):
        spec = Spec(5, 3, 4, 3)
        expected = {certificate(5, edges) for edges in brute_force(spec)}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = self.execute(spec, path, relabel_blocks=64)
            actual = {row['certificate'] for row in self.rows(path)}
        self.assertTrue(summary['complete'])
        self.assertEqual(actual, expected)

    def test_row_lex_preserves_all_small_isomorphism_classes(self):
        spec = Spec(6, 3, 6, 3, row_lex=True)
        expected = {certificate(6, edges) for edges in brute_force(spec)}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = self.execute(spec, path, relabel_blocks=64)
            actual = {row['certificate'] for row in self.rows(path)}
        self.assertTrue(summary['complete'])
        self.assertEqual(actual, expected)

    def test_fixed_degree_row_lex_preserves_classes(self):
        expected = {certificate(6, edges) for edges in brute_force(Spec(6, 3, 6, 3))}
        actual = set()
        for degrees in degree_sequences(6, 3, 6):
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)
                spec = Spec(6, 3, 6, 3, degrees=degrees, row_lex=True)
                summary = self.execute(spec, path)
                self.assertTrue(summary['complete'])
                actual.update(row['certificate'] for row in self.rows(path))
        self.assertEqual(actual, expected)

    def test_factor_critical_witnesses_do_not_duplicate_graphs(self):
        spec = Spec(5, 2, 5, 2, factor_critical=True, degree_order=False)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = self.execute(spec, path)
        self.assertTrue(summary['complete'])
        self.assertEqual(summary['labelled_solutions'], 12)  # The labelled five-cycles.
        self.assertEqual(summary['nonisomorphic_graphs'], 1)

    def test_resume_and_partial_status(self):
        spec = Spec(5, 2, 5, 2, factor_critical=True)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            first = self.execute(spec, path, max_solutions=1, relabel_blocks=64)
            self.assertFalse(first['complete'])
            self.assertEqual(first['stop_reason'], 'SOLUTION_LIMIT')
            second = self.execute(spec, path, resume=True, relabel_blocks=64)
            self.assertTrue(second['complete'])
            self.assertEqual(second['nonisomorphic_graphs'], 1)

    def test_infeasible_and_empty_graph(self):
        for spec, expected_count in [(Spec(3, 2, 3), 0), (Spec(3, 2, 0), 1)]:
            with tempfile.TemporaryDirectory() as tmp:
                summary = self.execute(spec, Path(tmp))
                self.assertTrue(summary['complete'])
                self.assertEqual(summary['labelled_solutions'], expected_count)

    def test_optimization_mantel_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = self.execute(Spec(5, 4, None), path)
            self.assertTrue(summary['complete'])
            self.assertEqual(summary['objective_upper_bound'], 6)
            self.assertEqual(self.rows(path)[0]['properties']['edges'], 6)

    def test_time_limit_is_not_completeness(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = run(Spec(5, 2, 5), tmp, seconds=1e-9)
            self.assertFalse(summary['complete'])
            self.assertEqual(summary['stop_reason'], 'TIME_LIMIT')

    def test_patterns_cover_every_small_degree_multiset(self):
        spec = Spec(5, 3, 4)
        actual = {tuple(sorted(properties(5, edges)['degrees'], reverse=True))
                  for edges in brute_force(spec)}
        self.assertTrue(actual <= set(degree_sequences(5, 3, 4)))
        self.assertEqual(len(degree_sequences(17, 7, 58)), 3)
        self.assertEqual(len(degree_sequences(17, 7, 57)), 7)

    def test_nauty_and_relabelled_cuts(self):
        graph = nx.cycle_graph(5)
        edges = tuple(sorted(tuple(sorted(e)) for e in graph.edges()))
        copies = list(relabelled_copies(5, edges, 100))
        self.assertEqual(len(copies), 12)
        self.assertEqual(len({certificate(5, e) for e in copies}), 1)
        self.assertNotEqual(certificate(5, edges), certificate(5, tuple(nx.path_graph(5).edges())))

    def test_corrupt_checkpoint_rejected(self):
        spec = Spec(5, 2, 5)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.execute(spec, path, max_solutions=1)
            with (path / 'solutions.jsonl').open('a') as sink:
                sink.write('{')
            with self.assertRaises(ValueError):
                self.execute(spec, path, resume=True)

    def test_seventeen_vertex_positive_and_negative_controls(self):
        sizes = (4, 4, 3, 3, 3)
        bags, start = [], 0
        for size in sizes:
            bags.append(range(start, start + size))
            start += size
        blowup = {tuple(sorted((u, v))) for i in range(5)
                  for u in bags[i] for v in bags[(i + 1) % 5]}
        construction = {(u, v) for u in range(8) for v in range(8, 16)
                        if not (u >= 6 and v >= 14) and not (u < 6 and v == u + 8)}
        construction.update((v, 16) for v in (6, 7, 14, 15))
        for edges, alpha, feasible in [(blowup, 7, True), (construction, 7, False),
                                        (construction, None, True)]:
            spec = Spec(17, 7, 58, alpha, factor_critical=True, degree_order=False)
            model, x = build_model(spec)
            for edge, var in x.items():
                model.add(var == (edge in edges))
            solver = cp_model.CpSolver()
            solver.parameters.num_search_workers = 1
            solver.parameters.max_time_in_seconds = 10
            status = solver.solve(model)
            self.assertEqual(status, cp_model.OPTIMAL if feasible else cp_model.INFEASIBLE)


if __name__ == '__main__':
    unittest.main()
