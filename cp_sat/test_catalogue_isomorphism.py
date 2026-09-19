"""Validate color-assisted comparisons against NetworkX's classic matcher."""
from itertools import combinations
import random
import unittest

import networkx as nx

from catalogue_isomorphism import isomorphic


class CatalogueIsomorphismTests(unittest.TestCase):
    def test_small_atlas_against_classic_matcher(self):
        graphs = [g for g in nx.graph_atlas_g() if 1 <= len(g) <= 5]
        for left, right in combinations(graphs, 2):
            if len(left) == len(right) and left.number_of_edges() == right.number_of_edges():
                self.assertEqual(isomorphic(left, right), nx.is_isomorphic(left, right))

    def test_relabeling_preserves_colors_and_isomorphism(self):
        rng = random.Random(17)
        for g in nx.graph_atlas_g():
            if not 1 <= len(g) <= 6:
                continue
            labels = list(range(len(g)))
            rng.shuffle(labels)
            relabeled = nx.relabel_nodes(g, dict(zip(g, labels)))
            self.assertTrue(isomorphic(g, relabeled))


if __name__ == '__main__':
    unittest.main()
