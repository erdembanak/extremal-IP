"""Exact NetworkX isomorphism accelerated by invariant vertex colors."""
import networkx as nx


def isomorphic(left, right):
    if len(left) != len(right) or left.number_of_edges() != right.number_of_edges():
        return False
    # Refinement colors are isomorphism invariants. Requiring them preserves
    # every actual isomorphism. Hash collisions only weaken this preprocessing;
    # VF2++ still checks a complete adjacency-preserving bijection.
    copies = []
    for graph in (left, right):
        graph = graph.copy()
        colors = nx.weisfeiler_lehman_subgraph_hashes(graph, iterations=5)
        nx.set_node_attributes(graph, {v: values[-1] for v, values in colors.items()},
                               '_catalogue_color')
        copies.append(graph)
    return nx.vf2pp_is_isomorphic(*copies, node_label='_catalogue_color')
