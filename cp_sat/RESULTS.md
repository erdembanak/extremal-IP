Completed enumeration of the 17-vertex Earth–Moon candidate family

Date: 19 September 2026.

**Subsequent mathematical proof.** The [uniqueness proof](../uniqueness-proof.md)
now gives a self-contained argument for the 58-edge layer with the independence
restriction below, without assuming factor-criticality or using solver results.
It also proves that every 17-vertex, 58-edge triangle-free graph with maximum
degree at most seven is homomorphic to C5. This is a new proof draft, not yet
independently refereed or incorporated into the Lean audit. The 57-edge
classification below remains computational.

**Result.** Up to isomorphism, the computational enumeration finds exactly
one graph with 58 edges and three graphs with 57 edges satisfying

\[
|V(F)|=17,\quad F\text{ triangle-free},\quad
\Delta(F)\le7,\quad\alpha(F)\le7,\quad F\text{ factor-critical}.
\]

All four have independence number seven and matching number eight.
The second encoding also omits factor-criticality and finds no additional
graphs, so the enumeration supports the same classification without that
assumption.

**Explicit description.** Let B be the complete independent-set blowup of C5
with consecutive bag sizes (4,4,3,3,3). The complete list is:

| Representative | Description | Degree multiset | Edges |
|---|---|---|---|
| F58 | B | 7^14, 6^3 | 58 |
| F57_44 | B minus an edge joining the two size-four bags | 7^12, 6^5 | 57 |
| F57_43 | B minus an edge joining a size-four and size-three bag | 7^12, 6^5 | 57 |
| F57_33 | B minus an edge joining two size-three bags | 7^13, 6^3, 5 | 57 |

Superscripts denote multiplicities. In particular the two graphs with degree
multiset 7^12,6^5 are not isomorphic. The 58 choices of a deleted edge in B
fall into exactly these three classes, with respectively 16, 24, and 18 choices.

Portable graph6 files are [candidates.g6](runs/full/candidates.g6) for F and
[cones.g6](runs/full/cones.g6) for K1 joined to its complement, in the table's
order. [explicit-representatives.json](runs/full/explicit-representatives.json)
contains bag labels, deleted edges, all edge lists, graph6 strings, and checked
properties. Regenerate these with:

```bash
.tools/cp-sat-venv/bin/python cp_sat/export_catalogue.py
```

**Why the enumeration is complete in this scope.** The degree deficiency
17·7−2e is three at e=58 and five at e=57. Integer partitions of these
deficiencies give three and seven degree multisets, respectively. Every one
of these ten cases reached residual INFEASIBLE after found graphs were saved
and excluded. No case ended on a time or solution limit.

The production encoding used fixed degree sequences, triangle exclusions,
constraints excluding independent eight-sets, and perfect-matching witnesses
after each vertex deletion. It used three concurrent jobs, with three CP-SAT
workers per job. Individual case times ranged from about 2 to 20 seconds.
The aggregate [status.json](runs/full/status.json) lists all ten completed
cases; per-case directories contain source hashes, package versions, saved
assignments, and terminal statuses. Solver logs are alongside these directories.

Symmetry restrictions compare the adjacency vector with vertex transpositions
within equal-degree classes. The lexicographically least vector in any such
orbit survives all comparisons. Extra exclusion cuts remove only relabelings
of an already saved graph, and nauty identifies duplicates in the output.
The resulting class sets were checked against brute-force small-graph
enumeration before the production run. All thirteen unit tests passed.

**Second encoding.** [audit_enumeration.py](audit_enumeration.py) independently
constructs a model using nonedge variables, reverses the lexicographic
symmetry orientation, omits every matching-witness constraint, and excludes
only individual returned assignments. It uses a separately implemented degree
sequence generator. Every case again reached INFEASIBLE, with the same four
isomorphism classes and no non-factor-critical graphs. Isomorphism comparison
in this audit uses NetworkX, independently of nauty's certificates.
See [audit.json](runs/full/audit.json).

Both encodings use OR-Tools CP-SAT 9.15.6755. This is a replicated computational
classification, not a separately checked formal UNSAT proof or an independent
solver replication. The mathematical completeness argument still relies on
the encoding, the justified symmetry restrictions, and the solver's
infeasibility results.

**Conditional consequence for graph thickness.** Define G0=K1 joined to the complement of B.
If F=B−e then

\[
K_1\vee\overline F=G_0+e.
\]

Thus every candidate cone contains G0 as a spanning subgraph. Since
biplanarity is inherited by subgraphs, proving G0 non-biplanar excludes
all four candidates. Conversely, if G0 is biplanar it is itself a ten-chromatic
witness. Together with the sender's structural reduction, the 18-vertex
independence-number-two question therefore reduces to this one graph.

This does not independently establish G0's non-biplanarity. The sender's
reported computational exclusion still needs its own verification. It also
does not settle the 19-vertex case or the full Earth–Moon problem.

This catalogue has the independence restriction α≤7. It is not the
unrestricted list of all (7,8) extremizers; in particular B_(7,8), whose
independence number is eight, is intentionally absent.
