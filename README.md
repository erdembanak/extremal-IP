# extremal-IP

Integer-programming and CP-SAT code for triangle-free graphs with bounded
maximum degree and matching number.

- `main/` and `archive/`: original C++/CPLEX implementations. Comment out
  callbacks to compare runs with and without orbital branching.
- [`cp_sat/`](cp_sat/README.md): Python/OR-Tools implementation with nauty
  canonicalization, complete enumeration, and a second encoding for auditing.
- [`cp_sat/UNRESTRICTED.md`](cp_sat/UNRESTRICTED.md): the complete 23-graph
  fixed-order catalogue, with reproduction instructions and audit results.
- [`cp_sat/data/unrestricted/`](cp_sat/data/unrestricted/): all 23 graphs as
  individual JSON files and graph6 files, with a machine-readable index.
- [`cp_sat/RESULTS.md`](cp_sat/RESULTS.md): the four-graph subset relevant to
  the independence-number-two Earth–Moon reduction; its graphs and
  solver-compatible cone inputs are in [`cp_sat/data/`](cp_sat/data/README.md).

## Catalogue without an independence restriction

The enumeration is complete for **triangle-free graphs on 17 vertices with
maximum degree at most 7 and exactly 57 or 58 edges**, up to isomorphism.
Neither independence number nor factor-criticality is constrained.

| Edges | Total classes | Independence number 7 | Independence number 8 |
|---|---:|---:|---:|
| 58 | **2** | 1 | 1 |
| 57 | **21** | 3 | 18 |

All 23 graphs turn out to be factor-critical and have matching number 8.
The two 58-edge graphs are B_(7,8) and the independent-bag C5 blowup with
consecutive bag sizes (4,4,3,3,3). Only six of the 21 size-57 classes arise
by deleting an edge from either size-58 graph; the other fifteen require
enumeration beyond that construction.

Both encodings exhausted all ten degree cases and agreed on the complete
catalogue. All 15 unit tests and the standalone graph checker passed.
See [results and reproduction](cp_sat/UNRESTRICTED.md). This catalogue fixes
the order at 17; it does not enumerate graphs of other orders.

## Completed restricted catalogue

For triangle-free graphs F on 17 vertices with maximum degree at most 7,
independence number at most 7, and 57 or 58 edges, the computational enumeration
finds exactly four isomorphism classes. All four are factor-critical:

| Edges | Classes | Description |
|---|---|---|
| 58 | 1 | Independent-bag C5 blowup B with consecutive sizes (4,4,3,3,3) |
| 57 | 3 | B with one edge removed, of bag-size types 4–4, 4–3, or 3–3 |

This four-graph table is the independence-at-most-seven subset of the larger
catalogue above. For both catalogues, completion relies on CP-SAT infeasibility
results. The two encodings use the same solver; no separately checkable formal
UNSAT proof is supplied.

## Quick check

From the repository root, using Python 3.10 and a C compiler:

```sh
python -m venv .venv
.venv/bin/python -m pip install -r cp_sat/requirements.txt
.venv/bin/python -m unittest discover -s cp_sat -v
.venv/bin/python cp_sat/verify_unrestricted.py
.venv/bin/python cp_sat/verify_catalogue.py
```

The two catalogue checkers need only NetworkX and do not run CP-SAT or nauty.
They verify individual graphs and deletion relationships, not exhaustiveness
among arbitrary graphs. Exhaustiveness is supported by the recorded solver runs.
See reproduction instructions for the [23-graph catalogue](cp_sat/UNRESTRICTED.md)
and [four-graph subset](cp_sat/REPRODUCE.md), and the
[mathematical formulation](cp_sat/FORMULATION.md).

## Connection to graph thickness

The four-graph restricted family is motivated by the structural reduction in
[Spogreev's earth-moon-alpha2 repository](https://github.com/ivanspog/earth-moon-alpha2),
inspected at commit `53af3680224e6ec5942e6527b964b90fc85312f4`.
Every candidate cone is G0 = K1 joined to the complement of B, or G0 with
one additional edge. Thus non-biplanarity of G0 would exclude all four.

That non-biplanarity is **reported by Spogreev, not independently established
here**. A local run without vertex-symmetry constraints hit its 300-second
limit without a verdict. The 18-vertex conclusion remains conditional on
that claim and the structural reduction. The 19-vertex case and the full
Earth–Moon problem are not settled by this catalogue.

## License

Repository code and data are released under the [MIT license](LICENSE).
External dependencies retain their respective licenses.
