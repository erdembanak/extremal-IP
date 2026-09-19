# extremal-IP

Integer-programming and CP-SAT code for triangle-free graphs with bounded
maximum degree and matching number.

- `main/` and `archive/`: original C++/CPLEX implementations. Comment out
  callbacks to compare runs with and without orbital branching.
- [`cp_sat/`](cp_sat/README.md): Python/OR-Tools implementation with nauty
  canonicalization, complete enumeration, and a second encoding for auditing.
- [`cp_sat/RESULTS.md`](cp_sat/RESULTS.md): the completed 17-vertex catalogue
  relevant to the independence-number-two Earth–Moon reduction.
- [`cp_sat/data/`](cp_sat/data/README.md): individual graph JSON files,
  solver-compatible cone inputs, and a machine-readable index.

## Completed restricted catalogue

For triangle-free graphs F on 17 vertices with maximum degree at most 7,
independence number at most 7, and 57 or 58 edges, the computational enumeration
finds exactly four isomorphism classes. All four are factor-critical:

| Edges | Classes | Description |
|---|---|---|
| 58 | 1 | Independent-bag C5 blowup B with consecutive sizes (4,4,3,3,3) |
| 57 | 3 | B with one edge removed, of bag-size types 4–4, 4–3, or 3–3 |

This is **not** the unrestricted catalogue of all (7,8) extremizers: the
independence-number restriction is essential. Completion relies on CP-SAT
infeasibility results. Two encodings agree, but both use the same solver;
no separately checkable formal UNSAT proof is supplied.

## Quick check

From the repository root, using Python 3.10 and a C compiler:

```sh
python -m venv .venv
.venv/bin/python -m pip install -r cp_sat/requirements.txt
.venv/bin/python -m unittest discover -s cp_sat -v
.venv/bin/python cp_sat/verify_catalogue.py
```

The catalogue checker needs only NetworkX and does not run CP-SAT or nauty.
It verifies the individual graphs and deletion classification, not exhaustiveness
among arbitrary graphs. Exhaustiveness is supported by the recorded solver runs.
See [reproduction instructions](cp_sat/REPRODUCE.md) and the
[mathematical formulation](cp_sat/FORMULATION.md).

## Connection to graph thickness

The family is motivated by the structural reduction in
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
