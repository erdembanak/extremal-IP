# Enumeration without an independence-number constraint

The scope here is **17 vertices**, triangle-free, maximum degree at most 7,
and exactly 57 or 58 edges. Neither independence number nor factor-criticality
is constrained. The matching upper bound of 8 follows from the order.

The completed production enumeration finds:

| Edges | Total classes | Independence 7 | Independence 8 |
|---|---:|---:|---:|
| 58 | 2 | 1 | 1 |
| 57 | 21 | 3 | 18 |

Every returned graph is factor-critical and has matching number 8. These
properties were checked after solving, not imposed on the model.

Independence number cannot exceed 8: an independent set of at least 9 leaves
at most 8 vertices meeting every edge, so the degree bound gives e ≤ 8·7=56.
The alpha-at-most-seven classes agree exactly with the earlier restricted
catalogue. The additional 58-edge class is isomorphic to B_(7,8), with degree
multiset 7^16,4; the other is the (4,4,3,3,3) independent-bag C5 blowup.

Only **six** of the 21 size-57 classes arise by deleting one edge from one of
the two size-58 graphs (three from each). The other **fifteen** require genuine
enumeration. The index records parent classes and deletion multiplicities.

## Files and commands

- [`data/unrestricted/index.json`](data/unrestricted/index.json): all 23 classes,
  independence numbers, factor-criticality, and edge-deletion relationships.
- [`data/unrestricted/`](data/unrestricted/): individual JSON files,
  `graphs-57.g6`, `graphs-58.g6`, and SHA-256 checksums.
- [`runs/unrestricted/`](runs/unrestricted/): manifests, solver logs,
  assignments, completion summaries, and audit records.

From the repository root, with `cp_sat/requirements.txt` installed:

```sh
python cp_sat/unrestricted_enumeration.py --round-seconds 120
python cp_sat/audit_unrestricted.py
python cp_sat/export_unrestricted.py
python cp_sat/verify_unrestricted.py
```

The supervisor resumes existing unfinished cases; completed cases are skipped.
To reproduce from scratch, first move the supplied `cp_sat/runs/unrestricted`
directory aside. Do not point this supervisor at the restricted `runs/full`
directory: incompatible filtered manifests are rejected.

The second encoding uses nonedge variables, opposite lexicographic orientation,
no independence/matching constraints, and only single-assignment exclusions.
It uses the same CP-SAT solver, with NetworkX for isomorphism comparisons.
Those comparisons use invariant vertex-color refinement followed by exact
VF2++ matching; no nauty certificates are used to decide audit matches.
Check `runs/unrestricted/audit.json` for `complete_agreement: true` before
claiming a completed replication. Its per-case time limit is 300 seconds.

**Audit result (19 September 2026): complete agreement.** All ten cases in
the second encoding reached residual INFEASIBLE, recovering exactly the same
23 isomorphism classes and no non-factor-critical graphs. The final audit
log is [`runs/unrestricted/audit.log`](runs/unrestricted/audit.log).

Completion remains a computational claim dependent on the encoding and solver;
neither run emits a separately checkable formal UNSAT proof. The standalone
catalogue checker uses only NetworkX and checks properties, distinctness, and
deletion relationships; it does not prove exhaustive coverage.

All 15 unit tests passed, including comparisons of the color-assisted matcher
against NetworkX's classic matcher on small atlas graphs, and invariance under
vertex relabeling. The standalone checker verified all 23 graph files,
independence and matching properties, graph6 encodings, checksums, and all
single-edge-deletion relationships. All production source hashes match the
saved manifests, and all ten production cases reached residual INFEASIBLE.

This catalogue does not enumerate graphs of other orders. It also does not
change the earlier four-candidate thickness reduction, which explicitly
requires independence number at most 7 in F.
