# Exact extremal formula: proof and Lean verification

A proof of the general formula for f_Δ(d,m), the maximum number of edges in
a triangle-free graph with maximum degree at most d and matching number at
most m, for all d ≥ 2 and m ≥ 0. It resolves Conjectures 2.1 and 2.2 of
Banak–Ekim–Taşkın (2023) and determines the threshold Z(d) for every d,
including Conjecture 6.5 of Ahanjideh–Ekim–Yıldız (2022):
Z(d) = ⌊5(d+1)/4⌋ for every odd d ≥ 9, and Z(7) = 9.

This is a research draft and has not been peer reviewed.

- [`paper/`](paper/): the manuscript
  ([PDF](paper/triangle-free-extremal-formula.pdf),
  [TeX](paper/triangle-free-extremal-formula.tex)).
- [`lean/`](lean/README.md): Lean 4 / Mathlib formalization. It proves the new
  fixed-order bound, the shortest-odd-cycle step, the packing argument, the
  Proposition 4.1 construction family, Z(7) = 9, and the odd threshold
  (`TriangleProof/OddThreshold.lean`). The final theorems take the cited
  published results as explicit hypotheses.
- [`literature-audit.md`](literature-audit.md): each Lean hypothesis and the
  published theorem it corresponds to.
- [`verify_proof.py`](verify_proof.py): standard-library finite checks that
  supplement the written proof (small graphs, capacity patterns, counting
  identities, the knapsack packing, and the construction family with every
  one-vertex deletion).

## Checking

```bash
python3 verify_proof.py

cd lean
lake exe cache get
lake build
lake env lean Audit.lean
```

On 2026-10-25 the build completed (3081 jobs), and all 165 audited theorems
depend only on `propext`, `Classical.choice`, and `Quot.sound`.
