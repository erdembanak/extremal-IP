# Lean verification of the proof draft

This formalization checks the fixed-order and packing arguments in
[the manuscript](../paper/triangle-free-extremal-formula.tex). It contains final theorems for
both conjectures conditional on an explicit structure of published inputs.
The standard initial-cycle principle and the Proposition 4.1 construction are
now proved locally. This is not yet an
unconditional, end-to-end formal proof of the remaining external results.

The integrated manuscript now leads with an explicit formula using `z_d` and
reproduces the short 2022 component and endpoint proofs. The current Lean
interface still takes those published results as hypotheses, and its final
formula uses `Z(d)`. The full identification `Z(d)=z_d` and the generalized
Hall argument for the construction family are written proofs; the degree-seven
case, the odd-degree lower bound `odd_Z_lower` (from `d+1 <= Z(d)` alone), and
the complete intermediate construction family are additionally checked in Lean. This manuscript reorganization does not
claim that the remaining formal inputs have been discharged.

**Checked status (2026-10-25, odd threshold formalized):** all 165 audited
theorems build successfully (3081 build jobs).
Their axiom audit reports only Lean's standard foundational axioms
`propext`, `Classical.choice`, and `Quot.sound` (some theorems use fewer).
There are no admitted proofs, project axioms, or native-evaluation axioms
in these theorem dependencies.

Lean and Mathlib are pinned to release `v4.33.0`. The generated
`lake-manifest.json` pins the exact dependency commits.
The Mathlib commit is `db584cd6d46c92f209a44c0f1c829460d327499d`.

With [elan](https://github.com/leanprover/elan) installed, from this directory:

```bash
lake exe cache get        # download prebuilt Mathlib; requires network access
lake build                # build the development
lake env lean Audit.lean  # print the axiom dependencies of every audited theorem
```

## Formalization scope

`TriangleProof/Obstructions.lean` states facts about the eight-vertex graphs
`H` and `M`: triangle-freeness, bounds on independent sets, the description
of independent triples in `M`, and its neighborhood-sum identity.
It also lifts the independence bounds to copies of these graphs inside
arbitrary triangle-free graphs.

`TriangleProof/Counting.lean` states general arithmetic implications from
the draft, and the finite-set pigeonhole lemma for finding a saturated vertex.
Their counting assumptions are explicit hypotheses. In particular,
an arithmetic implication about a saturated cycle is not itself a formal
proof that the requisite cycle exists in an arbitrary graph.

`TriangleProof/Capacity.lean` states the five-cycle capacity inequality for
every integer `t>=3` and all nonnegative integer five-tuples summing to three.
Its finite case split is exhaustive in the tuple entries; `t` remains arbitrary.
It also checks the analogous potential-sum-one inequality used for degrees
congruent to one modulo four.
The revised manuscript replaces the table with a weighted bipartition proof.
`odd_excess_weight_imbalance` formalizes that proof for arbitrary odd excess
whose positive support is bipartite, including the cross-term cancellation
and quadratic bound. The production proofs for excesses one and three now use
this theorem. The earlier finite pattern theorems remain as independent
regression checks.

`TriangleProof/GraphCounting.lean` connects the integer incidence model to
Mathlib's graph degree and edge count, proves the deficiency identity, and
provides the general finite-sum and independent-weight lemmas used later.

`TriangleProof/MobiusExclusion.lean` and
`TriangleProof/ObstructionExclusion.lean` prove the two eight-vertex
obstruction exclusions inside arbitrary finite triangle-free graphs. The
Möbius-ladder equality case, its induced partition, and the resulting
contradiction are included rather than assumed.

`TriangleProof/SaturatedCycle.lean` proves that any five-cycle in a
hypothetical counterexample can be chosen with all five vertices saturated.
`TriangleProof/CyclePartition.lean` and `TriangleProof/PartitionBounds.lean`
then construct the homomorphism to `C5` forced by the obstruction and degree
bounds.

`TriangleProof/Stability.lean` packages a general saturated-cycle stability
theorem for arbitrary order `n` and degree `d`: writing `T=8d-3n` and
`q=2n-5d`, the cycle coloring extends if `D<T`, or if `D≤T`, `T>0`, and
`q≥2`. Its equality case uses the generalized
`mobius_copy_rigidity_general` in `MobiusExclusion.lean`.

`TriangleProof/HardCase.lean` connects these results to prove the full
`d=4t-1` fixed-order edge bound, conditional on the existence of the initial
five-cycle. `TriangleProof/DegreeCases.lean` proves the same conditional
fixed-order bound for `d>=4t+3` and for `d=4t+1`. Both positive-excess cases
now invoke the same general stability theorem, as in the revised prose draft.

`TriangleProof/Prerequisites.lean` gives the initial-cycle principle an exact
type. `TriangleProof/FixedOrder.lean` uses its local proof to remove the `C5`
premise from all three degree cases, with no external graph-theoretic input.

`TriangleProof/ShortestOddCycle.lean` uses Mathlib's even-closed-walk
characterization of bipartiteness to obtain a shortest odd closed walk.
Splitting at a repeated vertex proves it is a genuine cycle. Applying the
same minimality to the two arcs between neighbors proves that each graph
vertex meets a shortest odd cycle of length at least five at most twice.
The cycle-existence theorem does not even require finiteness.

`TriangleProof/InitialCycle.lean` proves the bipartite edge bound from
Mathlib's bipartition degree-sum theorem and verifies the complete deficiency
calculation excluding every shortest odd cycle of length at least seven.
It proves `shortest_odd_cycle_alternative`, `dense_triangle_free_five_cycle`,
and `initial_cycle_principle`, with no external inputs.

`TriangleProof/Packing.lean` proves termination and the utility bound for the
full component-packing normalization, for an arbitrary list of component
weights. This replaces the earlier finite dynamic-programming check with a
general proof.

`TriangleProof/ExtremalDefinitions.lean` defines matching number,
factor-criticality, the degree profile used for `Z(d)`, `Z`, `fDelta`, and the
two conjectures. `TriangleProof/ExtremalFoundations.lean` proves that the
matching-size and attainable-edge-count sets are bounded, that `fDelta` is an
attained maximum, and that a published threshold witness makes `Z(d)` an
actual least candidate. `TriangleProof/LiteratureInputs.lean` states the cited
inputs as fields, proves Conjecture 2.1 from those fields and the checked
fixed-order argument, and proves Conjecture 2.2 using the checked packing
theorem.

`TriangleProof/SevenConstruction.lean` defines B(7,9) explicitly and proves its
triangle-freeness, degree profile, factor-criticality, and matching number.
Nineteen explicit partner functions supply the perfect matchings after vertex
deletion; their validity is checked by `decide`, with no external search
certificate or native-evaluation axiom. `Z_seven_le_nine` has no external
hypotheses. The exact equality `Z_seven_eq_nine` additionally uses the general
lower bound `8≤Z(7)`; the exclusion of order seventeen has no external hypotheses.
The manuscript uses Hall's theorem for factor-criticality instead of listing
the matching partners. The former `PublishedZResults.seven` field is removed.

`TriangleProof/IntermediateConstruction.lean` defines the full family
`B_(d,d+t)` from Proposition 4.1 using cyclic deleted matchings. Lean proves
its order, triangle-freeness, degree sequence, matching number, and edge count,
then transports the graph to a canonical `Fin` type and derives the required
lower bound for `fDelta`. Consequently `AuditedLiteratureInputs` no longer has
a `construction21` field.

`TriangleProof/AuditedLiterature.lean` replaces the original broad interface
with component-level statements aligned with named results in the papers. It
formally derives the needed parity regimes and threshold inequalities from the
published bounds on `Z(d)`, translates the structural component list into the
packing certificate, and derives both conjectures from
`AuditedLiteratureInputs`. See [literature-audit.md](../literature-audit.md) for
the exact field-to-source map and the data discrepancies in the 2023 PDF.

`Audit.lean` prints each theorem's axiom dependencies. There are no intended
project axioms, admitted lemmas, or assumed versions of the conjectures.
The small graph decisions use `decide`, not `native_decide`.

## Remaining external obligations

Construct `AuditedLiteratureInputs` by formalizing the cited component
reduction, endpoint results, remaining `Z(d)` results, the disjoint-union
assembly for Conjecture 2.2, and the previously solved cases. The
arithmetic consequences and their connection to the new proof are already
checked. The former `standardInitialCycle` field has been deleted from both
literature interfaces: shortest-cycle extraction and incidence are now theorems.

The $\mathbb{N}$-valued supremum and infimum now have their side conditions in
Lean. For `fDelta`, Lean constructs a maximal matching, proves that its
endpoints form a vertex cover, obtains a uniform edge bound, proves the
attainable set nonempty, and concludes that `fDelta` is attained. For `Z`, the
external boundary is stated precisely as `publishedZ.candidateExists`, sourced
to Proposition 4.3 of the 2022 paper; Lean then proves that the candidate set
is nonempty and that `Z(d)` itself satisfies `IsZCandidate` for `d >= 2`.

The narrow final theorems are `conjecture21_of_audited_literature` and
`conjecture22_of_audited_literature`. Their hypotheses make every remaining
external obligation visible; neither conjecture nor any new central lemma is
declared as an axiom. The older broad-interface theorems remain as internal
compatibility layers.
