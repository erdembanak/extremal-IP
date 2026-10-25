# Audit of the external results used by the Lean proof

**Integrated manuscript (2026-09-06):** The prose now proves the short
component-reduction and endpoint results within the paper, with attribution,
and proves the small threshold values via the common construction family.
The map below describes the Lean interface, whose published-result hypotheses
remain in place. Incorporating a written proof does not instantiate that field.
See the formalization section of the
[manuscript](paper/triangle-free-extremal-formula.tex) for the distinction
between manuscript and formal coverage.

The Lean development proves the new fixed-order graph argument and the final
packing argument. It represents results already established in the two source
papers by fields of `AuditedLiteratureInputs`. This is the usual formalization
boundary for a proof that cites prior work: Lean checks every deduction after
the cited statements, but it does not regard a bibliographic citation as a term
of the corresponding Lean type.

## Exact field-to-source map

| Lean field | Mathematical content | Source |
|---|---|---|
| `publishedZ.candidateExists` | For every `d >= 2`, the threshold-candidate set defining `Z(d)` is nonempty. | Ahanjideh--Ekim--Yildiz (2022), Proposition 4.3. |
| `publishedZ.lower` | `d+1 <= Z(d)` for `d >= 4`. | Ahanjideh--Ekim--Yildiz (2022), Lemma 4.6. |
| `publishedZ.even` | `Z(d)=floor(5d/4)` for even `d`. | Ahanjideh--Ekim--Yildiz (2022), Lemma 4.7. |
| `publishedZ.oddUpper` | `Z(d) <= floor(5(d+1)/4)` for odd `d`. | Ahanjideh--Ekim--Yildiz (2022), Lemma 4.7, using the construction in Proposition 4.3. |
| `leftEndpoint` | `f_Delta(d,d)=d^2+1`. | Ahanjideh--Ekim--Yildiz (2022), Theorem 3.4. |
| `thresholdEndpoint` | `f_Delta(d,Z(d))=d Z(d)+floor(d/2)`. | Ahanjideh--Ekim--Yildiz (2022), Theorem 4.8. |
| `structuralDecomposition` | An extremal graph can be reduced to `d`-stars and factor-critical extremal components with matching numbers in `[d,Z(d)]`, giving the displayed knapsack upper bound. | Ahanjideh--Ekim--Yildiz (2022), Corollary 3.5 and Lemmas 4.5 and 4.9; restated as Banak--Ekim--Taskin (2023), Lemma 2.3 and Equation (1). |
| `publishedSolvedCases` | The formula holds when `m <= d`, and for every `m` when `d <= 6`. | Ahanjideh--Ekim--Yildiz (2022), Theorem 5.1. |
| `construction22Assembly` | Given the locally established lower bound for every residual `B_(d,r)` component, disjoint unions of threshold components, one residual component when needed, and `d`-stars attain the conjectured value. | The disjoint-union assembly behind Ahanjideh--Ekim--Yildiz (2022), Proposition 4.3 and Theorem 4.8. |

The Lean theorems `belowZRegime_of_publishedZResults` and
`thresholdBounds_of_publishedZResults` derive all coarser numerical conditions
used by the new proof. Those conditions are therefore not extra assumptions.
Likewise, `knapsackDecomposition_of_structural` and
`packingCertificate_of_fixed_components` prove the translation from the
published component statement to the numerical packing theorem.

## Proposition 4.1 is now proved locally

`IntermediateConstruction.lean` defines `B_(d,d+t)` from its five vertex
classes and uses cyclic shifts `1,...,t` for the deleted perfect matchings.
Lean proves the vertex count, triangle-freeness, every degree, the edge count,
and matching number `d+t` using an explicit matching involution. It then
relabels the graph onto `Fin (2(d+t)+1)` and proves the lower bound required by
Conjecture 2.1. The former `AuditedLiteratureInputs.construction21` field has
therefore been removed. The `construction22Assembly` field receives that
locally proved lower bound as an argument and records only the separate
disjoint-union realization of the complete Conjecture 2.2 formula.

## Degree seven is now proved locally

The former `publishedZ.seven` field has been removed. In
`SevenConstruction.lean`, `Z_seven_le_nine` proves `Z(7)≤9` without external
hypotheses, using an explicitly defined B(7,9). Lean checks triangle-freeness,
all degrees, matching number nine, and perfect matchings after every vertex
deletion. The manuscript gives a short Hall argument for the latter property.

`seven_not_candidate_eight` uses the existing fixed-order bound on seventeen
vertices. `Z_seven_eq_nine` combines this with the general lower bound
`8≤Z(7)` and the explicit construction. This general numerical lower bound is
its only remaining premise; the initial-cycle argument is now proved locally.
The two arithmetic regime/packing lemmas now use `Z_seven_le_nine` directly.
The computation reported in the 2023 paper is no longer an input to the proof.

## Results deliberately not imported

Banak--Ekim--Taskin Proposition 2.1, and the corresponding proposition in the
2022 paper, are conditional on Conjecture 2.1. The Lean proof does not use
them. Its packing normalization is proved directly in `Packing.lean`, avoiding
any circular appeal to that conditional proposition.

The phrase in Banak--Ekim--Taskin Proposition 4.1 saying that the graph is
extremal if Conjecture 2.1 holds is not used as an upper bound. Lean uses only
the explicitly proved graph properties.

## Source discrepancies

The 2023 PDF says `Z(9)=13` once in its introductory discussion, but Section 4
and the computational table give `Z(9)=12`. The latter is consistent with the
2022 theorem `Z(9) <= floor(5(9+1)/4)=12`. The formalization never uses the
value `13`.

Table 5 also gives 148 edges at `(d,m)=(8,18)`, while Table 6 and the
conjectured formula give 149. Table 6 is consistent: its one threshold
component and one matching-number-eight component have 149 edges in total.

## The standard odd-cycle input is now proved locally

`ShortestOddCycle.lean` proves that a non-bipartite graph has a shortest odd
cycle, by minimizing odd closed walks and excluding repeated vertices.
Minimality also proves the two-neighbor bound required by the initial-cycle
argument. `InitialCycle.lean` consequently proves the full initial-cycle
principle. The former `standardInitialCycle` field is deleted from both
literature interfaces, and the fixed-order edge bound no longer has a cycle
existence premise.

## Remaining formal boundary

An unconditional Lean theorem now requires a value of
`AuditedLiteratureInputs`. Supplying it would mean formalizing the cited
papers' component-reduction, endpoint, remaining threshold, disjoint-union,
and solved-case proofs.
For a conventional paper proof these are
legitimate citations. For a foundational, single-kernel Lean proof they
remain explicit hypotheses, rather than hidden axioms or admitted lemmas.
