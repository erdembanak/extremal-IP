import TriangleProof.Stability
import TriangleProof.Counting

/-!
# The `d = 4t - 1` component case

This file connects the saturated-cycle, obstruction-exclusion, partition,
and capacity arguments, beginning with an initial five-cycle. The later
`InitialCycle` and `FixedOrder` modules prove and discharge that premise.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- In the last odd degree class, the counterexample deficiency bound is
incompatible with the presence of a five-cycle. -/
theorem hard_case_contradiction_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (t : ℤ) (ht : 3 ≤ t)
    (hn : (Fintype.card V : ℤ) = 10 * t - 1)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t - 1)
    (hD : deficiency G (4 * t - 1) ≤ 2 * t - 5)
    (hcycle : ∃ φ, IsCycleCopy G φ) : False := by
  have hn' : (Fintype.card V : ℤ) =
      2 * (4 * t - 1) + 2 * t + 1 := by omega
  have hD' : deficiency G (4 * t - 1) ≤
      (4 * t - 1) - 2 * t - 4 := by omega
  obtain ⟨φ, hφ, hsat⟩ := saturated_cycle_of_cycle G htriangle
    (4 * t - 1) t hn' hdeg hD' hcycle
  obtain ⟨c, hcφ, hhom⟩ := saturated_cycle_stability G htriangle φ hφ
    (4 * t - 1) hdeg hsat (Or.inr ⟨by omega, by omega, by omega⟩)
  have hsat' : ∀ j, ∃ v, c v = j ∧ degreeZ G v = 4 * t - 1 := by
    intro j
    exact ⟨φ j, hcφ j, hsat j⟩
  have hlower := five_cycle_coloring_deficiency G c hhom t ht hn hdeg hsat'
  omega

/-- Edge-count form of the hard component bound, parameterized by an initial
five-cycle; `FixedOrder.lean` discharges this premise. -/
theorem hard_case_edge_bound_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (t : ℤ) (ht : 3 ≤ t)
    (hn : (Fintype.card V : ℤ) = 10 * t - 1)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t - 1)
    (hcycle : ∃ φ, IsCycleCopy G φ) :
    edgeCountZ G ≤ (4 * t - 1) * (5 * t - 1) + t + 1 := by
  by_contra h
  have he : (4 * t - 1) * (5 * t - 1) + t + 1 < edgeCountZ G := by omega
  have hform : deficiency G (4 * t - 1) =
      (4 * t - 1) * (2 * (4 * t - 1) + 2 * t + 1) - 2 * edgeCountZ G := by
    rw [deficiency_eq G, hn]
    ring
  have he' : (4 * t - 1) * ((4 * t - 1) + t) + t + 1 < edgeCountZ G := by
    have hid : (4 * t - 1) + t = 5 * t - 1 := by ring
    rw [hid]
    exact he
  have hD' := counterexample_deficiency (4 * t - 1) t (edgeCountZ G)
    (deficiency G (4 * t - 1)) hform he'
  have hD : deficiency G (4 * t - 1) ≤ 2 * t - 5 := by omega
  exact hard_case_contradiction_from_cycle G htriangle t ht hn hdeg hD hcycle

end TriangleProof
