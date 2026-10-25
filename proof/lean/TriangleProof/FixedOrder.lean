import TriangleProof.InitialCycle

/-!
# Fixed-order component theorem

This file removes the `C5` hypothesis from the three checked degree cases,
using the proved initial-cycle theorem.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- The complete fixed-order upper bound in every numerical regime needed
below `Z(d)`. There are no external graph-theoretic premises. -/
theorem fixed_order_edge_bound
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (d t : ℤ) (hd : 7 ≤ d) (ht : 1 ≤ t) (htd : t < d)
    (hn : (Fintype.card V : ℤ) = 2 * d + 2 * t + 1)
    (hdeg : ∀ v, degreeZ G v ≤ d)
    (hregime : FixedOrderRegime d t) :
    edgeCountZ G ≤ d * (d + t) + t + 1 := by
  by_contra hnot
  have he : d * (d + t) + t + 1 < edgeCountZ G := by omega
  have hcycle := dense_triangle_free_five_cycle G d t hd ht htd hn htriangle hdeg he
  rcases hregime with hlarge | hboundary
  · exact hnot (large_degree_edge_bound_from_cycle G htriangle d t hn hdeg hlarge hcycle)
  rcases hboundary with ⟨a, ha, hda, hta⟩ | ⟨a, ha, hda, hta⟩
  · subst d
    subst t
    have hh := hard_case_edge_bound_from_cycle G htriangle a ha (by omega) hdeg hcycle
    exact hnot (by nlinarith [hh])
  · subst d
    subst t
    have hh := one_mod_four_edge_bound_from_cycle G htriangle a ha (by omega) hdeg hcycle
    exact hnot (by nlinarith [hh])

/-- Contradiction form, convenient when assembling an extremal theorem. -/
theorem no_fixed_order_counterexample
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (d t : ℤ) (hd : 7 ≤ d) (ht : 1 ≤ t) (htd : t < d)
    (hn : (Fintype.card V : ℤ) = 2 * d + 2 * t + 1)
    (hdeg : ∀ v, degreeZ G v ≤ d)
    (hregime : FixedOrderRegime d t)
    (he : d * (d + t) + t + 1 < edgeCountZ G) : False := by
  have hub := fixed_order_edge_bound G htriangle d t hd ht htd hn hdeg hregime
  omega

end TriangleProof
