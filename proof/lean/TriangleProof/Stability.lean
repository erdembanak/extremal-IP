import TriangleProof.PartitionBounds

/-! A saturated-cycle stability theorem with arbitrary order and degree.
The strict incidence-budget case and its rigid equality boundary share the
same obstruction exclusion and coloring extension. -/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

theorem saturated_cycle_stability
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d)
    (hsat : ∀ j, degreeZ G (φ j) = d)
    (hbudget : deficiency G d < 8 * d - 3 * (Fintype.card V : ℤ) ∨
      (deficiency G d ≤ 8 * d - 3 * (Fintype.card V : ℤ) ∧
       0 < 8 * d - 3 * (Fintype.card V : ℤ) ∧
       2 ≤ 2 * (Fintype.card V : ℤ) - 5 * d)) :
    ∃ c : V → CycleVertex, (∀ j, c (φ j) = j) ∧
      (∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w)) := by
  have hD : deficiency G d ≤ 8 * d - 3 * (Fintype.card V : ℤ) := by
    rcases hbudget with h | ⟨h, _, _⟩ <;> omega
  have hnoM : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, mobius.Adj i j → G.Adj (f i) (f j)) := by
    intro f hf hE
    rcases hbudget with h | ⟨h, hpos, hq⟩
    · have := mobius_copy_deficiency_lower_bound G htriangle f hf hE d hdeg
      omega
    · exact mobius_copy_rigidity_general G htriangle f hf hE d hdeg hq hpos h
  have hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)) := by
    intro f hf hE
    have := obstruction_copy_deficiency_lower_bound G htriangle hnoM f hf hE d hdeg
    omega
  exact saturated_cycle_coloring G htriangle hnoH φ hφ d hdeg hsat (by omega)

end TriangleProof
