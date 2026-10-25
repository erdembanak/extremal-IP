import TriangleProof.SaturatedCycle
import TriangleProof.Counting
import TriangleProof.Stability

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- Once a five-cycle is present, the saturated-cycle argument immediately
settles every parameter range with `d ≥ 4t+3`. -/
theorem large_degree_contradiction_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (d t : ℤ) (hn : (Fintype.card V : ℤ) = 2 * d + 2 * t + 1)
    (hdeg : ∀ v, degreeZ G v ≤ d)
    (hD : deficiency G d ≤ d - 2 * t - 4)
    (hd : 4 * t + 3 ≤ d)
    (hcycle : ∃ φ, IsCycleCopy G φ) : False := by
  obtain ⟨φ, hφ, hsat⟩ := saturated_cycle_of_cycle G htriangle d t hn hdeg hD hcycle
  let N := fun v => Finset.univ.filter (fun i : CycleVertex => G.Adj v (φ i))
  have hN : ∀ v, (N v).card ≤ 2 := by
    intro v
    exact cycle_independent_bound _ (by
      intro i hi j hj hij
      exact htriangle v (φ i) (φ j)
        ⟨(Finset.mem_filter.mp hi).2, hφ.2 i j hij,
          ((Finset.mem_filter.mp hj).2).symm⟩)
  have hcount : (∑ j, degreeZ G (φ j)) = ∑ v, ((N v).card : ℤ) := by
    rw [sum_degreeZ_on_copy]
    apply Finset.sum_congr rfl
    intro v _
    simp [N, incidence]
  have hupper : (∑ v, ((N v).card : ℤ)) ≤ 2 * (Fintype.card V : ℤ) := by
    have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hN v)
    norm_num only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hh
    have hh' : (∑ v, ((N v).card : ℤ)) ≤ (Fintype.card V : ℤ) * 2 := by
      exact_mod_cast hh
    simpa [mul_comm] using hh'
  simp only [hsat, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul] at hcount
  norm_num at hcount
  omega

theorem large_degree_edge_bound_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (d t : ℤ) (hn : (Fintype.card V : ℤ) = 2 * d + 2 * t + 1)
    (hdeg : ∀ v, degreeZ G v ≤ d)
    (hd : 4 * t + 3 ≤ d)
    (hcycle : ∃ φ, IsCycleCopy G φ) :
    edgeCountZ G ≤ d * (d + t) + t + 1 := by
  by_contra h
  have he : d * (d + t) + t + 1 < edgeCountZ G := by omega
  have hD : deficiency G d ≤ d - 2 * t - 4 := by
    apply counterexample_deficiency d t (edgeCountZ G) (deficiency G d)
    · rw [deficiency_eq G, hn]
    · exact he
  exact large_degree_contradiction_from_cycle G htriangle d t hn hdeg hD hd hcycle

/-- The `d = 4t+1` case follows from the same obstruction-driven coloring
used in the hard class.  Here the excess-potential sum is one. -/
theorem one_mod_four_contradiction_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (t : ℤ) (ht : 2 ≤ t)
    (hn : (Fintype.card V : ℤ) = 10 * t + 3)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t + 1)
    (hD : deficiency G (4 * t + 1) ≤ 2 * t - 3)
    (hcycle : ∃ φ, IsCycleCopy G φ) : False := by
  have hn' : (Fintype.card V : ℤ) =
      2 * (4 * t + 1) + 2 * t + 1 := by omega
  have hD' : deficiency G (4 * t + 1) ≤
      (4 * t + 1) - 2 * t - 4 := by omega
  obtain ⟨φ, hφ, hsat⟩ := saturated_cycle_of_cycle G htriangle
    (4 * t + 1) t hn' hdeg hD' hcycle
  obtain ⟨c, hcφ, hhom⟩ := saturated_cycle_stability G htriangle φ hφ
    (4 * t + 1) hdeg hsat (Or.inl (by omega))
  have hsat' : ∀ j, ∃ v, c v = j ∧ degreeZ G v = 4 * t + 1 := by
    intro j
    exact ⟨φ j, hcφ j, hsat j⟩
  have hlower := five_cycle_coloring_deficiency_one G c hhom t (by omega) hn hdeg hsat'
  omega

theorem one_mod_four_edge_bound_from_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (t : ℤ) (ht : 2 ≤ t)
    (hn : (Fintype.card V : ℤ) = 10 * t + 3)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t + 1)
    (hcycle : ∃ φ, IsCycleCopy G φ) :
    edgeCountZ G ≤ (4 * t + 1) * (5 * t + 1) + t + 1 := by
  by_contra h
  have he : (4 * t + 1) * (5 * t + 1) + t + 1 < edgeCountZ G := by omega
  have he' : (4 * t + 1) * ((4 * t + 1) + t) + t + 1 < edgeCountZ G := by
    have hid : (4 * t + 1) + t = 5 * t + 1 := by ring
    rw [hid]
    exact he
  have hform : deficiency G (4 * t + 1) =
      (4 * t + 1) * (2 * (4 * t + 1) + 2 * t + 1) - 2 * edgeCountZ G := by
    rw [deficiency_eq G, hn]
    ring
  have hD' := counterexample_deficiency (4 * t + 1) t (edgeCountZ G)
    (deficiency G (4 * t + 1)) hform he'
  have hD : deficiency G (4 * t + 1) ≤ 2 * t - 3 := by omega
  exact one_mod_four_contradiction_from_cycle G htriangle t ht hn hdeg hD hcycle

end TriangleProof
