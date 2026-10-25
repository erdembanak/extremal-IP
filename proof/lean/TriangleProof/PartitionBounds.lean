import TriangleProof.CyclePartition

namespace TriangleProof

set_option maxHeartbeats 0

theorem cycle_bag_upper_pattern : ∀ (s : Finset CycleVertex) (j : CycleVertex),
    cycleIndependent s →
    (if s = fiveCycle.neighborFinset j then (1 : ℤ) else 0) +
      (if j + 2 ∈ s then 1 else 0) + (if j + 3 ∈ s then 1 else 0) ≤ 1 := by decide

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

theorem saturated_cycle_partition_bounds
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (d : ℤ) (hsat : ∀ j, degreeZ G (φ j) = d) :
    ((cycleRemainder G φ).card : ℤ) ≤ 2 * (Fintype.card V : ℤ) - 5 * d ∧
    ∀ j, ((cycleBag G φ j).card : ℤ) ≤ (Fintype.card V : ℤ) - 2 * d := by
  classical
  let N := cycleNeighbors G φ
  have hi : ∀ v, cycleIndependent (N v) := cycleNeighbors_independent G htriangle φ hφ
  have hdeg : ∀ k, (∑ v, if k ∈ N v then (1 : ℤ) else 0) = degreeZ G (φ k) := by
    intro k
    simp only [N, cycleNeighbors, Finset.mem_filter, Finset.mem_univ, true_and]
    simp only [degreeZ, incidence, G.adj_comm]
  constructor
  · have hp : ∀ v, (if v ∈ cycleRemainder G φ then (1 : ℤ) else 0) ≤ 2 - (N v).card := by
      intro v
      have hc := cycle_independent_bound _ (hi v)
      by_cases hv : v ∈ cycleRemainder G φ
      · have hh := (Finset.mem_filter.mp hv).2
        change (N v).card < 2 at hh
        simp only [if_pos hv]
        omega
      · simp only [if_neg hv]
        omega
    have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hp v)
    have hsum : (∑ v, ((N v).card : ℤ)) = 5 * d := by
      have hc := sum_degreeZ_on_copy G φ
      simp only [hsat] at hc
      norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hc
      have heq : ∀ v, (∑ i, incidence G v (φ i)) = ((N v).card : ℤ) := by
        intro v
        simp [N, cycleNeighbors, incidence]
      simpa only [heq] using hc.symm
    simp only [Finset.sum_sub_distrib, hsum] at hh
    simpa [mul_comm] using hh
  · intro j
    have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => cycle_bag_upper_pattern (N v) j (hi v))
    simp only [Finset.sum_add_distrib, hdeg, hsat] at hh
    have hb : (∑ v, if N v = fiveCycle.neighborFinset j then (1 : ℤ) else 0) = (cycleBag G φ j).card := by
      simp [cycleBag, N]
    rw [hb] at hh
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one] at hh
    omega

theorem neighbors_outside_set_lower_bound (R : Finset V) (x : V) (hx : x ∈ R) :
    degreeZ G x - (R.card : ℤ) + 1 ≤
      ((Finset.univ.filter (fun v => G.Adj x v ∧ v ∉ R)).card : ℤ) := by
  classical
  have hp : ∀ v, incidence G x v ≤
      (if G.Adj x v ∧ v ∉ R then (1 : ℤ) else 0) + (if v ∈ R.erase x then 1 else 0) := by
    intro v
    by_cases hv : G.Adj x v
    · have hne := hv.ne.symm
      by_cases hr : v ∈ R <;> simp [incidence, hv, hr, hne]
    · simp [incidence, hv]
      split_ifs <;> omega
  have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hp v)
  simp only [Finset.sum_add_distrib] at hh
  have heq : (((R.erase x).card : ℕ) : ℤ) + 1 = (R.card : ℤ) := by
    exact_mod_cast Finset.card_erase_add_one hx
  have ho : (∑ v, if G.Adj x v ∧ v ∉ R then (1 : ℤ) else 0) =
      ((Finset.univ.filter (fun v => G.Adj x v ∧ v ∉ R)).card : ℤ) := by simp
  have hr : (∑ v, if v ∈ R.erase x then (1 : ℤ) else 0) = ((R.erase x).card : ℤ) := by
    rw [← Finset.sum_filter]
    have herase : Finset.univ.filter (fun v => v ∈ R.erase x) = R.erase x := by
      ext v
      simp
    rw [herase]
    simp
  rw [ho, hr] at hh
  change degreeZ G x ≤ _ at hh
  omega

omit [DecidableEq V] in
theorem degree_lower_bound_from_deficiency (d : ℤ)
    (hdeg : ∀ v, degreeZ G v ≤ d) (x : V) :
    d - deficiency G d ≤ degreeZ G x := by
  have hh := Finset.single_le_sum (f := fun v => d - degreeZ G v)
    (fun v (_ : v ∈ Finset.univ) => sub_nonneg.mpr (hdeg v)) (Finset.mem_univ x)
  change d - degreeZ G x ≤ deficiency G d at hh
  omega

/-- The partition extends whenever the numerical bounds force each exceptional
vertex to have more external neighbors than any one bag can hold. -/
theorem saturated_cycle_coloring
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d) (hsat : ∀ j, degreeZ G (φ j) = d)
    (hbudget : deficiency G d < 8 * d - 3 * (Fintype.card V : ℤ) + 1) :
    ∃ c : V → CycleVertex, (∀ j, c (φ j) = j) ∧
      (∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w)) := by
  have hb := saturated_cycle_partition_bounds G htriangle φ hφ d hsat
  apply cycle_partition_extension G htriangle hnoH φ hφ
  intro x hx
  apply cycle_support_large_of_neighbor_bound G htriangle φ hφ x
  intro j
  have hbag := hb.2 j
  have hlow := degree_lower_bound_from_deficiency G d hdeg x
  have hout := neighbors_outside_set_lower_bound G (cycleRemainder G φ) x hx
  have hR := hb.1
  have hh : ((cycleBag G φ j).card : ℤ) <
      ((Finset.univ.filter (fun v => G.Adj x v ∧ v ∉ cycleRemainder G φ)).card : ℤ) := by omega
  exact_mod_cast hh

end TriangleProof
