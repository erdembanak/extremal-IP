import TriangleProof.CycleCapacity

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

def cycleIndependent (s : Finset CycleVertex) : Prop :=
  ∀ i ∈ s, ∀ j ∈ s, ¬ fiveCycle.Adj i j

instance (s : Finset CycleVertex) : Decidable (cycleIndependent s) :=
  inferInstanceAs (Decidable (∀ i ∈ s, ∀ j ∈ s, ¬ fiveCycle.Adj i j))

theorem cycle_independent_bound : ∀ s, cycleIndependent s → s.card ≤ 2 := by decide

theorem cycle_neighborhood_card : ∀ j : CycleVertex, (fiveCycle.neighborFinset j).card = 2 := by decide

theorem cycle_neighborhood_injective : Function.Injective (fun j => fiveCycle.neighborFinset j) := by decide

theorem cycle_pattern_of_containing_neighborhood : ∀ s j,
    cycleIndependent s → (∀ k, fiveCycle.Adj j k → k ∈ s) →
    s = fiveCycle.neighborFinset j := by decide

theorem cycle_bag_cover : ∀ (s : Finset CycleVertex) (j : CycleVertex),
    cycleIndependent s →
    (1 : ℤ) ≤ (if s = fiveCycle.neighborFinset j then 1 else 0) +
      (if s.card < 2 then 1 else 0) +
      (if j + 2 ∈ s then 1 else 0) + (if j + 3 ∈ s then 1 else 0) := by decide

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

def IsCycleCopy (φ : CycleVertex → V) : Prop :=
  Function.Injective φ ∧ ∀ i j, fiveCycle.Adj i j → G.Adj (φ i) (φ j)

instance (φ : CycleVertex → V) : Decidable (IsCycleCopy G φ) :=
  inferInstanceAs (Decidable (Function.Injective φ ∧ ∀ i j, fiveCycle.Adj i j → G.Adj (φ i) (φ j)))

/-- Starting from any five-cycle, the counterexample deficiency bound forces
a five-cycle with every vertex saturated. `InitialCycle.lean` supplies the
initial cycle used by the final fixed-order theorem. -/
theorem saturated_cycle_of_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (d t : ℤ) (hn : (Fintype.card V : ℤ) = 2 * d + 2 * t + 1)
    (hdeg : ∀ v, degreeZ G v ≤ d)
    (hD : deficiency G d ≤ d - 2 * t - 4)
    (hcycle : ∃ φ, IsCycleCopy G φ) :
    ∃ φ, IsCycleCopy G φ ∧ ∀ j, degreeZ G (φ j) = d := by
  classical
  let cycles := Finset.univ.filter (IsCycleCopy G)
  have hnonempty : cycles.Nonempty := by
    obtain ⟨φ, hφ⟩ := hcycle
    exact ⟨φ, by simp [cycles, hφ]⟩
  obtain ⟨φ, hφmem, hmax⟩ := cycles.exists_max_image (fun φ => ∑ j, degreeZ G (φ j)) hnonempty
  have hφ : IsCycleCopy G φ := (Finset.mem_filter.mp hφmem).2
  let N := fun v => Finset.univ.filter (fun i => G.Adj v (φ i))
  have hN : ∀ v, cycleIndependent (N v) := by
    intro v i hi j hj hij
    exact htriangle v (φ i) (φ j)
      ⟨(Finset.mem_filter.mp hi).2, hφ.2 i j hij, ((Finset.mem_filter.mp hj).2).symm⟩
  have hNφ : ∀ j, N (φ j) = fiveCycle.neighborFinset j := by
    intro j
    exact cycle_pattern_of_containing_neighborhood _ j (hN _) (by
      intro k hk
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hφ.2 j k hk⟩)
  let R := Finset.univ.filter (fun v => (N v).card < 2)
  let B := fun j => Finset.univ.filter (fun v => N v = fiveCycle.neighborFinset j)
  let O := Finset.univ.filter (fun v => v ∉ Set.range φ)
  let δ := fun v => d - degreeZ G v
  let Dout := ∑ v ∈ O, δ v
  let T := fun v => ((N v).card : ℤ)
  let q := ∑ v, (2 - T v)
  have hδ : ∀ v, 0 ≤ δ v := fun v => sub_nonneg.mpr (hdeg v)
  have hTout : (∑ j, degreeZ G (φ j)) = ∑ v, T v := by
    rw [sum_degreeZ_on_copy]
    simp [T, N, incidence]
  have hq : q = 2 * (Fintype.card V : ℤ) - ∑ j, degreeZ G (φ j) := by
    dsimp [q]
    rw [Finset.sum_sub_distrib, hTout]
    simp [mul_comm]
  have hsplit : (∑ j : CycleVertex, δ (φ j)) + Dout = deficiency G d :=
    sum_copy_add_outside φ hφ.1 δ
  have hbudget : q + Dout ≤ 2 * t - 2 := by
    rw [Finset.sum_sub_distrib] at hsplit
    norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsplit
    change 5 * d - (∑ j : CycleVertex, degreeZ G (φ j)) + Dout = _ at hsplit
    omega
  have hrq : (R.card : ℤ) ≤ q := by
    have hh : ∀ v, (if v ∈ R then (1 : ℤ) else 0) ≤ 2 - T v := by
      intro v
      have hcard := cycle_independent_bound (N v) (hN v)
      dsimp [T]
      by_cases hv : v ∈ R
      · have hv' : (N v).card < 2 := (Finset.mem_filter.mp hv).2
        simp only [if_pos hv]
        omega
      · simp only [if_neg hv]
        omega
    have hsum := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hh v)
    simpa [q] using hsum
  have hBlarge : ∀ j, 2 * t + 1 - (R.card : ℤ) ≤ (B j).card := by
    intro j
    have hcover := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => cycle_bag_cover (N v) j (hN v))
    have hinc : ∀ k, (∑ v, if k ∈ N v then (1 : ℤ) else 0) = degreeZ G (φ k) := by
      intro k
      simp only [N, Finset.mem_filter, Finset.mem_univ, true_and]
      simp only [degreeZ, incidence, G.adj_comm]
    simp only [Finset.sum_add_distrib] at hcover
    rw [hinc, hinc] at hcover
    have hBsum : (∑ v, if N v = fiveCycle.neighborFinset j then (1 : ℤ) else 0) = (B j).card := by
      simp [B]
    have hRsum : (∑ v, if (N v).card < 2 then (1 : ℤ) else 0) = R.card := by simp [R]
    rw [hBsum, hRsum] at hcover
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one] at hcover
    have h₁ := hdeg (φ (j + 2))
    have h₂ := hdeg (φ (j + 3))
    omega
  refine ⟨φ, hφ, fun j => ?_⟩
  have hjB : φ j ∈ B j := Finset.mem_filter.mpr ⟨Finset.mem_univ _, hNφ j⟩
  have hBO : (B j).erase (φ j) ⊆ O := by
    intro v hv
    obtain ⟨hne, hB⟩ := Finset.mem_erase.mp hv
    refine Finset.mem_filter.mpr ⟨Finset.mem_univ _, ?_⟩
    rintro ⟨k, rfl⟩
    have hh := (Finset.mem_filter.mp hB).2
    rw [hNφ] at hh
    have hkj := cycle_neighborhood_injective hh
    exact hne (congrArg φ hkj)
  have hsmall : (∑ v ∈ (B j).erase (φ j), δ v) ≤ Dout :=
    Finset.sum_le_sum_of_subset_of_nonneg hBO (fun v _ _ => hδ v)
  have hsize : Dout < (((B j).erase (φ j)).card : ℤ) := by
    have hh := hBlarge j
    have hc := Finset.card_erase_add_one hjB
    have hc' : (((B j).erase (φ j)).card : ℤ) + 1 = ((B j).card : ℤ) := by exact_mod_cast hc
    omega
  have hex : ∃ v ∈ (B j).erase (φ j), δ v = 0 := by
    by_contra h
    push Not at h
    have hh : ∀ v ∈ (B j).erase (φ j), (1 : ℤ) ≤ δ v := by
      intro v hv
      have := h v hv
      have := hδ v
      omega
    have hsum := Finset.sum_le_sum hh
    simp only [Finset.sum_const, nsmul_eq_mul, mul_one] at hsum
    omega
  obtain ⟨v, hv, hvδ⟩ := hex
  have hvout : v ∉ Set.range φ := (Finset.mem_filter.mp (hBO hv)).2
  have hvB : N v = fiveCycle.neighborFinset j := (Finset.mem_filter.mp (Finset.mem_erase.mp hv).2).2
  have hpattern : ∀ k, G.Adj v (φ k) ↔ fiveCycle.Adj j k := by
    intro k
    have hh := Finset.ext_iff.mp hvB k
    simpa [N] using hh
  let ψ := Function.update φ j v
  have hψ : IsCycleCopy G ψ := by
    constructor
    · intro a b hab
      by_cases ha : a = j <;> by_cases hb : b = j
      · exact ha.trans hb.symm
      · subst a
        have hh : v = φ b := by simpa [ψ, hb] using hab
        exact (hvout ⟨b, hh.symm⟩).elim
      · subst b
        have hh : φ a = v := by simpa [ψ, ha] using hab
        exact (hvout ⟨a, hh⟩).elim
      · apply hφ.1
        simpa [ψ, ha, hb] using hab
    · intro a b hab
      by_cases ha : a = j <;> by_cases hb : b = j
      · subst a; subst b; exact (fiveCycle.loopless.irrefl _ hab).elim
      · subst a; simpa [ψ, hb] using (hpattern b).mpr hab
      · subst b; simpa [ψ, ha] using ((hpattern a).mpr hab.symm).symm
      · simpa [ψ, ha, hb] using hφ.2 a b hab
  have hmax' := hmax ψ (Finset.mem_filter.mpr ⟨Finset.mem_univ _, hψ⟩)
  have hsum : (∑ i, degreeZ G (ψ i)) =
      (∑ i, degreeZ G (φ i)) - degreeZ G (φ j) + degreeZ G v := by
    have h₁ := Finset.sum_erase_add Finset.univ (fun i => degreeZ G (ψ i)) (Finset.mem_univ j)
    have h₂ := Finset.sum_erase_add Finset.univ (fun i => degreeZ G (φ i)) (Finset.mem_univ j)
    have heq : (∑ i ∈ Finset.univ.erase j, degreeZ G (ψ i)) =
        ∑ i ∈ Finset.univ.erase j, degreeZ G (φ i) := by
      apply Finset.sum_congr rfl
      intro i hi
      simp [ψ, (Finset.mem_erase.mp hi).1]
    rw [heq] at h₁
    have hψj : ψ j = v := by simp [ψ]
    rw [hψj] at h₁
    omega
  rw [hsum] at hmax'
  have hbound := hdeg (φ j)
  dsimp [δ] at hvδ
  omega

end TriangleProof
