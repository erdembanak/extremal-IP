import TriangleProof.GraphCounting
import TriangleProof.Capacity

namespace TriangleProof

abbrev CycleVertex := Fin 5

def cycleRel (i j : CycleVertex) : Prop :=
  (i.val + 1) % 5 = j.val ∨ (j.val + 1) % 5 = i.val

instance (i j : CycleVertex) : Decidable (cycleRel i j) :=
  inferInstanceAs (Decidable ((i.val + 1) % 5 = j.val ∨ (j.val + 1) % 5 = i.val))

def fiveCycle : SimpleGraph CycleVertex where
  Adj := cycleRel
  symm := ⟨by decide⟩
  loopless := ⟨by decide⟩

instance : DecidableRel fiveCycle.Adj :=
  fun i j => inferInstanceAs (Decidable (cycleRel i j))

def cyclePotential (s : CycleVertex → ℤ) (j : CycleVertex) := s (j - 1) + s (j + 1)

theorem cyclePotential_eq_sum (s : CycleVertex → ℤ) (j : CycleVertex) :
    (∑ i, s i * incidence fiveCycle j i) = cyclePotential s j := by
  have hneg : (-1 : CycleVertex) = 4 := by decide
  fin_cases j <;>
    simp [Fin.sum_univ_succ, incidence, fiveCycle, cycleRel, cyclePotential, hneg] <;>
    ring

theorem cyclePotential_sum (s : CycleVertex → ℤ) :
    (∑ j, cyclePotential s j) = 2 * ∑ j, s j := by
  simp only [cyclePotential, Fin.sum_univ_succ]
  dsimp [Fin.add, Fin.sub, Fin.succ]
  ring

theorem cyclePotential_bag_identity (s : CycleVertex → ℤ) (j : CycleVertex) :
    s j + cyclePotential s (j + 2) + cyclePotential s (j + 3) = ∑ i, s i := by
  fin_cases j <;> simp only [cyclePotential, Fin.sum_univ_succ] <;>
    dsimp [Fin.add, Fin.sub, Fin.succ] <;> ring

/-- Positions carrying positive excess weight. -/
def cyclePositiveSupport (p : CycleVertex → ℕ) : Finset CycleVertex :=
  Finset.univ.filter fun i => p i ≠ 0

/-- The positive support can be split into two independent sets of the
five-cycle. This is the exact combinatorial premise of the odd-excess
imbalance argument in the manuscript. -/
def CycleSupportBipartite (p : CycleVertex → ℕ) : Prop :=
  ∃ I J : Finset CycleVertex,
    I ∪ J = cyclePositiveSupport p ∧ Disjoint I J ∧
    (∀ i ∈ I, ∀ j ∈ I, ¬ fiveCycle.Adj i j) ∧
    (∀ i ∈ J, ∀ j ∈ J, ¬ fiveCycle.Adj i j)

/-- The algebraic expansion behind the odd-excess lemma. For an independent
side of the support, cross terms cancel, terms on the other side are
nonnegative, and the same-side cost is at most `P(P-1)`. -/
theorem cycle_weight_bipartition_lower_bound
    (p : CycleVertex → ℕ) (h : ℤ) (I : Finset CycleVertex)
    (hI : ∀ i ∈ I, ∀ j ∈ I, ¬ fiveCycle.Adj i j) :
    let P := ∑ i ∈ I, (p i : ℤ)
    let Q := ∑ i, (p i : ℤ)
    let w := fun j => (p j : ℤ) *
      (h - (p (j + 2) : ℤ) - (p (j + 3) : ℤ))
    h * (2 * P - Q) - P * (P - 1) ≤
      2 * (∑ i ∈ I, w i) - ∑ i, w i := by
  classical
  dsimp only
  have hprod : ∀ i j : CycleVertex, 0 ≤ (p i : ℤ) * (p j : ℤ) := by
    intro i j
    positivity
  have hsq : ∀ i : CycleVertex, (p i : ℤ) ≤ (p i : ℤ) * (p i : ℤ) := by
    intro i
    rcases Nat.eq_zero_or_pos (p i) with hi | hi
    · simp [hi]
    · have hiZ : (1 : ℤ) ≤ (p i : ℤ) := by exact_mod_cast hi
      nlinarith [mul_nonneg (show 0 ≤ (p i : ℤ) by positivity)
        (show 0 ≤ (p i : ℤ) - 1 by omega)]
  fin_cases I <;>
    simp [Fin.sum_univ_succ, fiveCycle, cycleRel] at hI ⊢ <;>
    nlinarith [hprod 0 2, hprod 0 3, hprod 1 3, hprod 1 4, hprod 2 4,
      hsq 0, hsq 1, hsq 2, hsq 3, hsq 4]

/-- Every induced subgraph of the five-cycle on at most three vertices is
bipartite. The finite check concerns only the 32 subsets of `Fin 5`, not the
unbounded excess weights. -/
theorem cycle_support_bipartite_of_card_le_three
    (S : Finset CycleVertex) (hS : S.card ≤ 3) :
    ∃ I J : Finset CycleVertex,
      I ∪ J = S ∧ Disjoint I J ∧
      (∀ i ∈ I, ∀ j ∈ I, ¬ fiveCycle.Adj i j) ∧
      (∀ i ∈ J, ∀ j ∈ J, ¬ fiveCycle.Adj i j) := by
  fin_cases S <;>
    solve
    | exfalso; norm_num at hS
    | exact of_decide_eq_true rfl

theorem cyclePositiveSupport_card_le_sum (p : CycleVertex → ℕ) :
    (cyclePositiveSupport p).card ≤ ∑ i, p i := by
  calc
    (cyclePositiveSupport p).card =
        ∑ i ∈ cyclePositiveSupport p, 1 := by simp
    _ ≤ ∑ i ∈ cyclePositiveSupport p, p i := by
      apply Finset.sum_le_sum
      intro i hi
      exact Nat.one_le_iff_ne_zero.mpr (by
        simpa [cyclePositiveSupport] using hi)
    _ ≤ ∑ i, p i := by
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) (by
        simp)

theorem cycle_support_bipartite_of_sum_le_three
    (p : CycleVertex → ℕ) (hsum : (∑ i, p i) ≤ 3) :
    CycleSupportBipartite p := by
  exact cycle_support_bipartite_of_card_le_three (cyclePositiveSupport p)
    (le_trans (cyclePositiveSupport_card_le_sum p) hsum)

/-- General odd-excess imbalance for a bipartite positive support. This is
the Lean counterpart of the manuscript's weighted bipartition proof. -/
theorem odd_excess_weight_imbalance
    (p : CycleVertex → ℕ) (h k : ℤ)
    (hk : 1 ≤ k)
    (hsum : (∑ i, (p i : ℤ)) = 2 * k - 1)
    (hh : 2 * k - 1 ≤ h)
    (hbip : CycleSupportBipartite p) :
    ∃ I : Finset CycleVertex,
      (∀ i ∈ I, ∀ j ∈ I, ¬ fiveCycle.Adj i j) ∧
      h - k * (k - 1) ≤
        2 * (∑ i ∈ I, (p i : ℤ) *
          (h - (p (i + 2) : ℤ) - (p (i + 3) : ℤ))) -
        ∑ i, (p i : ℤ) *
          (h - (p (i + 2) : ℤ) - (p (i + 3) : ℤ)) := by
  classical
  obtain ⟨I, J, hcover, hdisjoint, hI, hJ⟩ := hbip
  let P := ∑ i ∈ I, (p i : ℤ)
  let Q := ∑ i ∈ J, (p i : ℤ)
  have hsupp : (∑ i ∈ cyclePositiveSupport p, (p i : ℤ)) =
      ∑ i, (p i : ℤ) := by
    apply Finset.sum_subset (Finset.subset_univ _)
    intro i _ hi
    simp only [cyclePositiveSupport, Finset.mem_filter, Finset.mem_univ,
      true_and, not_ne_iff] at hi
    simp [hi]
  have htotal : (∑ i, (p i : ℤ)) = P + Q := by
    rw [← hsupp, ← hcover, Finset.sum_union hdisjoint]
  have hP0 : 0 ≤ P := by
    dsimp [P]
    positivity
  have hQ0 : 0 ≤ Q := by
    dsimp [Q]
    positivity
  by_cases hQI : Q ≤ P
  · refine ⟨I, hI, ?_⟩
    have hkP : k ≤ P := by omega
    have hPsum : P ≤ 2 * k - 1 := by omega
    have hscalar := odd_excess_imbalance_scalar h k P hk hkP hPsum hh
    have halgebra := cycle_weight_bipartition_lower_bound p h I hI
    exact le_trans hscalar (by simpa [P, hsum] using halgebra)
  · have hPQ : P ≤ Q := by omega
    refine ⟨J, hJ, ?_⟩
    have hkQ : k ≤ Q := by omega
    have hQsum : Q ≤ 2 * k - 1 := by omega
    have hscalar := odd_excess_imbalance_scalar h k Q hk hkQ hQsum hh
    have halgebra := cycle_weight_bipartition_lower_bound p h J hJ
    exact le_trans hscalar (by simpa [Q, hsum] using halgebra)

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

omit [DecidableEq V] in
theorem coloring_capacity_inequality (c : V → CycleVertex)
    (hhom : ∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w))
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d)
    (I : Finset CycleVertex)
    (hI : ∀ i ∈ I, ∀ j ∈ I, ¬ fiveCycle.Adj i j) :
    let s := fun j => ((Finset.univ.filter (fun v => c v = j)).card : ℤ)
    let w := fun j => s j * (cyclePotential s j - d)
    2 * (∑ j ∈ I, w j) - (∑ j, w j) ≤ deficiency G d := by
  classical
  dsimp only
  let s := fun j => ((Finset.univ.filter (fun v => c v = j)).card : ℤ)
  let p := fun j => cyclePotential s j - d
  let M := fun v u => incidence fiveCycle (c v) (c u) - incidence G v u
  have hMnonneg : ∀ v u, 0 ≤ M v u := by
    intro v u
    have hh := hhom v u
    dsimp [M, incidence]
    split_ifs <;> first | omega | tauto
  have hMsymm : ∀ v u, M v u = M u v := by
    intro v u
    dsimp [M]
    rw [incidence_symm fiveCycle, incidence_symm G]
  have hMind : ∀ v u, c v ∈ I → c u ∈ I → M v u = 0 := by
    intro v u hv hu
    have hn := hI (c v) hv (c u) hu
    have hg : ¬ G.Adj v u := fun h => hn (hhom v u h)
    simp [M, incidence, hn, hg]
  have hpot : ∀ j, (∑ v, incidence fiveCycle j (c v)) = cyclePotential s j := by
    intro j
    rw [fiber_sum c]
    exact cyclePotential_eq_sum s j
  have hrow : ∀ v, (∑ u, M v u) = p (c v) + (d - degreeZ G v) := by
    intro v
    simp only [M, Finset.sum_sub_distrib, hpot]
    dsimp [p, degreeZ]
    ring
  have htotal : (∑ v, ∑ u, M v u) = (∑ j, s j * p j) + deficiency G d := by
    simp only [hrow, Finset.sum_add_distrib, deficiency]
    rw [fiber_sum c p]
  have hlow : (∑ j ∈ I, s j * p j) ≤ ∑ v, if c v ∈ I then ∑ u, M v u else 0 := by
    have hf : (∑ j ∈ I, s j * p j) = ∑ v, if c v ∈ I then p (c v) else 0 := by
      rw [fiber_sum c (fun j => if j ∈ I then p j else 0)]
      simp only [mul_ite, mul_zero, ← Finset.sum_filter]
      simp [s]
    rw [hf]
    apply Finset.sum_le_sum
    intro v _
    rw [hrow]
    split_ifs <;> have := hdeg v <;> omega
  have hh := independent_weight_bound M hMsymm hMnonneg (fun v => c v ∈ I) hMind
  rw [htotal] at hh
  change 2 * (∑ j ∈ I, s j * p j) - (∑ j, s j * p j) ≤ _
  omega

omit [DecidableEq V] in
/-- A graph mapping to the five-cycle, with a saturated vertex in every
color class, has the deficiency required by Section 6.4 for every t ≥ 3. -/
theorem five_cycle_coloring_deficiency (c : V → CycleVertex)
    (hhom : ∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w))
    (t : ℤ) (ht : 3 ≤ t) (hn : (Fintype.card V : ℤ) = 10 * t - 1)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t - 1)
    (hsat : ∀ j, ∃ v, c v = j ∧ degreeZ G v = 4 * t - 1) :
    2 * t - 1 ≤ deficiency G (4 * t - 1) := by
  classical
  let s := fun j => ((Finset.univ.filter (fun v => c v = j)).card : ℤ)
  let d := 4 * t - 1
  have hsum : (∑ j, s j) = 10 * t - 1 := by
    have hh := fiber_sum c (fun _ => (1 : ℤ))
    simpa [s, hn] using hh.symm
  have hpot : ∀ v, degreeZ G v ≤ cyclePotential s (c v) := by
    intro v
    have hh : degreeZ G v ≤ ∑ u, incidence fiveCycle (c v) (c u) := by
      apply Finset.sum_le_sum
      intro u _
      have hh := hhom v u
      dsimp [incidence]
      split_ifs <;> first | omega | tauto
    rw [fiber_sum c (incidence fiveCycle (c v)), cyclePotential_eq_sum s] at hh
    exact hh
  have hnonneg : ∀ j, 0 ≤ cyclePotential s j - d := by
    intro j
    obtain ⟨v, hv, hdv⟩ := hsat j
    have hh := hpot v
    rw [hv, hdv] at hh
    exact sub_nonneg.mpr hh
  let p : CycleVertex → ℕ := fun j => (cyclePotential s j - d).toNat
  have hp : ∀ j, (p j : ℤ) = cyclePotential s j - d := by
    intro j
    exact Int.toNat_of_nonneg (hnonneg j)
  have hpsum : (∑ j, p j) = 3 := by
    have hh : (∑ j, (p j : ℤ)) = 3 := by
      simp only [hp, Finset.sum_sub_distrib, cyclePotential_sum, hsum,
        Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
      dsimp [d]
      ring
    exact_mod_cast hh
  have hs : ∀ j, s j = 2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ) := by
    intro j
    have hh := cyclePotential_bag_identity s j
    rw [hsum] at hh
    rw [hp, hp]
    dsimp [d]
    omega
  have hbip : CycleSupportBipartite p :=
    cycle_support_bipartite_of_sum_le_three p (by omega)
  have hpsumZ : (∑ j, (p j : ℤ)) = 3 := by
    exact_mod_cast hpsum
  obtain ⟨I, hI, himbalance⟩ :=
    odd_excess_weight_imbalance p (2 * t + 1) (2 : ℤ) (by omega)
      (by norm_num; exact hpsumZ) (by omega) hbip
  have hcapacity := coloring_capacity_inequality G c hhom d hdeg I hI
  change 2 * (∑ j ∈ I, s j * (cyclePotential s j - d)) -
    (∑ j, s j * (cyclePotential s j - d)) ≤ deficiency G d at hcapacity
  have hcapacity' :
      2 * (∑ j ∈ I, (p j : ℤ) *
        (2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ))) -
      (∑ j, (p j : ℤ) *
        (2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ))) ≤
      deficiency G d := by
    simpa only [← hp, hs, mul_comm] using hcapacity
  nlinarith

omit [DecidableEq V] in
/-- Capacity bound for the `d = 4t+1` class, where the five excess
potentials sum to one. -/
theorem five_cycle_coloring_deficiency_one (c : V → CycleVertex)
    (hhom : ∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w))
    (t : ℤ) (ht : 1 ≤ t) (hn : (Fintype.card V : ℤ) = 10 * t + 3)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t + 1)
    (hsat : ∀ j, ∃ v, c v = j ∧ degreeZ G v = 4 * t + 1) :
    2 * t + 1 ≤ deficiency G (4 * t + 1) := by
  classical
  let s := fun j => ((Finset.univ.filter (fun v => c v = j)).card : ℤ)
  let d := 4 * t + 1
  have hsum : (∑ j, s j) = 10 * t + 3 := by
    have hh := fiber_sum c (fun _ => (1 : ℤ))
    simpa [s, hn] using hh.symm
  have hpot : ∀ v, degreeZ G v ≤ cyclePotential s (c v) := by
    intro v
    have hh : degreeZ G v ≤ ∑ u, incidence fiveCycle (c v) (c u) := by
      apply Finset.sum_le_sum
      intro u _
      have hh := hhom v u
      dsimp [incidence]
      split_ifs <;> first | omega | tauto
    rw [fiber_sum c (incidence fiveCycle (c v)), cyclePotential_eq_sum s] at hh
    exact hh
  have hnonneg : ∀ j, 0 ≤ cyclePotential s j - d := by
    intro j
    obtain ⟨v, hv, hdv⟩ := hsat j
    have hh := hpot v
    rw [hv, hdv] at hh
    exact sub_nonneg.mpr hh
  let p : CycleVertex → ℕ := fun j => (cyclePotential s j - d).toNat
  have hp : ∀ j, (p j : ℤ) = cyclePotential s j - d := by
    intro j
    exact Int.toNat_of_nonneg (hnonneg j)
  have hpsum : (∑ j, p j) = 1 := by
    have hh : (∑ j, (p j : ℤ)) = 1 := by
      simp only [hp, Finset.sum_sub_distrib, cyclePotential_sum, hsum,
        Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
      dsimp [d]
      ring
    exact_mod_cast hh
  have hs : ∀ j, s j = 2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ) := by
    intro j
    have hh := cyclePotential_bag_identity s j
    rw [hsum] at hh
    rw [hp, hp]
    dsimp [d]
    omega
  have hbip : CycleSupportBipartite p :=
    cycle_support_bipartite_of_sum_le_three p (by omega)
  have hpsumZ : (∑ j, (p j : ℤ)) = 1 := by
    exact_mod_cast hpsum
  obtain ⟨I, hI, himbalance⟩ :=
    odd_excess_weight_imbalance p (2 * t + 1) (1 : ℤ) (by omega)
      (by norm_num; exact hpsumZ) (by omega) hbip
  have hcapacity := coloring_capacity_inequality G c hhom d hdeg I hI
  change 2 * (∑ j ∈ I, s j * (cyclePotential s j - d)) -
    (∑ j, s j * (cyclePotential s j - d)) ≤ deficiency G d at hcapacity
  have hcapacity' :
      2 * (∑ j ∈ I, (p j : ℤ) *
        (2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ))) -
      (∑ j, (p j : ℤ) *
        (2 * t + 1 - (p (j + 2) : ℤ) - (p (j + 3) : ℤ))) ≤
      deficiency G d := by
    simpa only [← hp, hs, mul_comm] using hcapacity
  nlinarith

end TriangleProof
