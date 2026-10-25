import TriangleProof.MobiusExclusion

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

def extraPair (i j a b : Vertex) : Prop := (i = a ∧ j = b) ∨ (j = a ∧ i = b)

instance (i j a b : Vertex) : Decidable (extraPair i j a b) :=
  inferInstanceAs (Decidable ((i = a ∧ j = b) ∨ (j = a ∧ i = b)))

def augmentedRel (bits : Fin 4 → Bool) (i j : Vertex) : Prop :=
  obstruction.Adj i j ∨
  (bits 0 = true ∧ extraPair i j 4 6) ∨
  (bits 1 = true ∧ extraPair i j 4 7) ∨
  (bits 2 = true ∧ extraPair i j 5 6) ∨
  (bits 3 = true ∧ extraPair i j 5 7)

instance (bits : Fin 4 → Bool) (i j : Vertex) : Decidable (augmentedRel bits i j) :=
  inferInstanceAs (Decidable (obstruction.Adj i j ∨
    (bits 0 = true ∧ extraPair i j 4 6) ∨ (bits 1 = true ∧ extraPair i j 4 7) ∨
    (bits 2 = true ∧ extraPair i j 5 6) ∨ (bits 3 = true ∧ extraPair i j 5 7)))

def augmentation (bits : Fin 4 → Bool) : SimpleGraph Vertex where
  Adj := augmentedRel bits
  symm := ⟨(by decide : ∀ bits i j, augmentedRel bits i j → augmentedRel bits j i) bits⟩
  loopless := ⟨(by decide : ∀ bits i, ¬ augmentedRel bits i i) bits⟩

instance (bits : Fin 4 → Bool) : DecidableRel (augmentation bits).Adj :=
  fun i j => inferInstanceAs (Decidable (augmentedRel bits i j))

def mobiusPermutation₁ : Vertex → Vertex := ![0, 4, 5, 1, 2, 6, 7, 3]
def mobiusPermutation₂ : Vertex → Vertex := ![0, 4, 5, 1, 3, 7, 6, 2]

theorem mobiusPermutation₁_injective : Function.Injective mobiusPermutation₁ := by decide
theorem mobiusPermutation₂_injective : Function.Injective mobiusPermutation₂ := by decide

theorem obstruction_pair_classification : ∀ i j : Vertex,
    obstruction.Adj i j ∨ i = j ∨
    (∃ k, obstruction.Adj i k ∧ obstruction.Adj j k) ∨
    extraPair i j 4 6 ∨ extraPair i j 4 7 ∨ extraPair i j 5 6 ∨ extraPair i j 5 7 := by decide

/-- All sixteen choices of the four possible additional edges are checked
by the kernel. Two edges either create a triangle or one of the two explicit
Möbius-ladder copies. -/
theorem augmentation_degree_bound : ∀ bits : Fin 4 → Bool,
    (∀ a b c, ¬ ((augmentation bits).Adj a b ∧ (augmentation bits).Adj b c ∧
      (augmentation bits).Adj c a)) →
    ¬ (∀ i j, mobius.Adj i j → (augmentation bits).Adj (mobiusPermutation₁ i) (mobiusPermutation₁ j)) →
    ¬ (∀ i j, mobius.Adj i j → (augmentation bits).Adj (mobiusPermutation₂ i) (mobiusPermutation₂ j)) →
    (∑ i, degreeZ (augmentation bits) i) ≤ 22 := by decide

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

omit [Fintype V] [DecidableEq V] in
theorem obstruction_copy_induced_degree_bound
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoM : ∀ (ψ : Vertex → V), Function.Injective ψ →
      ¬ (∀ i j, mobius.Adj i j → G.Adj (ψ i) (ψ j)))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, obstruction.Adj i j → G.Adj (φ i) (φ j)) :
    (∑ i : Vertex, ∑ j : Vertex, incidence G (φ i) (φ j)) ≤ 22 := by
  classical
  let bits : Fin 4 → Bool := ![decide (G.Adj (φ 4) (φ 6)), decide (G.Adj (φ 4) (φ 7)),
    decide (G.Adj (φ 5) (φ 6)), decide (G.Adj (φ 5) (φ 7))]
  have heq : ∀ i j, G.Adj (φ i) (φ j) ↔ (augmentation bits).Adj i j := by
    intro i j
    change _ ↔ augmentedRel bits i j
    constructor
    · intro hij
      rcases obstruction_pair_classification i j with h | h | ⟨k, hik, hjk⟩ | h | h | h | h
      · exact Or.inl h
      · subst j; exact (G.loopless.irrefl _ hij).elim
      · exact (htriangle (φ i) (φ j) (φ k) ⟨hij, hedges j k hjk, (hedges i k hik).symm⟩).elim
      all_goals
        rcases h with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
          simp [augmentedRel, bits, extraPair, hij, hij.symm]
    · intro h
      rcases h with h | ⟨hb, hp⟩ | ⟨hb, hp⟩ | ⟨hb, hp⟩ | ⟨hb, hp⟩
      · exact hedges i j h
      all_goals
        rcases hp with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
          simpa [bits, G.adj_comm] using hb
  have htri : ∀ a b c, ¬ ((augmentation bits).Adj a b ∧ (augmentation bits).Adj b c ∧
      (augmentation bits).Adj c a) := by
    intro a b c ⟨hab, hbc, hca⟩
    exact htriangle (φ a) (φ b) (φ c) ⟨(heq a b).mpr hab, (heq b c).mpr hbc, (heq c a).mpr hca⟩
  have hn₁ : ¬ (∀ i j, mobius.Adj i j →
      (augmentation bits).Adj (mobiusPermutation₁ i) (mobiusPermutation₁ j)) := by
    intro h
    exact hnoM (φ ∘ mobiusPermutation₁) (hφ.comp mobiusPermutation₁_injective)
      (fun i j hij => (heq _ _).mpr (h i j hij))
  have hn₂ : ¬ (∀ i j, mobius.Adj i j →
      (augmentation bits).Adj (mobiusPermutation₂ i) (mobiusPermutation₂ j)) := by
    intro h
    exact hnoM (φ ∘ mobiusPermutation₂) (hφ.comp mobiusPermutation₂_injective)
      (fun i j hij => (heq _ _).mpr (h i j hij))
  have hh := augmentation_degree_bound bits htri hn₁ hn₂
  simpa only [degreeZ, incidence, heq] using hh

theorem mobius_copy_deficiency_lower_bound
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, mobius.Adj i j → G.Adj (φ i) (φ j))
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d) :
    8 * d - 3 * (Fintype.card V : ℤ) ≤ deficiency G d := by
  classical
  let T : V → ℤ := fun v => ∑ i : Vertex, incidence G v (φ i)
  have hT : ∀ v, T v ≤ 3 := by
    intro v
    have h := mobius_copy_neighbor_bound G htriangle φ hedges v
    have heq : T v = ((Finset.univ.filter (fun i => G.Adj v (φ i))).card : ℤ) := by
      simp [T, incidence]
    rw [heq]
    exact_mod_cast h
  have htotal : (∑ v, T v) ≤ 3 * (Fintype.card V : ℤ) := by
    simpa [mul_comm] using Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hT v)
  have hcopy : (∑ i, (d - degreeZ G (φ i))) = 8 * d - ∑ v, T v := by
    simp only [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul, sum_degreeZ_on_copy, T]
    norm_num
  have hsub := sum_copy_le_sum φ hφ (fun v => d - degreeZ G v)
    (fun v => sub_nonneg.mpr (hdeg v))
  change _ ≤ deficiency G d at hsub
  rw [hcopy] at hsub
  omega

theorem obstruction_copy_deficiency_lower_bound
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoM : ∀ (ψ : Vertex → V), Function.Injective ψ →
      ¬ (∀ i j, mobius.Adj i j → G.Adj (ψ i) (ψ j)))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, obstruction.Adj i j → G.Adj (φ i) (φ j))
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d) :
    8 * d - 3 * (Fintype.card V : ℤ) + 2 ≤ deficiency G d := by
  classical
  have hin := obstruction_copy_induced_degree_bound G htriangle hnoM φ hφ hedges
  let T : V → ℤ := fun v => ∑ i : Vertex, incidence G v (φ i)
  have hT : ∀ v, T v ≤ 3 := by
    intro v
    have h := obstruction_copy_neighbor_bound G htriangle φ hedges v
    have heq : T v = ((Finset.univ.filter (fun i => G.Adj v (φ i))).card : ℤ) := by
      simp [T, incidence]
    rw [heq]
    exact_mod_cast h
  have hgap : 2 ≤ ∑ v, (3 - T v) := by
    have hsub := sum_copy_le_sum φ hφ (fun v => 3 - T v)
      (fun v => sub_nonneg.mpr (hT v))
    have hin' : (∑ i, T (φ i)) ≤ 22 := hin
    rw [Finset.sum_sub_distrib] at hsub
    norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul] at hsub
    change 24 - (∑ i : Vertex, T (φ i)) ≤ _ at hsub
    omega
  have hcount : (∑ i, degreeZ G (φ i)) ≤ 3 * (Fintype.card V : ℤ) - 2 := by
    rw [Finset.sum_sub_distrib] at hgap
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hgap
    rw [sum_degreeZ_on_copy]
    change (∑ v, T v) ≤ _
    omega
  have hsub := sum_copy_le_sum φ hφ (fun v => d - degreeZ G v)
    (fun v => sub_nonneg.mpr (hdeg v))
  change _ ≤ deficiency G d at hsub
  rw [Finset.sum_sub_distrib] at hsub
  norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul] at hsub
  change 8 * d - (∑ i : Vertex, degreeZ G (φ i)) ≤ _ at hsub
  omega

/-- The ten-edge obstruction cannot occur under the counterexample bounds.
This includes its induced-edge bound and the count of edges leaving the copy. -/
theorem obstruction_copy_exclusion
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, obstruction.Adj i j → G.Adj (φ i) (φ j))
    (t : ℤ) (ht : 3 ≤ t) (hn : (Fintype.card V : ℤ) = 10 * t - 1)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t - 1)
    (hD : deficiency G (4 * t - 1) ≤ 2 * t - 5) : False := by
  classical
  have hnoM := fun ψ hψ hE => mobius_copy_rigidity G htriangle ψ hψ hE t ht hn hdeg hD
  have hin := obstruction_copy_induced_degree_bound G htriangle hnoM φ hφ hedges
  let T : V → ℤ := fun v => ∑ i : Vertex, incidence G v (φ i)
  have hT : ∀ v, T v ≤ 3 := by
    intro v
    have h := obstruction_copy_neighbor_bound G htriangle φ hedges v
    have heq : T v = ((Finset.univ.filter (fun i => G.Adj v (φ i))).card : ℤ) := by
      simp [T, incidence]
    rw [heq]
    exact_mod_cast h
  have hgap : 2 ≤ ∑ v, (3 - T v) := by
    have hsub := sum_copy_le_sum φ hφ (fun v => 3 - T v) (fun v => sub_nonneg.mpr (hT v))
    have hin' : (∑ i, T (φ i)) ≤ 22 := hin
    rw [Finset.sum_sub_distrib] at hsub
    norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsub
    change 24 - (∑ i : Vertex, T (φ i)) ≤ _ at hsub
    omega
  have hcount : (∑ i, degreeZ G (φ i)) ≤ 3 * (Fintype.card V : ℤ) - 2 := by
    rw [Finset.sum_sub_distrib] at hgap
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hgap
    rw [sum_degreeZ_on_copy]
    change (∑ v, T v) ≤ _
    omega
  have hsub := sum_copy_le_sum φ hφ (fun v => 4 * t - 1 - degreeZ G v)
    (fun v => sub_nonneg.mpr (hdeg v))
  change _ ≤ deficiency G (4 * t - 1) at hsub
  rw [Finset.sum_sub_distrib] at hsub
  norm_num only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsub
  change 8 * (4 * t - 1) - (∑ i : Vertex, degreeZ G (φ i)) ≤ _ at hsub
  omega

end TriangleProof
