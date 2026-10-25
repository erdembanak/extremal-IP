import TriangleProof.GraphCounting

/-! The Möbius ladder exclusion in arbitrary finite triangle-free graphs.
No structural or counting conclusion from the draft is a hypothesis here. -/

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

theorem mobius_neighborhood_inclusion : ∀ i j : Vertex,
    (∀ k, mobius.Adj i k → mobius.Adj j k) → i = j := by decide

theorem mobius_nonedge_common_neighbor : ∀ i j : Vertex,
    ¬ mobius.Adj i j → ∃ k, mobius.Adj i k ∧ mobius.Adj j k := by decide

theorem mobius_triple_incidence : ∀ i j : Vertex,
    incidence mobius (j - 1) i + incidence mobius j i + incidence mobius (j + 1) i
      = 1 + (if i = j then (1 : ℤ) else 0) := by decide

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- Equality in the eight-vertex incidence bound is impossible when each
induced bag has an outside representative and the total deficiency is positive. -/
theorem mobius_copy_rigidity_general
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, mobius.Adj i j → G.Adj (φ i) (φ j))
    (d : ℤ) (hdeg : ∀ v, degreeZ G v ≤ d)
    (hq : 2 ≤ 2 * (Fintype.card V : ℤ) - 5 * d)
    (hpos : 0 < 8 * d - 3 * (Fintype.card V : ℤ))
    (hD : deficiency G d ≤ 8 * d - 3 * (Fintype.card V : ℤ)) : False := by
  classical
  let δ : V → ℤ := fun v => d - degreeZ G v
  let T : V → ℤ := fun v => ∑ i, incidence G v (φ i)
  have hδ : ∀ v, 0 ≤ δ v := fun v => sub_nonneg.mpr (hdeg v)
  have hT : ∀ v, T v ≤ 3 := by
    intro v
    have h := mobius_copy_neighbor_bound G htriangle φ hedges v
    have heq : T v = ((Finset.univ.filter (fun i => G.Adj v (φ i))).card : ℤ) := by
      simp [T, incidence]
    rw [heq]
    exact_mod_cast h
  have hcount : (∑ i, δ (φ i)) = 8 * d - ∑ v, T v := by
    simp only [δ, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul, sum_degreeZ_on_copy, T]
    norm_num
  have hsub : (∑ i, δ (φ i)) ≤ ∑ v, δ v := sum_copy_le_sum φ hφ δ hδ
  have hbudget : (∑ v, δ v) ≤ 8 * d - 3 * (Fintype.card V : ℤ) := hD
  have htotal : (∑ v, T v) ≤ 3 * (Fintype.card V : ℤ) := by
    simpa [mul_comm] using Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hT v)
  have hS : (∑ i, δ (φ i)) = 8 * d - 3 * (Fintype.card V : ℤ) := by omega
  have hglobal : (∑ v, δ v) = 8 * d - 3 * (Fintype.card V : ℤ) := by omega
  have hTeq : ∀ v, T v = 3 := by
    have hz : (∑ v, (3 - T v)) = 0 := by
      rw [Finset.sum_sub_distrib]
      simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
      omega
    intro v
    have := zero_of_sum_eq_zero (fun v => 3 - T v) (fun v => sub_nonneg.mpr (hT v)) hz v
    omega
  have hout : ∀ v, v ∉ Set.range φ → δ v = 0 := by
    have hh := sum_copy_add_outside φ hφ δ
    rw [hS, hglobal] at hh
    have hz : (∑ v ∈ Finset.univ.filter (fun v => v ∉ Set.range φ), δ v) = 0 := by omega
    intro v hv
    exact (Finset.sum_eq_zero_iff_of_nonneg (fun w _ => hδ w)).mp hz v
      (Finset.mem_filter.mpr ⟨Finset.mem_univ _, hv⟩)
  have hpatterns : ∀ v, ∃ j : Vertex, ∀ i, G.Adj v (φ i) ↔ mobius.Adj j i := by
    intro v
    let s := Finset.univ.filter (fun i => G.Adj v (φ i))
    have hind : Independent mobius s := by
      intro i hi j hj hij
      exact htriangle v (φ i) (φ j)
        ⟨(Finset.mem_filter.mp hi).2, hedges i j hij, ((Finset.mem_filter.mp hj).2).symm⟩
    have hcard : s.card = 3 := by
      have heq : (s.card : ℤ) = T v := by simp [s, T, incidence]
      have := hTeq v
      omega
    obtain ⟨j, hj⟩ := (mobius_independent_triples s).mp ⟨hind, hcard⟩
    refine ⟨j, fun i => ?_⟩
    have := Finset.ext_iff.mp hj i
    simpa [s] using this
  choose p hp using hpatterns
  have hpφ : ∀ j, p (φ j) = j := by
    intro j
    exact (mobius_neighborhood_inclusion j (p (φ j))
      (fun i hi => (hp (φ j) i).mp (hedges j i hi))).symm
  have hhom : ∀ v w, G.Adj v w → mobius.Adj (p v) (p w) := by
    intro v w hvw
    by_contra h
    obtain ⟨k, hvk, hwk⟩ := mobius_nonedge_common_neighbor (p v) (p w) h
    exact htriangle v w (φ k) ⟨hvw, (hp w k).mpr hwk, ((hp v k).mpr hvk).symm⟩
  have hdom : ∀ v, degreeZ G v ≤ degreeZ G (φ (p v)) := by
    intro v
    apply Finset.sum_le_sum
    intro w _
    have hh : G.Adj v w → G.Adj (φ (p v)) w := by
      intro h
      exact ((hp w (p v)).mpr (hhom v w h).symm).symm
    dsimp [incidence]
    split_ifs <;> first | omega | tauto
  have hinc : ∀ j v, incidence G (φ j) v = incidence mobius j (p v) := by
    intro j v
    have hh : G.Adj (φ j) v ↔ mobius.Adj j (p v) := by
      rw [G.adj_comm, hp, mobius.adj_comm]
    simp only [incidence, hh]
  have hbag : ∀ j, degreeZ G (φ (j - 1)) + degreeZ G (φ j) + degreeZ G (φ (j + 1))
      = (Fintype.card V : ℤ) + ((Finset.univ.filter (fun v => p v = j)).card : ℤ) := by
    intro j
    simp only [degreeZ, hinc, ← Finset.sum_add_distrib, mobius_triple_incidence]
    simp [Finset.sum_add_distrib]
  have hthree : ∀ j, δ (φ (j - 1)) + δ (φ j) + δ (φ (j + 1)) ≤ ∑ v, δ v := by
    intro j
    let ψ : Fin 3 → V := fun k => φ (if k = 0 then j - 1 else if k = 1 then j else j + 1)
    have hψ : Function.Injective ψ := by
      have hh : ∀ j : Vertex, Function.Injective
          (fun k : Fin 3 => if k = 0 then j - 1 else if k = 1 then j else j + 1) := by decide
      exact hφ.comp (hh j)
    have hh := sum_copy_le_sum ψ hψ δ hδ
    simpa [ψ, Fin.sum_univ_succ, add_assoc] using hh
  have hlarge : ∀ j, 2 ≤ (Finset.univ.filter (fun v => p v = j)).card := by
    intro j
    have hb := hbag j
    have hh := hthree j
    rw [hglobal] at hh
    dsimp [δ] at hh
    have : (2 : ℤ) ≤ ((Finset.univ.filter (fun v => p v = j)).card : ℤ) := by omega
    exact_mod_cast this
  have hsat : ∀ j, δ (φ j) = 0 := by
    intro j
    have hex : ∃ v, p v = j ∧ v ∉ Set.range φ := by
      by_contra h
      push Not at h
      have hsubset : Finset.univ.filter (fun v => p v = j) ⊆ {φ j} := by
        intro v hv
        have hpv := (Finset.mem_filter.mp hv).2
        obtain ⟨k, rfl⟩ := h v hpv
        rw [hpφ] at hpv
        simp [hpv]
      have hh := Finset.card_le_card hsubset
      have := hlarge j
      simp only [Finset.card_singleton] at hh
      omega
    obtain ⟨v, hpv, hv⟩ := hex
    have hvδ := hout v hv
    have hh := hdom v
    rw [hpv] at hh
    dsimp [δ] at hvδ ⊢
    have := hdeg (φ j)
    omega
  have hz : (∑ i, δ (φ i)) = 0 := by simp [hsat]
  omega

/-- The original fixed-order specialization of the general rigidity lemma. -/
theorem mobius_copy_rigidity
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V) (hφ : Function.Injective φ)
    (hedges : ∀ i j, mobius.Adj i j → G.Adj (φ i) (φ j))
    (t : ℤ) (ht : 3 ≤ t) (hn : (Fintype.card V : ℤ) = 10 * t - 1)
    (hdeg : ∀ v, degreeZ G v ≤ 4 * t - 1)
    (hD : deficiency G (4 * t - 1) ≤ 2 * t - 5) : False := by
  exact mobius_copy_rigidity_general G htriangle φ hφ hedges (4 * t - 1)
    hdeg (by omega) (by omega) (by omega)

end TriangleProof
