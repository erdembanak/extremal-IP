import TriangleProof.ExtremalFoundations
import TriangleProof.FixedOrder

/-! An explicit B(7,9) threshold graph. Its adjacency rule is the published
construction specialized to (7,9), without any premise involving Z(7).
The displayed partner functions witness a perfect matching after each
vertex deletion; all finite facts are checked by the Lean kernel. -/

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

private def sevenCross (a b : Fin 19) : Prop :=
  a.val < 9 ∧ 9 ≤ b.val ∧ b.val < 18 ∧
  ¬ (6 ≤ a.val ∧ 15 ≤ b.val) ∧
  ¬ (a.val < 6 ∧ b.val < 15 ∧
    (b.val - 9 = a.val ∨ b.val - 9 = (a.val + 1) % 6))

private def sevenH (a : Fin 19) : Prop :=
  (6 ≤ a.val ∧ a.val < 9) ∨ (15 ≤ a.val ∧ a.val < 18)

private def sevenRel (a b : Fin 19) : Prop :=
  sevenCross a b ∨ sevenCross b a ∨
  (a = 18 ∧ sevenH b) ∨ (b = 18 ∧ sevenH a)

private instance (a b : Fin 19) : Decidable (sevenCross a b) :=
  inferInstanceAs (Decidable (_ ∧ _ ∧ _ ∧ _ ∧ _))

private instance (a : Fin 19) : Decidable (sevenH a) :=
  inferInstanceAs (Decidable ((_ ∧ _) ∨ (_ ∧ _)))

private instance (a b : Fin 19) : Decidable (sevenRel a b) :=
  inferInstanceAs (Decidable (_ ∨ _ ∨ (_ ∧ _) ∨ (_ ∧ _)))

def sevenGraph : SimpleGraph (Fin 19) where
  Adj := sevenRel
  symm := ⟨by intro a b; unfold sevenRel; tauto⟩
  loopless := ⟨by decide⟩

instance : DecidableRel sevenGraph.Adj := fun a b =>
  inferInstanceAs (Decidable (sevenRel a b))

theorem sevenGraph_triangle_free : IsTriangleFree sevenGraph := by
  change ∀ a b c : Fin 19, ¬ (sevenRel a b ∧ sevenRel b c ∧ sevenRel c a)
  decide

theorem sevenGraph_degree (v : Fin 19) :
    degreeN sevenGraph v = if v = 18 then 6 else 7 := by
  rw [degreeN_eq_degree]
  exact (by decide : ∀ v : Fin 19, sevenGraph.degree v = if v = 18 then 6 else 7) v

/-- Row r lists the partner of every vertex after r is deleted.
The unused entry at r is fixed to r. -/
def sevenPartner : Fin 19 → Fin 19 → Fin 19 :=
  ![
    ![0, 12, 13, 14, 17, 16, 11, 10, 9, 8, 7, 6, 1, 2, 3, 18, 5, 4, 15],
    ![12, 1, 13, 14, 17, 16, 11, 10, 9, 8, 7, 6, 0, 2, 3, 18, 5, 4, 15],
    ![16, 12, 2, 14, 17, 13, 11, 10, 9, 8, 7, 6, 1, 5, 3, 18, 0, 4, 15],
    ![12, 16, 14, 3, 17, 13, 11, 10, 9, 8, 7, 6, 0, 5, 2, 18, 1, 4, 15],
    ![12, 13, 14, 16, 4, 17, 11, 10, 9, 8, 7, 6, 0, 1, 2, 18, 3, 5, 15],
    ![17, 13, 14, 16, 12, 5, 11, 10, 9, 8, 7, 6, 4, 1, 2, 18, 3, 0, 15],
    ![14, 17, 13, 11, 16, 12, 6, 10, 9, 8, 7, 3, 5, 2, 0, 18, 4, 1, 15],
    ![14, 17, 13, 11, 16, 12, 10, 7, 9, 8, 6, 3, 5, 2, 0, 18, 4, 1, 15],
    ![14, 17, 13, 11, 16, 12, 10, 9, 8, 7, 6, 3, 5, 2, 0, 18, 4, 1, 15],
    ![17, 16, 13, 14, 15, 12, 18, 11, 10, 9, 8, 7, 5, 2, 3, 4, 1, 0, 6],
    ![12, 17, 15, 14, 16, 13, 18, 11, 9, 8, 10, 7, 0, 5, 3, 2, 4, 1, 6],
    ![12, 15, 17, 14, 16, 13, 18, 10, 9, 8, 7, 11, 0, 5, 3, 1, 4, 2, 6],
    ![15, 14, 17, 11, 16, 13, 18, 10, 9, 8, 7, 3, 12, 5, 1, 0, 4, 2, 6],
    ![16, 14, 17, 11, 15, 12, 18, 10, 9, 8, 7, 3, 5, 13, 1, 4, 0, 2, 6],
    ![16, 13, 17, 11, 15, 12, 18, 10, 9, 8, 7, 3, 5, 1, 14, 4, 0, 2, 6],
    ![14, 17, 13, 11, 16, 12, 18, 10, 9, 8, 7, 3, 5, 2, 0, 15, 4, 1, 6],
    ![14, 17, 13, 11, 15, 12, 18, 10, 9, 8, 7, 3, 5, 2, 0, 4, 16, 1, 6],
    ![14, 16, 13, 11, 15, 12, 18, 10, 9, 8, 7, 3, 5, 2, 0, 4, 1, 17, 6],
    ![12, 17, 15, 14, 16, 13, 11, 10, 9, 8, 7, 6, 0, 5, 3, 2, 4, 1, 18]
  ]

theorem sevenPartner_valid : ∀ r x : Fin 19, x ≠ r →
    sevenPartner r x ≠ r ∧
    sevenPartner r (sevenPartner r x) = x ∧
    sevenGraph.Adj x (sevenPartner r x) := by decide

/-- A fixed-point-free involution along graph edges is a perfect matching. -/
theorem hasPerfectMatching_of_partner {W : Type*} (H : SimpleGraph W)
    (p : W → W) (hinv : ∀ x, p (p x) = x)
    (hadj : ∀ x, H.Adj x (p x)) : HasPerfectMatching H := by
  let M : H.Subgraph := {
    verts := Set.univ
    Adj := fun x y => p x = y
    adj_sub := by intro x y h; rw [← h]; exact hadj x
    edge_vert := by intro x y h; trivial
    symm := ⟨by intro x y h; rw [← h]; exact hinv x⟩ }
  refine ⟨M, SimpleGraph.Subgraph.isPerfectMatching_iff.mpr ?_⟩
  intro x
  exact ⟨p x, rfl, fun y hy => (show p x = y from hy).symm⟩

theorem sevenGraph_factor_critical : IsFactorCritical sevenGraph := by
  intro r
  let p : {x : Fin 19 // x ≠ r} → {x : Fin 19 // x ≠ r} :=
    fun x => ⟨sevenPartner r x, (sevenPartner_valid r x x.property).1⟩
  apply hasPerfectMatching_of_partner _ p
  · intro x
    apply Subtype.ext
    exact (sevenPartner_valid r x x.property).2.1
  · intro x
    exact (sevenPartner_valid r x x.property).2.2

def sevenMatching : sevenGraph.Subgraph where
  verts := {x | x ≠ 18}
  Adj x y := x ≠ 18 ∧ sevenPartner 18 x = y
  adj_sub := by
    intro x y h
    rw [← h.2]
    exact (sevenPartner_valid 18 x h.1).2.2
  edge_vert := by intro x y h; exact h.1
  symm := ⟨by
    intro x y h
    rw [← h.2]
    exact ⟨(sevenPartner_valid 18 x h.1).1,
      (sevenPartner_valid 18 x h.1).2.1⟩⟩

theorem sevenMatching_isMatching : sevenMatching.IsMatching := by
  intro x hx
  exact ⟨sevenPartner 18 x, ⟨hx, rfl⟩, fun y hy => hy.2.symm⟩

theorem sevenGraph_matchingNumber : matchingNumber sevenGraph = 9 := by
  apply le_antisymm
  · unfold matchingNumber
    apply csSup_le'
    rintro k ⟨M, hM, rfl⟩
    have hc : M.verts.ncard ≤ 19 := by
      simpa [Nat.card_eq_fintype_card] using Set.ncard_le_card M.verts
    omega
  · have h := matching_size_le_matchingNumber sevenGraph sevenMatching_isMatching
    have hc : sevenMatching.verts.ncard = 18 := by
      change ({x : Fin 19 | x ≠ 18} : Set (Fin 19)).ncard = 18
      rw [Set.ncard_eq_toFinset_card']
      change (Finset.univ.filter (fun x : Fin 19 => x ≠ 18)).card = 18
      decide
    rw [hc] at h
    exact h

theorem sevenGraph_threshold_profile : HasThresholdDegreeProfile sevenGraph 7 := by
  right
  refine ⟨by decide, ?_, ?_⟩
  · intro v
    rw [sevenGraph_degree]
    split_ifs <;> omega
  · simp_rw [sevenGraph_degree]
    decide

theorem seven_isZCandidate : IsZCandidate 7 9 := by
  exact ⟨sevenGraph, sevenGraph_triangle_free, sevenGraph_factor_critical,
    sevenGraph_matchingNumber, sevenGraph_threshold_profile⟩

/-- This upper bound has no literature or initial-cycle hypotheses. -/
theorem Z_seven_le_nine : Z 7 ≤ 9 := by
  exact Nat.sInf_le seven_isZCandidate

/-- The fixed-order bound on seventeen vertices excludes the preceding
threshold candidate, with no external hypotheses. -/
theorem seven_not_candidate_eight :
    ¬ IsZCandidate 7 8 := by
  rintro ⟨G, htri, _hcrit, _hmatch, hprofile⟩
  let _ : DecidableRel G.Adj := Classical.decRel G.Adj
  rcases hprofile with ⟨heven, _⟩ | ⟨_hodd, hmax, hsum⟩
  · norm_num at heven
  have hdeg : ∀ v, degreeZ G v ≤ 7 := by
    intro v
    rw [degreeZ_eq_degree, ← degreeN_eq_degree]
    exact_mod_cast hmax v
  have hD : deficiency G 7 = 1 := by
    calc
      deficiency G 7 = ∑ v, ((7 - degreeN G v : ℕ) : ℤ) := by
        apply Finset.sum_congr rfl
        intro v _
        rw [Nat.cast_sub (hmax v), degreeN_eq_degree, degreeZ_eq_degree]
        norm_num
      _ = 1 := by rw [← Nat.cast_sum, hsum]; norm_num
  have he := fixed_order_edge_bound G htri 7 1
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    hdeg (Or.inl (by norm_num))
  have hid := deficiency_eq G 7
  norm_num only [Fintype.card_fin] at hid
  omega

/-- The exact value uses only the published general lower bound.
Its upper bound is the explicit construction above. -/
theorem Z_seven_eq_nine (hlower : 8 ≤ Z 7) :
    Z 7 = 9 := by
  have hupper := Z_seven_le_nine
  by_contra h
  have heq : Z 7 = 8 := by omega
  have hc := Z_isZCandidate 7 ⟨9, seven_isZCandidate⟩
  rw [heq] at hc
  exact seven_not_candidate_eight hc

end TriangleProof
