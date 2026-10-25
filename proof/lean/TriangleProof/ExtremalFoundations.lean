import TriangleProof.ExtremalDefinitions
import Mathlib.Combinatorics.SimpleGraph.VertexCover

/-!
# Foundational facts for the extremal definitions

This file proves the finiteness facts needed by the natural-number `sSup`
in `fDelta`.  It also provides the small interface used to turn a published
existence construction for `Z` into the statement that `Z d` is itself a
candidate.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V]
variable (G : SimpleGraph V)

/-- The set of cardinalities used to define the matching number is bounded. -/
theorem matchingSizes_bddAbove :
    BddAbove {k | ∃ M : G.Subgraph, M.IsMatching ∧ k = M.verts.ncard / 2} := by
  refine ⟨Fintype.card V / 2, ?_⟩
  rintro k ⟨M, _hM, rfl⟩
  exact Nat.div_le_div_right (by
    simpa [Nat.card_eq_fintype_card] using Set.ncard_le_card M.verts)

/-- Every matching has at most `matchingNumber G` edges. -/
theorem matching_size_le_matchingNumber {M : G.Subgraph} (hM : M.IsMatching) :
    M.verts.ncard / 2 ≤ matchingNumber G := by
  unfold matchingNumber
  exact le_csSup (matchingSizes_bddAbove G) ⟨M, hM, rfl⟩

/-- A finite graph has a matching maximal under subgraph inclusion. -/
theorem exists_maximal_matching :
    ∃ M : G.Subgraph, M.IsMatching ∧
      Maximal (fun N : G.Subgraph => N.IsMatching) M := by
  obtain ⟨M, _hbotM, hmax⟩ :=
    @Finite.exists_le_maximal G.Subgraph _ _
      (fun N : G.Subgraph => N.IsMatching) ⊥ (by
        intro v hv
        simp at hv)
  exact ⟨M, hmax.prop, hmax⟩

omit [Fintype V] in
/-- The endpoints of a maximal matching form a vertex cover. -/
theorem maximal_matching_isVertexCover {M : G.Subgraph}
    (hM : M.IsMatching)
    (hmax : Maximal (fun N : G.Subgraph => N.IsMatching) M) :
    G.IsVertexCover M.verts := by
  intro u v huv
  by_contra hnot
  push Not at hnot
  let E := G.subgraphOfAdj huv
  have hE : E.IsMatching := SimpleGraph.Subgraph.IsMatching.subgraphOfAdj huv
  have hdisj : Disjoint M.support E.support := by
    rw [hM.support_eq_verts, hE.support_eq_verts, Set.disjoint_left]
    intro x hxM hxE
    simp only [E, SimpleGraph.subgraphOfAdj_verts, Set.mem_insert_iff,
      Set.mem_singleton_iff] at hxE
    rcases hxE with rfl | rfl
    · exact hnot.1 hxM
    · exact hnot.2 hxM
  have hsup : (M ⊔ E).IsMatching := hM.sup hE hdisj
  have heq : M = M ⊔ E := hmax.eq_of_le hsup le_sup_left
  have hEM : E ≤ M := by rw [heq]; exact le_sup_right
  have hMuv : M.Adj u v := (SimpleGraph.subgraphOfAdj_le_iff huv M).mp hEM
  exact hnot.1 hMuv.fst_mem

/-- `degreeN` agrees with Mathlib's finite degree. -/
theorem degreeN_eq_degree [DecidableRel G.Adj] (v : V) :
    degreeN G v = G.degree v := by
  rw [degreeN, Set.ncard_eq_toFinset_card']
  change (G.neighborFinset v).card = G.degree v
  rfl

/-- A vertex cover of size `c` in a graph of maximum degree `d` covers at
most `c*d` edges. -/
theorem edge_count_le_cover_mul_degree [DecidableEq V] [DecidableRel G.Adj]
    (c : Set V) (hc : G.IsVertexCover c) (d : ℕ)
    (hdeg : ∀ v, degreeN G v ≤ d) :
    G.edgeSet.ncard ≤ c.ncard * d := by
  let _ : Fintype c := Fintype.ofFinite _
  have hsubset : G.edgeFinset ⊆
      c.toFinset.biUnion (fun v => G.incidenceFinset v) := by
    intro e he
    revert he
    refine Sym2.inductionOn e ?_
    intro u v he
    have huv : G.Adj u v := by simpa using he
    rcases hc huv with hu | hv
    · simp only [Finset.mem_biUnion]
      refine ⟨u, by simpa using hu, ?_⟩
      rw [G.incidenceFinset_eq_filter]
      simp [he]
    · simp only [Finset.mem_biUnion]
      refine ⟨v, by simpa using hv, ?_⟩
      rw [G.incidenceFinset_eq_filter]
      simp [he]
  rw [Set.ncard_eq_toFinset_card']
  change G.edgeFinset.card ≤ c.ncard * d
  calc
    G.edgeFinset.card ≤
        (c.toFinset.biUnion (fun v => G.incidenceFinset v)).card :=
      Finset.card_le_card hsubset
    _ ≤ ∑ v ∈ c.toFinset, (G.incidenceFinset v).card :=
      Finset.card_biUnion_le
    _ ≤ ∑ _v ∈ c.toFinset, d := by
      gcongr with v hv
      rw [G.card_incidenceFinset_eq_degree, ← degreeN_eq_degree G v]
      exact hdeg v
    _ = c.ncard * d := by
      rw [Finset.sum_const_nat (fun _ _ => rfl)]
      rw [← Set.ncard_eq_toFinset_card' c]

/-- The endpoint cover of a maximal matching gives a uniform edge bound. -/
theorem edge_count_le_of_matchingNumber_le (d m : ℕ)
    (hdeg : HasMaximumDegreeAtMost G d)
    (hmatch : matchingNumber G ≤ m) :
    G.edgeSet.ncard ≤ 2 * m * d := by
  let _ : DecidableEq V := Classical.decEq V
  let _ : DecidableRel G.Adj := Classical.decRel G.Adj
  obtain ⟨M, hM, hmax⟩ := exists_maximal_matching G
  let _ : Fintype M.verts := Fintype.ofFinite _
  have hcover := maximal_matching_isVertexCover G hM hmax
  have hmatchingSize : M.verts.ncard / 2 ≤ m :=
    (matching_size_le_matchingNumber G hM).trans hmatch
  have hEven : Even M.verts.ncard := by
    rw [Set.ncard_eq_toFinset_card']
    exact hM.even_card
  obtain ⟨q, hq⟩ := hEven
  have hcoverCard : M.verts.ncard ≤ 2 * m := by omega
  exact (edge_count_le_cover_mul_degree G M.verts hcover d hdeg).trans
    (Nat.mul_le_mul_right d hcoverCard)

/-- The attainable edge counts in the definition of `fDelta` are bounded. -/
theorem attainableEdgeCounts_bddAbove (d m : ℕ) :
    BddAbove (attainableEdgeCounts d m) := by
  refine ⟨2 * m * d, ?_⟩
  rintro e ⟨n, G, hG, rfl⟩
  exact edge_count_le_of_matchingNumber_le G d m hG.2.1 hG.2.2

/-- Every attainable value lies below `fDelta`. -/
theorem attainableEdgeCount_le_fDelta {d m e : ℕ}
    (he : e ∈ attainableEdgeCounts d m) : e ≤ fDelta d m := by
  unfold fDelta
  exact le_csSup (attainableEdgeCounts_bddAbove d m) he

/-- A concrete universal upper bound for the extremal function. -/
theorem fDelta_le_general_bound (d m : ℕ) :
    fDelta d m ≤ 2 * m * d := by
  unfold fDelta
  apply csSup_le'
  rintro e ⟨n, G, hG, rfl⟩
  exact edge_count_le_of_matchingNumber_le G d m hG.2.1 hG.2.2

omit [Fintype V] in
/-- The empty graph has matching number zero. -/
theorem matchingNumber_bot : matchingNumber (⊥ : SimpleGraph V) = 0 := by
  have hle : matchingNumber (⊥ : SimpleGraph V) ≤ 0 := by
    unfold matchingNumber
    apply csSup_le'
    rintro k ⟨M, hM, rfl⟩
    have hverts : M.verts = ∅ := by
      ext v
      simp only [Set.mem_empty_iff_false, iff_false]
      intro hv
      obtain ⟨w, hw, _⟩ := hM hv
      have := M.adj_sub hw
      simp at this
    simp [hverts]
  omega

/-- Zero is always attainable, using the empty graph. -/
theorem attainableEdgeCounts_nonempty (d m : ℕ) :
    (attainableEdgeCounts d m).Nonempty := by
  refine ⟨0, 0, (⊥ : SimpleGraph (Fin 0)), ?_, by simp⟩
  refine ⟨?_, ?_, ?_⟩
  · intro a
    exact Fin.elim0 a
  · intro v
    exact Fin.elim0 v
  · rw [matchingNumber_bot]
    omega

/-- Thus `fDelta` is an attained maximum, rather than a fallback value of
the natural-number supremum. -/
theorem fDelta_is_attainable (d m : ℕ) :
    fDelta d m ∈ attainableEdgeCounts d m := by
  unfold fDelta
  exact Nat.sSup_mem (attainableEdgeCounts_nonempty d m)
    (attainableEdgeCounts_bddAbove d m)

/-- A published threshold construction makes the set defining `Z` nonempty. -/
theorem zCandidateSet_nonempty (d : ℕ)
    (hexists : ∃ i, IsZCandidate d i) :
    {i | IsZCandidate d i}.Nonempty := hexists

/-- Once a threshold construction exists, the natural-number infimum is
itself a candidate and therefore is the intended minimum. -/
theorem Z_isZCandidate (d : ℕ) (hexists : ∃ i, IsZCandidate d i) :
    IsZCandidate d (Z d) := by
  unfold Z
  exact Nat.sInf_mem (zCandidateSet_nonempty d hexists)

end TriangleProof
