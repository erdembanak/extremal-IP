import TriangleProof.ExtremalFoundations

/-!
# The intermediate construction from Proposition 4.1

This file formalizes the graph `B_{d,d+t}` from Proposition 4.1 of
Banak--Ekim--Taşkın (2023).  It starts with two parts `F ∪ H`, where
`|F| = d - 1` and `|H| = t + 1`, deletes `t` cyclic perfect matchings
between the two copies of `F`, deletes all `H`--`H` edges, and joins a new
apex to both copies of `H`.
-/

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

/-- The five vertex classes of `B_{d,d+t}`. -/
inductive IntermediateVertex (d t : ℕ) where
  | apex
  | leftF (x : Fin (d - 1))
  | leftH (x : Fin (t + 1))
  | rightF (x : Fin (d - 1))
  | rightH (x : Fin (t + 1))
  deriving DecidableEq, Fintype

namespace IntermediateConstruction

/-- The nonzero cyclic shifts `1,...,t` are the deleted matchings. -/
def IsDeletedShift {d t : ℕ} (x : Fin (d - 1)) : Prop :=
  0 < x.val ∧ x.val ≤ t

instance {d t : ℕ} (x : Fin (d - 1)) : Decidable (IsDeletedShift (t := t) x) :=
  inferInstanceAs (Decidable (_ ∧ _))

/-- An `F`--`F` edge belongs to one of the deleted cyclic matchings. -/
def IsDeletedEdge {d t : ℕ} (a b : Fin (d - 1)) : Prop :=
  IsDeletedShift (t := t) (b - a)

instance {d t : ℕ} (a b : Fin (d - 1)) : Decidable (IsDeletedEdge (t := t) a b) :=
  inferInstanceAs (Decidable (IsDeletedShift (t := t) (b - a)))

private def Rel (d t : ℕ) : IntermediateVertex d t → IntermediateVertex d t → Prop
  | .leftF a, .rightF b => ¬ IsDeletedEdge (t := t) a b
  | .rightF b, .leftF a => ¬ IsDeletedEdge (t := t) a b
  | .leftF _, .rightH _ => True
  | .rightH _, .leftF _ => True
  | .leftH _, .rightF _ => True
  | .rightF _, .leftH _ => True
  | .apex, .leftH _ => True
  | .leftH _, .apex => True
  | .apex, .rightH _ => True
  | .rightH _, .apex => True
  | _, _ => False

private instance (d t : ℕ) : DecidableRel (Rel d t) := by
  intro a b
  cases a <;> cases b <;> simp only [Rel] <;> infer_instance

/-- The graph `B_{d,d+t}` used in Proposition 4.1. -/
def graph (d t : ℕ) : SimpleGraph (IntermediateVertex d t) where
  Adj := Rel d t
  symm := ⟨by intro a b; cases a <;> cases b <;> simp [Rel]⟩
  loopless := ⟨by intro a; cases a <;> simp [Rel]⟩

instance (d t : ℕ) : DecidableRel (graph d t).Adj :=
  inferInstanceAs (DecidableRel (Rel d t))

@[simp] theorem graph_adj_iff (d t : ℕ) (a b : IntermediateVertex d t) :
    (graph d t).Adj a b ↔ Rel d t a b := Iff.rfl

private def vertexEquiv (d t : ℕ) : IntermediateVertex d t ≃
    Unit ⊕ (Fin (d - 1) ⊕ (Fin (t + 1) ⊕ (Fin (d - 1) ⊕ Fin (t + 1)))) where
  toFun
    | .apex => .inl ()
    | .leftF x => .inr (.inl x)
    | .leftH x => .inr (.inr (.inl x))
    | .rightF x => .inr (.inr (.inr (.inl x)))
    | .rightH x => .inr (.inr (.inr (.inr x)))
  invFun
    | .inl _ => .apex
    | .inr (.inl x) => .leftF x
    | .inr (.inr (.inl x)) => .leftH x
    | .inr (.inr (.inr (.inl x))) => .rightF x
    | .inr (.inr (.inr (.inr x))) => .rightH x
  left_inv := by intro x; cases x <;> rfl
  right_inv := by intro x; rcases x with _ | (_ | (_ | (_ | _))) <;> rfl

/-- The vertex count is `2(d+t)+1`. -/
theorem card_vertices {d t : ℕ} (hd : 1 ≤ d) :
    Fintype.card (IntermediateVertex d t) = 2 * (d + t) + 1 := by
  rw [Fintype.card_congr (vertexEquiv d t)]
  simp
  omega

/-- There are exactly `t` deleted shifts. -/
theorem card_deleted_shifts {d t : ℕ} (ht : t < d - 1) :
    Fintype.card {x : Fin (d - 1) // IsDeletedShift (t := t) x} = t := by
  let e : Fin t ≃ {x : Fin (d - 1) // IsDeletedShift (t := t) x} := {
    toFun := fun x => ⟨⟨x.val + 1, by omega⟩, by
      constructor <;> simp⟩
    invFun := fun x => ⟨x.val.val - 1, by
      have := x.property
      simp only [IsDeletedShift] at this
      omega⟩
    left_inv := by
      intro x
      apply Fin.ext
      change (x.val + 1 - 1) = x.val
      omega
    right_inv := by
      intro x
      apply Subtype.ext
      apply Fin.ext
      have := x.property
      simp only [IsDeletedShift] at this
      change (x.val.val - 1 + 1) = x.val.val
      omega }
  simpa using (Fintype.card_congr e).symm

/-- Translation shows that each `F` vertex loses exactly `t` neighbors. -/
theorem card_deleted_neighbors {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1)
    (a : Fin (d - 1)) :
    Fintype.card {b : Fin (d - 1) // IsDeletedEdge (t := t) a b} = t := by
  letI : NeZero (d - 1) := ⟨by omega⟩
  let e := (Equiv.subRight a).subtypeEquiv
    (p := fun b : Fin (d - 1) => IsDeletedEdge (t := t) a b)
    (q := fun x : Fin (d - 1) => IsDeletedShift (t := t) x)
    (by intro b; rfl)
  rw [Fintype.card_congr e]
  exact card_deleted_shifts ht

/-- Consequently every `F` vertex retains `d-1-t` neighbors in the
opposite `F` part. -/
theorem card_retained_neighbors {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1)
    (a : Fin (d - 1)) :
    Fintype.card {b : Fin (d - 1) // ¬ IsDeletedEdge (t := t) a b} = d - 1 - t := by
  rw [Fintype.card_subtype_compl]
  rw [card_deleted_neighbors hd ht a]
  simp

/-- The same count holds when the right endpoint is fixed. -/
theorem card_deleted_left_neighbors {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1)
    (b : Fin (d - 1)) :
    Fintype.card {a : Fin (d - 1) // IsDeletedEdge (t := t) a b} = t := by
  letI : NeZero (d - 1) := ⟨by omega⟩
  let e := (Equiv.subLeft b).subtypeEquiv
    (p := fun a : Fin (d - 1) => IsDeletedEdge (t := t) a b)
    (q := fun x : Fin (d - 1) => IsDeletedShift (t := t) x)
    (by intro a; rfl)
  rw [Fintype.card_congr e]
  exact card_deleted_shifts ht

theorem card_retained_left_neighbors {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1)
    (b : Fin (d - 1)) :
    Fintype.card {a : Fin (d - 1) // ¬ IsDeletedEdge (t := t) a b} = d - 1 - t := by
  rw [Fintype.card_subtype_compl]
  rw [card_deleted_left_neighbors hd ht b]
  simp

private def apexNeighborEquiv (d t : ℕ) :
    {v // (graph d t).Adj (.apex) v} ≃ Fin (t + 1) ⊕ Fin (t + 1) where
  toFun
    | ⟨.leftH x, _⟩ => .inl x
    | ⟨.rightH x, _⟩ => .inr x
  invFun
    | .inl x => ⟨.leftH x, by simp [graph, Rel]⟩
    | .inr x => ⟨.rightH x, by simp [graph, Rel]⟩
  left_inv := by rintro ⟨v, hv⟩; cases v <;> simp [graph, Rel] at hv ⊢
  right_inv := by intro x; cases x <;> rfl

private def leftHNeighborEquiv (d t : ℕ) :
    {v // (graph d t).Adj (.leftH (0 : Fin (t + 1))) v} ≃
      Option (Fin (d - 1)) where
  toFun
    | ⟨.apex, _⟩ => none
    | ⟨.rightF x, _⟩ => some x
  invFun
    | none => ⟨.apex, by simp [graph, Rel]⟩
    | some x => ⟨.rightF x, by simp [graph, Rel]⟩
  left_inv := by rintro ⟨v, hv⟩; cases v <;> simp [graph, Rel] at hv ⊢
  right_inv := by intro x; cases x <;> rfl

private def rightHNeighborEquiv (d t : ℕ) :
    {v // (graph d t).Adj (.rightH (0 : Fin (t + 1))) v} ≃
      Option (Fin (d - 1)) where
  toFun
    | ⟨.apex, _⟩ => none
    | ⟨.leftF x, _⟩ => some x
  invFun
    | none => ⟨.apex, by simp [graph, Rel]⟩
    | some x => ⟨.leftF x, by simp [graph, Rel]⟩
  left_inv := by rintro ⟨v, hv⟩; cases v <;> simp [graph, Rel] at hv ⊢
  right_inv := by intro x; cases x <;> rfl

private def leftFNeighborEquiv (d t : ℕ) (a : Fin (d - 1)) :
    {v // (graph d t).Adj (.leftF a) v} ≃
      Fin (t + 1) ⊕ {b : Fin (d - 1) // ¬ IsDeletedEdge (t := t) a b} where
  toFun
    | ⟨.rightF b, hb⟩ => .inr ⟨b, hb⟩
    | ⟨.rightH h, _⟩ => .inl h
  invFun
    | .inl h => ⟨.rightH h, by simp [graph, Rel]⟩
    | .inr b => ⟨.rightF b, b.property⟩
  left_inv := by rintro ⟨v, hv⟩; cases v <;> simp [graph, Rel] at hv ⊢
  right_inv := by intro x; rcases x with x | x <;> rfl

private def rightFNeighborEquiv (d t : ℕ) (b : Fin (d - 1)) :
    {v // (graph d t).Adj (.rightF b) v} ≃
      Fin (t + 1) ⊕ {a : Fin (d - 1) // ¬ IsDeletedEdge (t := t) a b} where
  toFun
    | ⟨.leftF a, ha⟩ => .inr ⟨a, ha⟩
    | ⟨.leftH h, _⟩ => .inl h
  invFun
    | .inl h => ⟨.leftH h, by simp [graph, Rel]⟩
    | .inr a => ⟨.leftF a, a.property⟩
  left_inv := by rintro ⟨v, hv⟩; cases v <;> simp [graph, Rel] at hv ⊢
  right_inv := by intro x; rcases x with x | x <;> rfl

/-- All non-apex vertices have degree `d`; the apex has degree `2(t+1)`. -/
theorem graph_degree {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1)
    (v : IntermediateVertex d t) :
    degreeN (graph d t) v = match v with
      | .apex => 2 * (t + 1)
      | _ => d := by
  rw [degreeN, ← Set.fintypeCard_eq_ncard]
  change Fintype.card {w // (graph d t).Adj v w} = _
  cases v with
  | apex =>
      rw [Fintype.card_congr (apexNeighborEquiv d t)]
      simp
      omega
  | leftF a =>
      rw [Fintype.card_congr (leftFNeighborEquiv d t a)]
      simp [card_retained_neighbors hd ht a]
      omega
  | leftH h =>
      let e : {v // (graph d t).Adj (.leftH h) v} ≃
          {v // (graph d t).Adj (.leftH (0 : Fin (t + 1))) v} :=
        Equiv.subtypeEquiv (Equiv.refl _) (by
          intro v
          cases v <;> simp [graph, Rel])
      rw [Fintype.card_congr e, Fintype.card_congr (leftHNeighborEquiv d t)]
      simp
      omega
  | rightF b =>
      rw [Fintype.card_congr (rightFNeighborEquiv d t b)]
      simp [card_retained_left_neighbors hd ht b]
      omega
  | rightH h =>
      let e : {v // (graph d t).Adj (.rightH h) v} ≃
          {v // (graph d t).Adj (.rightH (0 : Fin (t + 1))) v} :=
        Equiv.subtypeEquiv (Equiv.refl _) (by
          intro v
          cases v <;> simp [graph, Rel])
      rw [Fintype.card_congr e, Fintype.card_congr (rightHNeighborEquiv d t)]
      simp
      omega

/-- Under the numerical hypothesis in Proposition 4.1, the maximum degree
is at most `d`. -/
theorem graph_maximum_degree {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    HasMaximumDegreeAtMost (graph d t) d := by
  intro v
  rw [graph_degree hd (by omega) v]
  cases v <;> simp <;> omega

/-- The maximum degree is exactly `d`, as stated in Proposition 4.1. -/
theorem graph_maxDegree {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    (graph d t).maxDegree = d := by
  apply le_antisymm
  · apply SimpleGraph.maxDegree_le_of_forall_degree_le
    intro v
    rw [← degreeN_eq_degree]
    exact graph_maximum_degree hd ht v
  · have h := SimpleGraph.degree_le_maxDegree (graph d t)
      (.leftH (0 : Fin (t + 1)))
    rw [← degreeN_eq_degree, graph_degree hd (by omega)] at h
    exact h

/-- The edge count stated in Proposition 4.1. -/
theorem graph_edge_count {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    (graph d t).edgeSet.ncard = d * (d + t) + t + 1 := by
  let G := graph d t
  let a : IntermediateVertex d t := .apex
  have hdegree (v : IntermediateVertex d t) :
      G.degree v = if v = a then 2 * (t + 1) else d := by
    rw [← degreeN_eq_degree]
    rw [graph_degree hd (by omega) v]
    cases v <;> simp [a]
  have hsum : ∑ v, G.degree v = 2 * (d * (d + t) + t + 1) := by
    classical
    rw [← Finset.sum_erase_add Finset.univ (fun v => G.degree v)
      (Finset.mem_univ a)]
    have hnon : ∀ v ∈ (Finset.univ.erase a), G.degree v = d := by
      intro v hv
      rw [hdegree]
      have hvne : v ≠ a := by simpa using hv
      simp [hvne]
    rw [Finset.sum_congr rfl hnon]
    rw [Finset.sum_const_nat (fun _ _ => rfl)]
    change (Finset.univ.erase a).card * d + G.degree a = _
    rw [hdegree]
    simp only [if_pos, Finset.card_erase_of_mem, Finset.mem_univ,
      Finset.card_univ]
    rw [card_vertices (by omega)]
    simp
    ring
  have htwice := SimpleGraph.sum_degrees_eq_twice_card_edges G
  rw [hsum] at htwice
  rw [Set.ncard_eq_toFinset_card']
  change G.edgeFinset.card = _
  omega

/-- The order alone bounds every matching by `d+t`. -/
theorem graph_matchingNumber_le {d t : ℕ} (hd : 1 ≤ d) :
    matchingNumber (graph d t) ≤ d + t := by
  unfold matchingNumber
  apply csSup_le'
  rintro k ⟨M, _hM, rfl⟩
  have hc : M.verts.ncard ≤ 2 * (d + t) + 1 := by
    rw [← card_vertices hd]
    simpa [Nat.card_eq_fintype_card] using Set.ncard_le_card M.verts
  omega

private def smallF {d t : ℕ} (ht : t < d - 1) (x : Fin t) : Fin (d - 1) :=
  Fin.castLE (Nat.le_of_lt ht) x

/-- A fixed-point-free involution on all vertices except one `H` vertex.
Its pairs give the matching of size `d+t` described in the proof of
Proposition 4.1. -/
def partner {d t : ℕ} (ht : t < d - 1) :
    IntermediateVertex d t → IntermediateVertex d t
  | .apex => .leftH 0
  | .leftH h => Fin.cases .apex (fun x => .rightF (smallF ht x)) h
  | .rightH h => Fin.cases (.rightH 0) (fun x => .leftF (smallF ht x)) h
  | .leftF a => if ha : a.val < t then
      .rightH (Fin.succ ⟨a.val, ha⟩) else .rightF a
  | .rightF b => if hb : b.val < t then
      .leftH (Fin.succ ⟨b.val, hb⟩) else .leftF b

def unmatched (d t : ℕ) : IntermediateVertex d t := .rightH 0

theorem partner_involutive {d t : ℕ} (ht : t < d - 1) :
    ∀ x, partner ht (partner ht x) = x := by
  intro x
  cases x with
  | apex => simp [partner]
  | leftH h =>
      refine Fin.cases ?_ (fun x => ?_) h
      · simp [partner]
      · simp [partner, smallF, x.isLt]
  | rightH h =>
      refine Fin.cases ?_ (fun x => ?_) h
      · simp [partner]
      · simp [partner, smallF, x.isLt]
  | leftF a =>
      by_cases ha : a.val < t
      · simp [partner, ha, smallF]
      · simp [partner, ha]
  | rightF b =>
      by_cases hb : b.val < t
      · simp [partner, hb, smallF]
      · simp [partner, hb]

theorem partner_avoids_unmatched {d t : ℕ} (ht : t < d - 1) :
    ∀ x, x ≠ unmatched d t → partner ht x ≠ unmatched d t := by
  intro x hx
  cases x with
  | apex => simp [partner, unmatched]
  | leftF a =>
      by_cases ha : a.val < t <;> simp [partner, unmatched, ha]
  | leftH h =>
      refine Fin.cases ?_ (fun x => ?_) h <;> simp [partner, unmatched]
  | rightF b =>
      by_cases hb : b.val < t <;> simp [partner, unmatched, hb]
  | rightH h =>
      revert hx
      refine Fin.cases ?_ (fun x => ?_) h
      · intro hx
        exact (hx rfl).elim
      · intro _hx
        simp [partner, unmatched]

theorem partner_adjacent {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1) :
    ∀ x, x ≠ unmatched d t → (graph d t).Adj x (partner ht x) := by
  intro x hx
  cases x with
  | apex => simp [partner, graph, Rel]
  | leftH h =>
      refine Fin.cases ?_ (fun x => ?_) h <;> simp [partner, graph, Rel]
  | rightH h =>
      revert hx
      refine Fin.cases ?_ (fun x => ?_) h
      · intro hx
        exact (hx rfl).elim
      · intro _hx
        simp [partner, graph, Rel]
  | leftF a =>
      by_cases ha : a.val < t
      · simp [partner, ha, graph, Rel]
      · haveI : NeZero (d - 1) := ⟨by omega⟩
        simp [partner, ha, graph, Rel, IsDeletedEdge, IsDeletedShift]
  | rightF b =>
      by_cases hb : b.val < t
      · simp [partner, hb, graph, Rel]
      · haveI : NeZero (d - 1) := ⟨by omega⟩
        simp [partner, hb, graph, Rel, IsDeletedEdge, IsDeletedShift]

def largeMatching {d t : ℕ} (ht : t < d - 1) : (graph d t).Subgraph where
  verts := {x | x ≠ unmatched d t}
  Adj x y := x ≠ unmatched d t ∧ partner ht x = y
  adj_sub := by
    intro x y h
    rw [← h.2]
    exact partner_adjacent (by omega) ht x h.1
  edge_vert := by intro x y h; exact h.1
  symm := ⟨by
    intro x y h
    rw [← h.2]
    exact ⟨partner_avoids_unmatched ht x h.1, partner_involutive ht x⟩⟩

theorem largeMatching_isMatching {d t : ℕ} (ht : t < d - 1) :
    (largeMatching ht).IsMatching := by
  intro x hx
  exact ⟨partner ht x, ⟨hx, rfl⟩, fun y hy => hy.2.symm⟩

theorem largeMatching_vertex_count {d t : ℕ} (hd : 1 ≤ d) (ht : t < d - 1) :
    (largeMatching ht).verts.ncard = 2 * (d + t) := by
  change ({x | x ≠ unmatched d t} : Set (IntermediateVertex d t)).ncard = _
  letI : Fintype {x | x ≠ unmatched d t} := Fintype.ofFinite _
  rw [Set.ncard_eq_toFinset_card']
  rw [Set.toFinset_ofPred]
  rw [Finset.filter_ne']
  rw [Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ]
  rw [card_vertices hd]
  omega

/-- The matching number part of Proposition 4.1. -/
theorem graph_matchingNumber {d t : ℕ} (hd : 2 ≤ d) (ht : t < d - 1) :
    matchingNumber (graph d t) = d + t := by
  apply le_antisymm (graph_matchingNumber_le (by omega))
  have h := matching_size_le_matchingNumber (graph d t) (largeMatching_isMatching ht)
  rw [largeMatching_vertex_count (by omega) ht] at h
  omega

/-- No three vertices of the construction form a triangle. -/
theorem graph_triangle_free (d t : ℕ) : IsTriangleFree (graph d t) := by
  intro a b c
  cases a <;> cases b <;> cases c <;> simp [graph, Rel]

/-- The graph-theoretic statement of Proposition 4.1, before relabeling the
vertices onto a canonical `Fin` type. -/
theorem proposition41_graph {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    IsTriangleFree (graph d t) ∧
    (graph d t).maxDegree = d ∧
    matchingNumber (graph d t) = d + t ∧
    (graph d t).edgeSet.ncard = d * (d + t) + t + 1 := by
  exact ⟨graph_triangle_free d t, graph_maxDegree hd ht,
    graph_matchingNumber hd (by omega), graph_edge_count hd ht⟩

/-- A canonical relabeling onto the `Fin` type used in the definition of
`fDelta`. -/
noncomputable def vertexFinEquiv {d t : ℕ} (hd : 1 ≤ d) :
    IntermediateVertex d t ≃ Fin (2 * (d + t) + 1) :=
  (Fintype.equivFin (IntermediateVertex d t)).trans (finCongr (card_vertices hd))

noncomputable def finGraph {d t : ℕ} (hd : 1 ≤ d) :
    SimpleGraph (Fin (2 * (d + t) + 1)) :=
  (graph d t).map (vertexFinEquiv hd)

theorem finGraph_triangle_free {d t : ℕ} (hd : 1 ≤ d) :
    IsTriangleFree (finGraph (t := t) hd) := by
  let e := SimpleGraph.Iso.map (vertexFinEquiv hd) (graph d t)
  intro a b c habc
  apply graph_triangle_free d t (e.symm a) (e.symm b) (e.symm c)
  have hab : (graph d t).Adj (e.symm a) (e.symm b) := by
    apply e.map_adj_iff.mp
    simpa [e, finGraph] using habc.1
  have hbc : (graph d t).Adj (e.symm b) (e.symm c) := by
    apply e.map_adj_iff.mp
    simpa [e, finGraph] using habc.2.1
  have hca : (graph d t).Adj (e.symm c) (e.symm a) := by
    apply e.map_adj_iff.mp
    simpa [e, finGraph] using habc.2.2
  exact ⟨hab, hbc, hca⟩

theorem finGraph_maximum_degree {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    HasMaximumDegreeAtMost (finGraph (t := t) (Nat.le_trans (by decide) hd)) d := by
  classical
  unfold finGraph
  let e := SimpleGraph.Iso.map
    (vertexFinEquiv (t := t) (Nat.le_trans (by decide) hd)) (graph d t)
  intro v
  obtain ⟨u, rfl⟩ := e.surjective v
  rw [degreeN_eq_degree, e.degree_eq, ← degreeN_eq_degree]
  exact graph_maximum_degree hd ht u

theorem finGraph_edge_count {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    (finGraph (t := t) (Nat.le_trans (by decide) hd)).edgeSet.ncard =
      d * (d + t) + t + 1 := by
  classical
  unfold finGraph
  let e := (vertexFinEquiv (t := t) (Nat.le_trans (by decide) hd)).toEmbedding
  let _ : DecidableRel ((graph d t).map e).Adj := inferInstance
  rw [Set.ncard_eq_toFinset_card']
  change ((graph d t).map e).edgeFinset.card = _
  calc
    _ = (graph d t).edgeFinset.card := by
      exact SimpleGraph.card_edgeFinset_map e (graph d t)
    _ = d * (d + t) + t + 1 := by
      have hs := graph_edge_count hd ht
      rw [Set.ncard_eq_toFinset_card'] at hs
      exact hs

theorem finGraph_matchingNumber_le {d t : ℕ} (hd : 1 ≤ d) :
    matchingNumber (finGraph (t := t) hd) ≤ d + t := by
  unfold matchingNumber
  apply csSup_le'
  rintro k ⟨M, _hM, rfl⟩
  have hc : M.verts.ncard ≤ 2 * (d + t) + 1 := by
    simpa [Nat.card_eq_fintype_card] using Set.ncard_le_card M.verts
  omega

/-- Proposition 4.1 as the lower bound needed by Conjecture 2.1.  The graph,
triangle-free proof, degree proof, matching bound, and edge count are all
constructed above; no existence statement from the paper is assumed. -/
theorem proposition41_lower_bound {d t : ℕ} (hd : 2 ≤ d) (ht : 2 * t + 2 ≤ d) :
    d * (d + t) + (d + t) - d + 1 ≤ fDelta d (d + t) := by
  apply attainableEdgeCount_le_fDelta
  let hd1 : 1 ≤ d := Nat.le_trans (by decide) hd
  refine ⟨2 * (d + t) + 1, finGraph (t := t) hd1, ?_, ?_⟩
  · exact ⟨finGraph_triangle_free (t := t) hd1,
      finGraph_maximum_degree hd ht,
      finGraph_matchingNumber_le (t := t) hd1⟩
  · rw [finGraph_edge_count hd ht]
    omega

end IntermediateConstruction

end TriangleProof
