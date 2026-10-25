import TriangleProof.SaturatedCycle
import Mathlib.Combinatorics.SimpleGraph.Coloring.Constructions

/-! A shortest odd closed walk is a cycle. Its minimality also controls
the neighbors of every vertex on that cycle. -/

namespace TriangleProof

open SimpleGraph

variable {V : Type*} {G : SimpleGraph V}

def IsShortestOddWalk {u : V} (p : G.Walk u u) : Prop :=
  Odd p.length ∧ ∀ v (w : G.Walk v v), Odd w.length → p.length ≤ w.length

theorem exists_shortest_odd_walk (h : ¬ G.IsBipartite) :
    ∃ u, ∃ p : G.Walk u u, IsShortestOddWalk p := by
  classical
  have hex : ∃ n : ℕ, ∃ u, ∃ p : G.Walk u u, Odd p.length ∧ p.length = n := by
    rw [SimpleGraph.IsBipartite, two_colorable_iff_forall_loop_even] at h
    push Not at h
    obtain ⟨u, p, hp⟩ := h
    exact ⟨p.length, u, p, Nat.not_even_iff_odd.mp hp, rfl⟩
  obtain ⟨u, p, hp, heq⟩ := Nat.find_spec hex
  refine ⟨u, p, hp, ?_⟩
  intro v w hw
  rw [heq]
  exact Nat.find_min' hex ⟨v, w, hw, rfl⟩

private def segment {u v : V} (p : G.Walk u v) (i j : ℕ) (hij : i ≤ j) :
    G.Walk (p.getVert i) (p.getVert j) :=
  ((p.drop i).take (j - i)).copy rfl (by rw [Walk.drop_getVert]; congr 1; omega)

private theorem segment_length {u v : V} (p : G.Walk u v)
    (i j : ℕ) (hij : i ≤ j) (hj : j ≤ p.length) :
    (segment p i j hij).length = j - i := by
  simp [segment, Walk.take_length, Walk.drop_length, Nat.min_eq_left (by omega : j-i ≤ p.length-i)]

private def complement {u : V} (p : G.Walk u u) (i j : ℕ) :
    G.Walk (p.getVert j) (p.getVert i) := (p.drop j).append (p.take i)

private theorem complement_length {u : V} (p : G.Walk u u)
    (i j : ℕ) (hij : i ≤ j) (hj : j ≤ p.length) :
    (complement p i j).length = p.length - (j - i) := by
  simp only [complement, Walk.length_append, Walk.drop_length, Walk.take_length]
  rw [Nat.min_eq_left (by omega : i ≤ p.length)]
  omega

theorem shortest_odd_walk_repeat {u : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p) {i j : ℕ} (hij : i ≤ j) (hj : j ≤ p.length)
    (heq : p.getVert i = p.getVert j) : i = j ∨ (i = 0 ∧ j = p.length) := by
  have ho := Nat.odd_iff.mp hp.1
  by_cases ha : Odd (j-i)
  · have hw := hp.2 _ ((segment p i j hij).copy rfl heq.symm)
      (by simpa only [Walk.length_copy, segment_length p i j hij hj] using ha)
    simp only [Walk.length_copy, segment_length p i j hij hj] at hw
    omega
  · have ha' : (j-i)%2 ≠ 1 := by simpa only [Nat.odd_iff] using ha
    have hb : Odd (p.length-(j-i)) := Nat.odd_iff.mpr (by omega)
    have hw := hp.2 _ ((complement p i j).copy rfl heq)
      (by simpa only [Walk.length_copy, complement_length p i j hij hj] using hb)
    simp only [Walk.length_copy, complement_length p i j hij hj] at hw
    omega

theorem shortest_odd_walk_isCycle {u : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p) : p.IsCycle := by
  have ho := Nat.odd_iff.mp hp.1
  have hthree : 3 ≤ p.length := by
    by_contra h
    have hl : p.length = 1 := by omega
    have ha := p.adj_getVert_succ (i := 0) (by omega)
    have he : p.getVert 1 = u := by rw [← hl, Walk.getVert_length]
    simp [he] at ha
  apply Walk.isCycle_iff_isPath_tail_and_le_length.mpr
  refine ⟨(Walk.IsPath.getVert_injOn_iff _).mp ?_, hthree⟩
  intro i hi j hj heq
  simp only [Set.mem_ofPred_eq, Walk.length_tail] at hi hj
  simp only [Walk.getVert_tail] at heq
  rcases le_total i j with hij | hji
  · have hh := shortest_odd_walk_repeat hp (by omega : i+1 ≤ j+1) (by omega) heq
    omega
  · have hh := shortest_odd_walk_repeat hp (by omega : j+1 ≤ i+1) (by omega) heq.symm
    omega

/-- Non-bipartiteness supplies an odd cycle of minimum length, in fact
minimum among all odd closed walks. No finiteness assumption is needed. -/
theorem exists_shortest_odd_cycle (h : ¬ G.IsBipartite) :
    ∃ u, ∃ p : G.Walk u u, p.IsCycle ∧ IsShortestOddWalk p := by
  obtain ⟨u, p, hp⟩ := exists_shortest_odd_walk h
  exact ⟨u, p, shortest_odd_walk_isCycle hp, hp⟩

theorem shortest_odd_walk_neighbor_gap {u v : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p) {i j : ℕ} (hij : i < j) (hj : j < p.length)
    (hiadj : G.Adj v (p.getVert i)) (hjadj : G.Adj v (p.getVert j)) :
    j - i = 2 ∨ j - i = p.length - 2 := by
  have ho := Nat.odd_iff.mp hp.1
  by_cases ha : Odd (j-i)
  · have hw := hp.2 v (Walk.cons hiadj ((segment p i j (by omega)).concat hjadj.symm))
      (by simpa only [Walk.length_cons, Walk.length_concat,
        segment_length p i j (by omega) (by omega)] using ha.add_even (by decide : Even 2))
    simp only [Walk.length_cons, Walk.length_concat,
      segment_length p i j (by omega) (by omega)] at hw
    have := Nat.odd_iff.mp ha
    omega
  · have ha' : (j-i)%2 ≠ 1 := by simpa only [Nat.odd_iff] using ha
    have hb : Odd (p.length-(j-i)+2) := Nat.odd_iff.mpr (by omega)
    have hw := hp.2 v (Walk.cons hjadj ((complement p i j).concat hiadj.symm))
      (by simpa only [Walk.length_cons, Walk.length_concat,
        complement_length p i j (by omega) (by omega)] using hb)
    simp only [Walk.length_cons, Walk.length_concat,
      complement_length p i j (by omega) (by omega)] at hw
    omega

theorem shortest_odd_cycle_neighbor_bound [DecidableRel G.Adj]
    {u : V} {p : G.Walk u u} (hp : IsShortestOddWalk p) (hl : 5 ≤ p.length)
    (v : V) :
    (Finset.univ.filter fun i : Fin p.length => G.Adj v (p.getVert i.val)).card ≤ 2 := by
  have hthree : ∀ (i j k : Fin p.length), i < j → j < k →
      G.Adj v (p.getVert i.val) → G.Adj v (p.getVert j.val) →
      G.Adj v (p.getVert k.val) → False := by
    intro i j k hij hjk hi hj hk
    have h1 := shortest_odd_walk_neighbor_gap hp hij j.isLt hi hj
    have h2 := shortest_odd_walk_neighbor_gap hp hjk k.isLt hj hk
    have h3 := shortest_odd_walk_neighbor_gap hp (lt_trans hij hjk) k.isLt hi hk
    have ho := Nat.odd_iff.mp hp.1
    rcases h1 with h1 | h1 <;> rcases h2 with h2 | h2 <;>
      rcases h3 with h3 | h3 <;> omega
  by_contra hn
  obtain ⟨i, j, k, hi, hj, hk, hij, hik, hjk⟩ :=
    Finset.two_lt_card_iff.mp (by omega : 2 < (Finset.univ.filter
      fun i : Fin p.length => G.Adj v (p.getVert i.val)).card)
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hi hj hk
  have hord : (i < j ∧ j < k) ∨ (i < k ∧ k < j) ∨ (j < i ∧ i < k) ∨
      (j < k ∧ k < i) ∨ (k < i ∧ i < j) ∨ (k < j ∧ j < i) := by omega
  rcases hord with h | h | h | h | h | h
  · exact hthree i j k h.1 h.2 hi hj hk
  · exact hthree i k j h.1 h.2 hi hk hj
  · exact hthree j i k h.1 h.2 hj hi hk
  · exact hthree j k i h.1 h.2 hj hk hi
  · exact hthree k i j h.1 h.2 hk hi hj
  · exact hthree k j i h.1 h.2 hk hj hi

theorem shortest_odd_cycle_length_five_le {u : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p)
    (htriangle : ∀ a b c, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)) :
    5 ≤ p.length := by
  have hthree := (shortest_odd_walk_isCycle hp).three_le_length
  have ho := Nat.odd_iff.mp hp.1
  by_contra h
  have hl : p.length = 3 := by omega
  have ha := p.adj_getVert_succ (i := 0) (by omega)
  have hb := p.adj_getVert_succ (i := 1) (by omega)
  have hc := p.adj_getVert_succ (i := 2) (by omega)
  have he : p.getVert 3 = p.getVert 0 := by rw [← hl, Walk.getVert_length, Walk.getVert_zero]
  exact htriangle _ _ _ ⟨ha, hb, by simpa [he] using hc⟩

theorem shortest_odd_cycle_injective {u : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p) :
    Function.Injective (fun i : Fin p.length => p.getVert i.val) := by
  intro i j hij
  apply Fin.ext
  exact (shortest_odd_walk_isCycle hp).getVert_injOn'
    (by simp only [Set.mem_ofPred_eq]; omega)
    (by simp only [Set.mem_ofPred_eq]; omega) hij

theorem five_cycle_copy_of_shortest_walk {u : V} {p : G.Walk u u}
    (hp : IsShortestOddWalk p) (hl : p.length = 5) :
    ∃ φ, IsCycleCopy G φ := by
  let φ : CycleVertex → V := fun i => p.getVert i.val
  have hinj : Function.Injective φ := by
    intro i j hij
    have hh := shortest_odd_cycle_injective hp
      (a₁ := ⟨i.val, by omega⟩) (a₂ := ⟨j.val, by omega⟩) hij
    exact Fin.ext (congrArg (fun x : Fin p.length => x.val) hh)
  have hadj : ∀ i : CycleVertex, G.Adj (φ i) (φ (i+1)) := by
    intro i
    have hh := p.adj_getVert_succ (i := i.val) (by omega)
    have he : p.getVert 5 = p.getVert 0 := by
      rw [← hl, Walk.getVert_length, Walk.getVert_zero]
    fin_cases i <;> simpa [φ, he] using hh
  refine ⟨φ, hinj, ?_⟩
  intro i j hij
  change cycleRel i j at hij
  rcases hij with hij | hji
  · have he : i+1 = j := Fin.ext (by change (i.val+1)%5 = j.val; exact hij)
    simpa [he] using hadj i
  · have he : j+1 = i := Fin.ext (by change (j.val+1)%5 = i.val; exact hji)
    simpa [he] using (hadj j).symm

end TriangleProof
