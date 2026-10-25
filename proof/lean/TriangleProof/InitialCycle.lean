import TriangleProof.Prerequisites
import TriangleProof.ShortestOddCycle
import Mathlib.Combinatorics.SimpleGraph.Bipartite

/-!
# The initial-cycle step, with both standard graph facts proved

The deficiency calculation and bipartite edge bound are proved here. The
shortest odd cycle is supplied by `ShortestOddCycle.lean`.
It says either that cycle is a `C5`, or that it has at least seven vertices
and every graph vertex has at most two neighbors on it.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

def LongCycleIncidenceAlternative : Prop :=
  (∃ φ, IsCycleCopy G φ) ∨
    ∃ (l : ℕ) (φ : Fin l → V),
      7 ≤ l ∧ Function.Injective φ ∧
        ∀ v, (Finset.univ.filter fun i => G.Adj v (φ i)).card ≤ 2

/-- The degree-sum bound on the smaller side of a bipartition. -/
def BipartiteComponentBound : Prop :=
  ∀ (d t : ℤ),
    0 ≤ d → (Fintype.card V : ℤ) = 2 * d + 2 * t + 1 →
    (∀ v, degreeZ G v ≤ d) → G.IsBipartite →
    edgeCountZ G ≤ d * (d + t)

omit [DecidableEq V] in
/-- The bipartite component bound is already available from Mathlib's
bipartition degree-sum theorem. -/
theorem bipartite_component_bound : BipartiteComponentBound G := by
  intro d t hd hn hdeg hbi
  obtain ⟨s, u, hsu⟩ := hbi.exists_isBipartiteWith
  let _ : Fintype s := Fintype.ofFinite s
  let _ : Fintype u := Fintype.ofFinite u
  have hcardsN : s.ncard + u.ncard ≤ Fintype.card V := by
    rw [← Set.ncard_union_eq hsu.disjoint]
    simpa using Set.ncard_mono (Set.subset_univ (s ∪ u))
  have hcards : (s.ncard : ℤ) + (u.ncard : ℤ) ≤ (Fintype.card V : ℤ) := by
    exact_mod_cast hcardsN
  have hsmall : (s.ncard : ℤ) ≤ d + t ∨ (u.ncard : ℤ) ≤ d + t := by
    omega
  rcases hsmall with hs | hu
  · have hedge : edgeCountZ G = ∑ v ∈ s.toFinset, degreeZ G v := by
      simp_rw [degreeZ_eq_degree]
      rw [← Nat.cast_sum,
        SimpleGraph.isBipartiteWith_sum_degrees_eq_card_edges
          (s := s.toFinset) (t := u.toFinset) (by simpa using hsu)]
      simp [edgeCountZ]
    have hsum := Finset.sum_le_sum (fun v (_ : v ∈ s.toFinset) => hdeg v)
    simp only [Finset.sum_const, nsmul_eq_mul] at hsum
    rw [hedge]
    have hmul : (s.ncard : ℤ) * d ≤ (d + t) * d := by
      have := mul_nonneg (sub_nonneg.mpr hs) hd
      nlinarith
    have hmul' : (s.toFinset.card : ℤ) * d ≤ d * (d + t) := by
      simpa [Set.ncard_eq_toFinset_card, mul_comm] using hmul
    exact hsum.trans hmul'
  · have hedge : edgeCountZ G = ∑ v ∈ u.toFinset, degreeZ G v := by
      simp_rw [degreeZ_eq_degree]
      rw [← Nat.cast_sum,
        SimpleGraph.isBipartiteWith_sum_degrees_eq_card_edges'
          (s := s.toFinset) (t := u.toFinset) (by simpa using hsu)]
      simp [edgeCountZ]
    have hsum := Finset.sum_le_sum (fun v (_ : v ∈ u.toFinset) => hdeg v)
    simp only [Finset.sum_const, nsmul_eq_mul] at hsum
    rw [hedge]
    have hmul : (u.ncard : ℤ) * d ≤ (d + t) * d := by
      have := mul_nonneg (sub_nonneg.mpr hu) hd
      nlinarith
    have hmul' : (u.toFinset.card : ℤ) * d ≤ d * (d + t) := by
      simpa [Set.ncard_eq_toFinset_card, mul_comm] using hmul
    exact hsum.trans hmul'

/-- Standard shortest-odd-cycle output.  Triangle-freeness rules out length
three; minimality rules out three neighbors on a cycle of length at least
seven. -/
def ShortestOddCycleAlternative : Prop :=
  (∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)) →
  ¬ G.IsBipartite → LongCycleIncidenceAlternative G

omit [Fintype V] [DecidableEq V] in
/-- The shortest-cycle alternative is a theorem, not an external input. -/
theorem shortest_odd_cycle_alternative : ShortestOddCycleAlternative G := by
  intro htriangle hnonbip
  obtain ⟨u, p, _hcycle, hp⟩ := exists_shortest_odd_cycle hnonbip
  have hl := shortest_odd_cycle_length_five_le hp htriangle
  by_cases hfive : p.length = 5
  · exact Or.inl (five_cycle_copy_of_shortest_walk hp hfive)
  · have ho := Nat.odd_iff.mp hp.1
    exact Or.inr ⟨p.length, fun i => p.getVert i.val, by omega,
      shortest_odd_cycle_injective hp, shortest_odd_cycle_neighbor_bound hp hl⟩

/-- All numerical and incidence counting in the initial-cycle argument. -/
theorem dense_triangle_free_five_cycle_of_standard_facts
    (hbip : BipartiteComponentBound G)
    (hodd : ShortestOddCycleAlternative G) :
    DenseTriangleFreeFiveCycle G := by
  intro d t hd ht htd hn htriangle hdeg he
  have hnonbip : ¬ G.IsBipartite := by
    intro hbi
    have hub := hbip d t (by omega) hn hdeg hbi
    nlinarith
  rcases hodd htriangle hnonbip with hfive | ⟨l, φ, hl, hφ, hinc⟩
  · exact hfive
  exfalso
  let N := fun v => Finset.univ.filter (fun i : Fin l => G.Adj v (φ i))
  have hcount : (∑ i, degreeZ G (φ i)) = ∑ v, ((N v).card : ℤ) := by
    rw [sum_degreeZ_on_copy]
    apply Finset.sum_congr rfl
    intro v _
    simp [N, incidence]
  have hupper : (∑ v, ((N v).card : ℤ)) ≤ 2 * (Fintype.card V : ℤ) := by
    have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) => hinc v)
    norm_num only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hh
    have hh' : (∑ v, ((N v).card : ℤ)) ≤ (Fintype.card V : ℤ) * 2 := by
      exact_mod_cast hh
    simpa [mul_comm] using hh'
  let δ := fun v => d - degreeZ G v
  have hδ : ∀ v, 0 ≤ δ v := fun v => sub_nonneg.mpr (hdeg v)
  have hcopy_le : (∑ i, δ (φ i)) ≤ deficiency G d := by
    exact sum_copy_le_sum φ hφ δ hδ
  have hcopy : (∑ i, δ (φ i)) = (l : ℤ) * d - ∑ i, degreeZ G (φ i) := by
    simp only [δ, Finset.sum_sub_distrib, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hDform : deficiency G d =
      d * (2 * d + 2 * t + 1) - 2 * edgeCountZ G := by
    rw [deficiency_eq G, hn]
  have hD : deficiency G d ≤ d - 2 * t - 4 :=
    counterexample_deficiency d t (edgeCountZ G) (deficiency G d) hDform he
  rw [hcount] at hcopy
  have hlZ : (7 : ℤ) ≤ (l : ℤ) := by exact_mod_cast hl
  nlinarith

/-- The graph-level initial-cycle statement, with no external premises. -/
theorem dense_triangle_free_five_cycle : DenseTriangleFreeFiveCycle G :=
  dense_triangle_free_five_cycle_of_standard_facts G
    (bipartite_component_bound G) (shortest_odd_cycle_alternative G)

/-- The initial-cycle principle has no external graph-theoretic premises. -/
theorem initial_cycle_principle : InitialCyclePrinciple := by
  intro n G
  let _ : DecidableRel G.Adj := Classical.decRel G.Adj
  exact dense_triangle_free_five_cycle G

end TriangleProof
