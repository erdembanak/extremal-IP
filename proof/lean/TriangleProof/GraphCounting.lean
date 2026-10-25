import TriangleProof.Obstructions
import Mathlib.Combinatorics.SimpleGraph.DegreeSum

/-! Integer-valued incidence counting for arbitrary finite simple graphs.
The degree used here is proved equal to Mathlib's graph degree. -/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

def incidence (v w : V) : ℤ := if G.Adj v w then 1 else 0

def degreeZ (v : V) : ℤ := ∑ w, incidence G v w

def deficiency (d : ℤ) : ℤ := ∑ v, (d - degreeZ G v)

def edgeCountZ : ℤ := (G.edgeFinset.card : ℤ)

omit [Fintype V] [DecidableEq V] in
theorem incidence_symm (v w : V) : incidence G v w = incidence G w v := by
  simp only [incidence, G.adj_comm]

omit [Fintype V] [DecidableEq V] in
theorem incidence_nonneg (v w : V) : 0 ≤ incidence G v w := by
  unfold incidence
  split_ifs <;> norm_num

omit [DecidableEq V] in
theorem degreeZ_eq_degree (v : V) : degreeZ G v = (G.degree v : ℤ) := by
  classical
  rw [← G.card_neighborFinset_eq_degree]
  simp [degreeZ, incidence, SimpleGraph.neighborFinset_def]

omit [DecidableEq V] in
theorem sum_degreeZ_eq_twice_edgeCountZ :
    (∑ v, degreeZ G v) = 2 * edgeCountZ G := by
  simp_rw [degreeZ_eq_degree]
  rw [← Nat.cast_sum, G.sum_degrees_eq_twice_card_edges]
  simp [edgeCountZ]

omit [DecidableEq V] in
theorem deficiency_eq (d : ℤ) :
    deficiency G d = d * (Fintype.card V : ℤ) - 2 * edgeCountZ G := by
  simp only [deficiency, Finset.sum_sub_distrib, Finset.sum_const,
    Finset.card_univ, nsmul_eq_mul, sum_degreeZ_eq_twice_edgeCountZ]
  ring

omit [DecidableEq V] in
theorem sum_degreeZ_on_copy {k : ℕ} (φ : Fin k → V) :
    (∑ i, degreeZ G (φ i)) = ∑ v, ∑ i, incidence G v (φ i) := by
  simp only [degreeZ]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro v _
  apply Finset.sum_congr rfl
  intro i _
  exact incidence_symm G _ _

theorem sum_copy_le_sum {k : ℕ} (φ : Fin k → V) (hφ : Function.Injective φ)
    (f : V → ℤ) (hf : ∀ v, 0 ≤ f v) :
    (∑ i, f (φ i)) ≤ ∑ v, f v := by
  rw [← Finset.sum_image (fun i _ j _ h => hφ h)]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) (by
    intro v _ _
    exact hf v)

theorem sum_copy_add_outside {k : ℕ} (φ : Fin k → V) (hφ : Function.Injective φ)
    (f : V → ℤ) :
    (∑ i, f (φ i)) + (∑ v ∈ Finset.univ.filter (fun v => v ∉ Set.range φ), f v)
      = ∑ v, f v := by
  classical
  rw [← Finset.sum_image (fun i _ j _ h => hφ h)]
  have hset : Finset.univ.filter (fun v => v ∉ Set.range φ) =
      Finset.univ \ Finset.univ.image φ := by ext v; simp
  rw [hset]
  exact Finset.sum_add_sum_compl (Finset.univ.image φ) f

omit [DecidableEq V] in
theorem zero_of_sum_eq_zero (f : V → ℤ) (hf : ∀ v, 0 ≤ f v)
    (hs : ∑ v, f v = 0) (v : V) : f v = 0 := by
  exact (Finset.sum_eq_zero_iff_of_nonneg (fun w _ => hf w)).mp hs v (Finset.mem_univ _)

omit [DecidableEq V] in
theorem fiber_sum {k : ℕ} (p : V → Fin k) (f : Fin k → ℤ) :
    (∑ v, f (p v)) = ∑ j, ((Finset.univ.filter (fun v => p v = j)).card : ℤ) * f j := by
  classical
  calc
    _ = ∑ j, ∑ v ∈ Finset.univ.filter (fun v => p v = j), f (p v) := by
      symm
      simpa using Finset.sum_fiberwise_eq_sum_filter Finset.univ Finset.univ p (fun v => f (p v))
    _ = _ := by
      apply Finset.sum_congr rfl
      intro j _
      calc
        _ = ∑ _v ∈ Finset.univ.filter (fun v => p v = j), f j := by
          apply Finset.sum_congr rfl
          intro v hv
          rw [(Finset.mem_filter.mp hv).2]
        _ = _ := by simp

omit [DecidableEq V] in
/-- An independent vertex set receives at most half of the total incidence
weight of any symmetric nonnegative matrix supported on graph edges. -/
theorem independent_weight_bound (w : V → V → ℤ)
    (hsymm : ∀ v u, w v u = w u v) (hnonneg : ∀ v u, 0 ≤ w v u)
    (s : V → Prop) [DecidablePred s]
    (hind : ∀ v u, s v → s u → w v u = 0) :
    2 * (∑ v, if s v then ∑ u, w v u else 0) ≤ ∑ v, ∑ u, w v u := by
  have hp : ∀ v u, (if s v then w v u else 0) + (if s u then w v u else 0) ≤ w v u := by
    intro v u
    by_cases hv : s v <;> by_cases hu : s u <;> simp [hv, hu, hind, hnonneg]
  have hh := Finset.sum_le_sum (fun v (_ : v ∈ Finset.univ) =>
    Finset.sum_le_sum (fun u (_ : u ∈ Finset.univ) => hp v u))
  simp only [Finset.sum_add_distrib] at hh
  have hleft : (∑ v, ∑ u, if s v then w v u else 0) =
      ∑ v, if s v then ∑ u, w v u else 0 := by
    apply Finset.sum_congr rfl
    intro v _
    by_cases hv : s v <;> simp [hv]
  have hright : (∑ v, ∑ u, if s u then w v u else 0) =
      ∑ v, if s v then ∑ u, w v u else 0 := by
    rw [Finset.sum_comm]
    simpa only [hsymm] using hleft
  rw [hleft, hright] at hh
  omega

end TriangleProof
