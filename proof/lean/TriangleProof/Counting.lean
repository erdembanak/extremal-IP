import Mathlib.Tactic

/-!
# General arithmetic implications used in the proof draft

These lemmas are quantified over arbitrary integers. Later graph modules
establish the incidence and deficiency hypotheses needed by the main proof.
-/

namespace TriangleProof

/-- The pigeonhole step used to find a saturated vertex in each cycle bag. -/
theorem exists_zero_deficiency {V : Type*} (s : Finset V) (δ : V → ℕ)
    (hbudget : (∑ v ∈ s, δ v) < s.card) : ∃ v ∈ s, δ v = 0 := by
  by_contra h
  push Not at h
  have hsum : (∑ _v ∈ s, (1 : ℕ)) ≤ ∑ v ∈ s, δ v := by
    apply Finset.sum_le_sum
    intro v hv
    have := h v hv
    omega
  simp only [Finset.sum_const, smul_eq_mul, mul_one] at hsum
  omega

theorem counterexample_deficiency (d t e D : ℤ)
    (hD : D = d * (2 * d + 2 * t + 1) - 2 * e)
    (he : d * (d + t) + t + 1 < e) :
    D ≤ d - 2 * t - 4 := by
  have he' : d * (d + t) + t + 2 ≤ e := by omega
  nlinarith

theorem enough_saturated_representatives (t r q D_out : ℤ)
    (hr : r ≤ q) (hbudget : q + D_out ≤ 2 * t - 2) :
    D_out < 2 * t - r := by
  omega

theorem saturated_cycle_contradiction (d t : ℤ)
    (hd : 4 * t + 3 ≤ d)
    (hcycle : 5 * d ≤ 2 * (2 * d + 2 * t + 1)) : False := by
  omega

theorem mobius_bag_lower_bound (d n t D a δ₁ δ₂ δ₃ : ℤ)
    (hd : d = 4 * t - 1) (hn : n = 10 * t - 1)
    (hD : D ≤ 2 * t - 5) (hδ : δ₁ + δ₂ + δ₃ ≤ D)
    (ha : a = 3 * d - n - (δ₁ + δ₂ + δ₃)) : 3 ≤ a := by
  omega

theorem small_obstruction_contradiction (d n t D D_S degreeSum : ℤ)
    (hd : d = 4 * t - 1) (hn : n = 10 * t - 1)
    (hD : D ≤ 2 * t - 5) (hsub : D_S ≤ D)
    (hidentity : D_S = 8 * d - degreeSum)
    (hcount : degreeSum ≤ 3 * n - 2) : False := by
  omega

theorem one_exception_positive (h a u v : ℤ)
    (ha : 0 ≤ a) (hu : 1 ≤ u) (hv : 1 ≤ v) :
    h - 1 ≤ h + 1 + 2 * a - u - v + u * (v - 1) + v * (u - 1) := by
  have hprod : 0 ≤ (u - 1) * (v - 1) := mul_nonneg (by omega) (by omega)
  nlinarith

theorem packing_pair_large (d z b a c : ℤ)
    (hz : z ≤ d + b - 1) :
    (a - d + 1) + (c - d + 1) ≤ b + (a + c - z - d + 1) := by
  omega

theorem packing_pair_small (d z b a c : ℤ)
    (hz : z ≤ d + b - 1) (hsmall : a + c - z < d) :
    (a - d + 1) + (c - d + 1) ≤ b := by
  omega

end TriangleProof
