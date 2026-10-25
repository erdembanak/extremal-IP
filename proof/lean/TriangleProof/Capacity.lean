import Mathlib.Tactic

/-!
# The capacity calculation from Section 6.4

The two pattern theorems cover every admissible five-tuple for excess three
and one. They are retained as independent regression checks. The production
proof in `CycleCapacity.lean` now uses the general weighted-bipartition lemma
formalizing the revised manuscript argument.
-/

namespace TriangleProof

/-- The scalar inequality used after bipartitioning the positive support of
an odd excess. `CycleCapacity.lean` proves the graph-to-weight reduction. -/
theorem odd_excess_imbalance_scalar (h k P : ℤ)
    (_hk : 1 ≤ k) (hP : k ≤ P) (hPq : P ≤ 2 * k - 1)
    (hh : 2 * k - 1 ≤ h) :
    h - k * (k - 1) ≤ h * (2 * P - (2 * k - 1)) - P * (P - 1) := by
  have hfactor : 0 ≤ 2 * h - P - k + 1 := by omega
  have hprod := mul_nonneg (sub_nonneg.mpr hP) hfactor
  nlinarith

set_option maxHeartbeats 0 in
theorem capacity_pattern_bound (t : ℤ) (ht : 3 ≤ t)
    (p0 p1 p2 p3 p4 : ℕ) (hs : p0 + p1 + p2 + p3 + p4 = 3) :
    let w0 := (2 * t + 1 - (p2 : ℤ) - (p3 : ℤ)) * (p0 : ℤ)
    let w1 := (2 * t + 1 - (p3 : ℤ) - (p4 : ℤ)) * (p1 : ℤ)
    let w2 := (2 * t + 1 - (p4 : ℤ) - (p0 : ℤ)) * (p2 : ℤ)
    let w3 := (2 * t + 1 - (p0 : ℤ) - (p1 : ℤ)) * (p3 : ℤ)
    let w4 := (2 * t + 1 - (p1 : ℤ) - (p2 : ℤ)) * (p4 : ℤ)
    let W := w0 + w1 + w2 + w3 + w4
    2 * t - 1 ≤ 2 * (w0 + w2) - W ∨
    2 * t - 1 ≤ 2 * (w0 + w3) - W ∨
    2 * t - 1 ≤ 2 * (w1 + w3) - W ∨
    2 * t - 1 ≤ 2 * (w1 + w4) - W ∨
    2 * t - 1 ≤ 2 * (w2 + w4) - W := by
  have h0 : p0 ≤ 3 := by omega
  have h1 : p1 ≤ 3 := by omega
  have h2 : p2 ≤ 3 := by omega
  have h3 : p3 ≤ 3 := by omega
  have h4 : p4 ≤ 3 := by omega
  dsimp only
  interval_cases p0 <;> interval_cases p1 <;> interval_cases p2 <;>
    interval_cases p3 <;> interval_cases p4 <;> norm_num at hs <;> norm_num <;> omega

set_option maxHeartbeats 0 in
theorem capacity_pattern_bound_one (t : ℤ) (_ht : 1 ≤ t)
    (p0 p1 p2 p3 p4 : ℕ) (hs : p0 + p1 + p2 + p3 + p4 = 1) :
    let w0 := (2 * t + 1 - (p2 : ℤ) - (p3 : ℤ)) * (p0 : ℤ)
    let w1 := (2 * t + 1 - (p3 : ℤ) - (p4 : ℤ)) * (p1 : ℤ)
    let w2 := (2 * t + 1 - (p4 : ℤ) - (p0 : ℤ)) * (p2 : ℤ)
    let w3 := (2 * t + 1 - (p0 : ℤ) - (p1 : ℤ)) * (p3 : ℤ)
    let w4 := (2 * t + 1 - (p1 : ℤ) - (p2 : ℤ)) * (p4 : ℤ)
    let W := w0 + w1 + w2 + w3 + w4
    2 * t + 1 ≤ 2 * (w0 + w2) - W ∨
    2 * t + 1 ≤ 2 * (w0 + w3) - W ∨
    2 * t + 1 ≤ 2 * (w1 + w3) - W ∨
    2 * t + 1 ≤ 2 * (w1 + w4) - W ∨
    2 * t + 1 ≤ 2 * (w2 + w4) - W := by
  have h0 : p0 ≤ 1 := by omega
  have h1 : p1 ≤ 1 := by omega
  have h2 : p2 ≤ 1 := by omega
  have h3 : p3 ≤ 1 := by omega
  have h4 : p4 ≤ 1 := by omega
  dsimp only
  interval_cases p0 <;> interval_cases p1 <;> interval_cases p2 <;>
    interval_cases p3 <;> interval_cases p4 <;> norm_num at hs <;> norm_num <;> omega

end TriangleProof
