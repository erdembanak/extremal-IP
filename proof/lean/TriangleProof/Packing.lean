import TriangleProof.Counting

/-!
# Arithmetic layer of the component packing argument

The graph decomposition is a published input.  These lemmas record the
utility calculation used after that decomposition and do not assume either
conjecture.
-/

namespace TriangleProof

/-- Excess over the baseline `d*j` contributed by a non-star component of
matching weight `j`. -/
def componentUtility (d z b j : ℤ) : ℤ :=
  if j = z then b else j - d + 1

theorem componentUtility_at_threshold (d z b : ℤ) :
    componentUtility d z b z = b := by
  simp [componentUtility]

theorem componentUtility_below_threshold (d z b j : ℤ) (hj : j < z) :
    componentUtility d z b j = j - d + 1 := by
  simp [componentUtility, ne_of_lt hj]

/-- Replacing two sub-threshold components by a threshold component and a
second legal component does not lower utility. -/
theorem packing_replace_two_large (d z b a c : ℤ)
    (ha : a < z) (hc : c < z) (_hremainder : d ≤ a + c - z)
    (hz : z ≤ d + b - 1) :
    componentUtility d z b a + componentUtility d z b c ≤
      componentUtility d z b z + componentUtility d z b (a + c - z) := by
  rw [componentUtility_below_threshold d z b a ha,
    componentUtility_below_threshold d z b c hc,
    componentUtility_at_threshold]
  have hrem_lt : a + c - z < z := by omega
  rw [componentUtility_below_threshold d z b _ hrem_lt]
  exact packing_pair_large d z b a c hz

/-- If the leftover weight is too small to support another component, the
threshold component alone has at least the utility of the old pair. -/
theorem packing_replace_two_small (d z b a c : ℤ)
    (ha : a < z) (hc : c < z) (hremainder : a + c - z < d)
    (hz : z ≤ d + b - 1) :
    componentUtility d z b a + componentUtility d z b c ≤
      componentUtility d z b z := by
  rw [componentUtility_below_threshold d z b a ha,
    componentUtility_below_threshold d z b c hc,
    componentUtility_at_threshold]
  exact packing_pair_small d z b a c hz hremainder

/-- Utility of the one possible residual component after pair replacement. -/
def residualUtility (d r : ℤ) : ℤ := if r < d then 0 else r - d + 1

theorem residual_component_bound (d z b r j : ℤ)
    (hj0 : d ≤ j) (hjr : j ≤ r) (hrz : r < z) :
    componentUtility d z b j ≤ residualUtility d r := by
  have hjz : j < z := lt_of_le_of_lt hjr hrz
  rw [componentUtility_below_threshold d z b j hjz]
  simp only [residualUtility]
  split_ifs with hrd
  · omega
  · omega

/-- The utility formula once pair replacement has produced exactly `k`
threshold components and at most one residual component. -/
theorem normalized_packing_bound (d b k r q : ℤ)
    (hq : q ≤ residualUtility d r) :
    k * b + q ≤ k * b + residualUtility d r := by
  omega

theorem residualUtility_eq (d r : ℤ) :
    residualUtility d r = if r < d then 0 else r - d + 1 := rfl

/-- Repeated pair replacement terminates with threshold components and at
most one sub-threshold component.  `s = 0` denotes that there is no residual
component. -/
theorem normalize_small_components (d z b : ℤ)
    (_hd : 1 ≤ d) (hz2 : z < 2 * d) (hz : z ≤ d + b - 1)
    (xs : List ℤ) (hvalid : ∀ j ∈ xs, d ≤ j ∧ j < z) :
    ∃ q : ℕ, ∃ s : ℤ,
      (s = 0 ∨ (d ≤ s ∧ s < z)) ∧
      (q : ℤ) * z + s ≤ xs.sum ∧
      (xs.map fun j => j - d + 1).sum ≤
        (q : ℤ) * b + (if s = 0 then 0 else s - d + 1) := by
  induction xs using List.twoStepInduction with
  | nil =>
      exact ⟨0, 0, Or.inl rfl, by simp, by simp⟩
  | singleton a =>
      have ha := hvalid a (by simp)
      refine ⟨0, a, Or.inr ha, ?_, ?_⟩
      · simp
      · have hane : a ≠ 0 := by omega
        simp [hane]
  | cons_cons a c ys ih0 ih1 =>
      have ha := hvalid a (by simp)
      have hc := hvalid c (by simp)
      have hys : ∀ j ∈ ys, d ≤ j ∧ j < z := by
        intro j hj
        exact hvalid j (by simp [hj])
      by_cases hrem : d ≤ a + c - z
      · have hremz : a + c - z < z := by omega
        have hvalid' : ∀ j ∈ (a + c - z) :: ys, d ≤ j ∧ j < z := by
          intro j hj
          simp only [List.mem_cons] at hj
          rcases hj with rfl | hj
          · exact ⟨hrem, hremz⟩
          · exact hys j hj
        obtain ⟨q, s, hs, hw, hu⟩ := ih1 (a + c - z) hvalid'
        refine ⟨q + 1, s, hs, ?_, ?_⟩
        · push_cast
          simp only [List.sum_cons] at hw ⊢
          nlinarith
        · have hp := packing_pair_large d z b a c hz
          push_cast
          simp only [List.map_cons, List.sum_cons] at hu ⊢
          nlinarith
      · have hsmall : a + c - z < d := by omega
        obtain ⟨q, s, hs, hw, hu⟩ := ih0 hys
        refine ⟨q + 1, s, hs, ?_, ?_⟩
        · push_cast
          simp only [List.sum_cons] at hw ⊢
          nlinarith
        · have hp := packing_pair_small d z b a c hz hsmall
          push_cast
          simp only [List.map_cons, List.sum_cons] at hu ⊢
          nlinarith

/-- Full numerical packing bound.  `q₀` counts components of weight `z` and
`xs` lists all smaller non-star components.  Stars consume unused matching
capacity and have zero utility. -/
theorem packing_utility_bound (d z b k r : ℤ) (q₀ : ℕ) (xs : List ℤ)
    (hd : 1 ≤ d) (hzpos : 1 ≤ z) (hb : 0 ≤ b)
    (hz2 : z < 2 * d) (hz : z ≤ d + b - 1) (_hr0 : 0 ≤ r) (hrz : r < z)
    (hvalid : ∀ j ∈ xs, d ≤ j ∧ j < z)
    (hweight : (q₀ : ℤ) * z + xs.sum ≤ k * z + r) :
    (q₀ : ℤ) * b + (xs.map fun j => j - d + 1).sum ≤
      k * b + residualUtility d r := by
  obtain ⟨q₁, s, hs, hw, hu⟩ :=
    normalize_small_components d z b hd hz2 hz xs hvalid
  let Q : ℕ := q₀ + q₁
  have hs0 : 0 ≤ s := by
    rcases hs with rfl | hs
    · norm_num
    · omega
  have hQweight : (Q : ℤ) * z + s ≤ k * z + r := by
    dsimp [Q]
    nlinarith
  have hQk : (Q : ℤ) ≤ k := by
    by_contra h
    have hkQ : k < (Q : ℤ) := lt_of_not_ge h
    nlinarith
  have hutil : (q₀ : ℤ) * b + (xs.map fun j => j - d + 1).sum ≤
      (Q : ℤ) * b + (if s = 0 then 0 else s - d + 1) := by
    dsimp [Q]
    nlinarith
  by_cases hQeq : (Q : ℤ) = k
  · have hsr : s ≤ r := by nlinarith
    have hsres : (if s = 0 then 0 else s - d + 1) ≤ residualUtility d r := by
      rcases hs with hs | hs
      · subst s
        simp only
        simp only [residualUtility]
        split_ifs <;> omega
      · have hcomp := residual_component_bound d z b r s hs.1 hsr hrz
        rw [componentUtility_below_threshold d z b s hs.2] at hcomp
        simpa [ne_of_gt (lt_of_lt_of_le (by omega : 0 < d) hs.1)] using hcomp
    nlinarith
  · have hQlt : (Q : ℤ) < k := by
      omega
    have hsutil : (if s = 0 then 0 else s - d + 1) ≤ b := by
      rcases hs with hs | hs
      · simp [hs, hb]
      · simp only [if_neg (ne_of_gt (lt_of_lt_of_le (by omega : 0 < d) hs.1))]
        omega
    have hkb : (Q : ℤ) * b + b ≤ k * b := by
      have hnonneg : 0 ≤ (k - (Q : ℤ) - 1) * b :=
        mul_nonneg (by omega) hb
      nlinarith
    have hrnonneg : 0 ≤ residualUtility d r := by
      simp only [residualUtility]
      split_ifs <;> omega
    nlinarith

/-- Once the published component reduction supplies a packing certificate,
the numerical normalization determines the extremal value exactly. -/
theorem extremal_value_from_packing (d z b k r F : ℤ)
    (hd : 1 ≤ d) (hzpos : 1 ≤ z) (hb : 0 ≤ b)
    (hz2 : z < 2 * d) (hz : z ≤ d + b - 1)
    (hr0 : 0 ≤ r) (hrz : r < z)
    (hlower : d * (k * z + r) + k * b + residualUtility d r ≤ F)
    (hcertificate : ∃ (q₀ : ℕ) (xs : List ℤ),
      (∀ j ∈ xs, d ≤ j ∧ j < z) ∧
      (q₀ : ℤ) * z + xs.sum ≤ k * z + r ∧
      F ≤ d * (k * z + r) +
        ((q₀ : ℤ) * b + (xs.map fun j => j - d + 1).sum)) :
    F = d * (k * z + r) + k * b + residualUtility d r := by
  obtain ⟨q₀, xs, hvalid, hweight, hupper⟩ := hcertificate
  have hpacking := packing_utility_bound d z b k r q₀ xs
    hd hzpos hb hz2 hz hr0 hrz hvalid hweight
  apply le_antisymm
  · nlinarith
  · exact hlower

end TriangleProof
