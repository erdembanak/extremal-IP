import TriangleProof.AuditedLiterature

/-!
# The odd-degree threshold

For odd `d ≥ 9`, the fixed-order bound excludes every threshold candidate
with matching number strictly between `d` and `⌊5(d+1)/4⌋`. With the published
lower bound `d+1 ≤ Z(d)` this gives `⌊5(d+1)/4⌋ ≤ Z(d)`; with the published
construction it identifies `Z(d)` exactly. The degree-seven case is
`Z_seven_eq_nine`.
-/

namespace TriangleProof

/-- Below the conjectured threshold, the shift `t = i - d` lies in one of the
three numerical regimes of the fixed-order theorem. -/
theorem odd_belowThreshold_regime (d i : ℕ) (hd : 9 ≤ d) (hodd : Odd d)
    (hdi : d < i) (hi : i < 5 * (d + 1) / 4) :
    FixedOrderRegime (d : ℤ) ((i : ℤ) - (d : ℤ)) := by
  obtain ⟨k, rfl⟩ := hodd
  by_cases hlarge : 4 * ((i : ℤ) - ((2 * k + 1 : ℕ) : ℤ)) + 3 ≤ ((2 * k + 1 : ℕ) : ℤ)
  · exact Or.inl hlarge
  have hmod : k % 2 = 0 ∨ k % 2 = 1 := by omega
  rcases hmod with h0 | h1
  · right; right
    refine ⟨(i : ℤ) - ((2 * k + 1 : ℕ) : ℤ), ?_, ?_, rfl⟩ <;> push_cast at hlarge ⊢ <;> omega
  · right; left
    refine ⟨(i : ℤ) - ((2 * k + 1 : ℕ) : ℤ), ?_, ?_, rfl⟩ <;> push_cast at hlarge ⊢ <;> omega

/-- The fixed-order bound on `2i+1` vertices excludes a threshold candidate
below `⌊5(d+1)/4⌋`, with no external hypotheses. -/
theorem odd_not_candidate_below (d i : ℕ) (hd : 9 ≤ d) (hodd : Odd d)
    (hdi : d < i) (hi : i < 5 * (d + 1) / 4) :
    ¬ IsZCandidate d i := by
  rintro ⟨G, htri, _hcrit, _hmatch, hprofile⟩
  let _ : DecidableRel G.Adj := Classical.decRel G.Adj
  rcases hprofile with ⟨heven, _⟩ | ⟨_hodd, hmax, hsum⟩
  · exact (Nat.not_even_iff_odd.mpr hodd) heven
  have hdeg : ∀ v, degreeZ G v ≤ (d : ℤ) := by
    intro v
    rw [degreeZ_eq_degree, ← degreeN_eq_degree]
    exact_mod_cast hmax v
  have hD : deficiency G (d : ℤ) = 1 := by
    calc
      deficiency G (d : ℤ) = ∑ v, ((d - degreeN G v : ℕ) : ℤ) := by
        apply Finset.sum_congr rfl
        intro v _
        rw [Nat.cast_sub (hmax v), degreeN_eq_degree, degreeZ_eq_degree]
      _ = 1 := by rw [← Nat.cast_sum, hsum]; norm_num
  have hn : (Fintype.card (Fin (2 * i + 1)) : ℤ) =
      2 * (d : ℤ) + 2 * ((i : ℤ) - (d : ℤ)) + 1 := by
    rw [Fintype.card_fin]; push_cast; ring
  have he := fixed_order_edge_bound G htri (d : ℤ) ((i : ℤ) - (d : ℤ))
    (by exact_mod_cast (by omega : 7 ≤ d)) (by omega) (by omega) hn hdeg
    (odd_belowThreshold_regime d i hd hodd hdi hi)
  have hid := deficiency_eq G (d : ℤ)
  rw [hn, hD] at hid
  -- `d ≤ 2t+3`, while `4t ≤ d+1` and `d ≥ 9`.
  have hsmall : (d : ℤ) ≤ 2 * ((i : ℤ) - (d : ℤ)) + 3 := by nlinarith [hid, he]
  omega

/-- The lower bound uses only existence of a candidate and the published
general lower bound `d+1 ≤ Z(d)`. -/
theorem odd_Z_lower (d : ℕ) (hd : 9 ≤ d) (hodd : Odd d)
    (hexists : ∃ i, IsZCandidate d i) (hlower : d + 1 ≤ Z d) :
    5 * (d + 1) / 4 ≤ Z d := by
  by_contra h
  exact odd_not_candidate_below d (Z d) hd hodd (by omega) (by omega)
    (Z_isZCandidate d hexists)

/-- The paper's odd threshold value, from the published Z results. -/
theorem odd_Z_eq_of_publishedZResults (hz : PublishedZResults)
    (d : ℕ) (hd : 9 ≤ d) (hodd : Odd d) :
    Z d = 5 * (d + 1) / 4 := by
  have hlow := odd_Z_lower d hd hodd (hz.candidateExists d (by omega))
    (hz.lower d (by omega))
  have hup := hz.oddUpper d hodd
  omega

end TriangleProof
