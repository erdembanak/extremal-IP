import TriangleProof.LiteratureInputs
import TriangleProof.SevenConstruction
import TriangleProof.IntermediateConstruction

/-!
# Literature interface audited against the cited theorem statements

This file replaces broad consequence-shaped assumptions by a component-level
interface matching Ahanjideh--Ekim--Yildiz (2022), Corollary 3.5 and Lemmas
4.5/4.9.  Endpoint fields correspond to Theorems 3.4 and 4.8.  The lower
construction for the open component range, Proposition 4.1 of
Banak--Ekim--Taskin (2023), is proved in `IntermediateConstruction`.
-/

namespace TriangleProof

/-- A non-star component furnished by the structural reduction.  Its graph is
retained as a witness, so the new fixed-order theorem can be applied directly. -/
noncomputable def ComponentWitness (d j : ℕ) (e : ℤ) : Prop :=
  ∃ (n : ℕ) (G : SimpleGraph (Fin n)),
    let _ : DecidableRel G.Adj := Classical.decRel G.Adj
    (Fintype.card (Fin n) : ℤ) = 2 * (j : ℤ) + 1 ∧
    (∀ a b c, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)) ∧
    (∀ v, degreeZ G v ≤ (d : ℤ)) ∧
    edgeCountZ G = e ∧ e ≤ (fDelta d j : ℤ)

/-- List form of the published component decomposition. `q₀` counts threshold
components; each pair in `cs` records the matching number and edge count of
one remaining factor-critical component. -/
noncomputable def StructuralDecomposition (d m : ℕ) : Prop :=
  let z : ℤ := (Z d : ℕ)
  let b : ℤ := (d / 2 : ℕ)
  ∃ (q₀ : ℕ) (cs : List (ℕ × ℤ)),
    (∀ p ∈ cs, d ≤ p.1 ∧ p.1 < Z d ∧ ComponentWitness d p.1 p.2) ∧
    (q₀ : ℤ) * z + (cs.map fun p => (p.1 : ℤ)).sum ≤ (m : ℤ) ∧
    (fDelta d m : ℤ) ≤ (d : ℤ) * (m : ℤ) + (q₀ : ℤ) * b +
      (cs.map fun p => p.2 - (d : ℤ) * (p.1 : ℤ)).sum

/-- Regrouping the component-level structural theorem gives the numerical
knapsack decomposition used by the final arithmetic layer. -/
theorem knapsackDecomposition_of_structural
    (d m : ℕ) (hthreshold : fDelta d (Z d) = d * Z d + d / 2)
    (hstruct : StructuralDecomposition d m) : KnapsackDecomposition d m := by
  dsimp [StructuralDecomposition] at hstruct
  obtain ⟨q₀, cs, hvalid, hweight, hupper⟩ := hstruct
  let xs : List ℤ := cs.map fun p => (p.1 : ℤ)
  refine ⟨q₀, xs, ?_, ?_, ?_⟩
  · intro j hj
    simp only [xs, List.mem_map] at hj
    obtain ⟨p, hp, rfl⟩ := hj
    have hh := hvalid p hp
    exact ⟨by exact_mod_cast hh.1, by exact_mod_cast hh.2.1⟩
  · simpa [xs] using hweight
  · have hthresholdZ : (fDelta d (Z d) : ℤ) -
        (d : ℤ) * (Z d : ℤ) = (d / 2 : ℕ) := by
      rw [hthreshold]
      push_cast
      ring
    have hcomponent : ∀ p ∈ cs,
        p.2 - (d : ℤ) * (p.1 : ℤ) ≤
          (fDelta d p.1 : ℤ) - (d : ℤ) * (p.1 : ℤ) := by
      intro p hp
      have hw := (hvalid p hp).2.2
      dsimp [ComponentWitness] at hw
      obtain ⟨n, G, hn, htri, hdeg, he, hle⟩ := hw
      omega
    have hsum := List.sum_le_sum hcomponent
    rw [hthresholdZ]
    dsimp [xs]
    have hfinal : (fDelta d m : ℤ) ≤ (d : ℤ) * (m : ℤ) +
        (q₀ : ℤ) * (d / 2 : ℕ) +
        (cs.map fun p => (fDelta d p.1 : ℤ) - (d : ℤ) * (p.1 : ℤ)).sum := by
      nlinarith
    have hmap :
        ((cs.map fun p => (p.1 : ℤ)).map
          (fun j : ℤ => (fDelta d j.toNat : ℤ) - (d : ℤ) * j)).sum =
        (cs.map fun p => (fDelta d p.1 : ℤ) - (d : ℤ) * (p.1 : ℤ)).sum := by
      rw [List.map_map]
      exact congrArg List.sum (List.map_congr_left (by
        intro p hp
        simp))
    rw [hmap]
    push_cast at hfinal
    simpa only [add_assoc] using hfinal

/-- The structural component list becomes the simple utility certificate once
the left endpoint and every open-range fixed-order component bound are known. -/
theorem packingCertificate_of_fixed_components
    (d m : ℕ) (hleft : fDelta d d = d * d + 1)
    (hfixed : ∀ j : ℕ, d < j → j < Z d → FixedOrderComponentUpper d j)
    (hstruct : StructuralDecomposition d m) : PackingCertificate d m := by
  dsimp [StructuralDecomposition] at hstruct
  obtain ⟨q₀, cs, hvalid, hweight, hupper⟩ := hstruct
  let xs : List ℤ := cs.map fun p => (p.1 : ℤ)
  refine ⟨q₀, xs, ?_, ?_, ?_⟩
  · intro j hj
    simp only [xs, List.mem_map] at hj
    obtain ⟨p, hp, rfl⟩ := hj
    have hh := hvalid p hp
    exact ⟨by exact_mod_cast hh.1, by exact_mod_cast hh.2.1⟩
  · simpa [xs] using hweight
  · have hcomponent : ∀ p ∈ cs,
        p.2 - (d : ℤ) * (p.1 : ℤ) ≤ (p.1 : ℤ) - (d : ℤ) + 1 := by
      intro p hp
      have hv := hvalid p hp
      have hw := hv.2.2
      dsimp [ComponentWitness] at hw
      obtain ⟨n, G, hn, htri, hdeg, he, hle⟩ := hw
      rcases p with ⟨j, e⟩
      simp only at hv he hle ⊢
      by_cases hjd : j = d
      · subst j
        rw [hleft] at hle
        push_cast at hle
        nlinarith
      · have hdj : d < j := by omega
        have hub := hfixed j hdj hv.2.1
        have heupper := hub n G hn htri hdeg
        rw [he] at heupper
        nlinarith
    have hsum := List.sum_le_sum hcomponent
    dsimp [xs]
    have hmap :
        ((cs.map fun p => (p.1 : ℤ)).map fun j => j - (d : ℤ) + 1).sum =
        (cs.map fun p => (p.1 : ℤ) - (d : ℤ) + 1).sum := by
      simp [List.map_map, Function.comp_def]
    rw [hmap]
    nlinarith

theorem openComponentValue_natCast (d i : ℕ) (hdi : d ≤ i) :
    ((d * i + i - d + 1 : ℕ) : ℤ) =
      (d : ℤ) * (i : ℤ) + (i : ℤ) - (d : ℤ) + 1 := by
  have hsub : d ≤ d * i + i := by nlinarith
  rw [Nat.cast_add, Nat.cast_sub hsub, Nat.cast_add, Nat.cast_mul]
  norm_num

/-- The exact information about `Z` used below. `candidateExists` records the
threshold construction in Proposition 4.3 of Ahanjideh--Ekim--Yildiz (2022).
The remaining three fields are their Lemmas 4.6 and 4.7. The degree-seven
upper bound is proved by `SevenConstruction`, not supplied as an input. -/
structure PublishedZResults : Prop where
  candidateExists : ∀ d : ℕ, 2 ≤ d → ∃ i, IsZCandidate d i
  lower : ∀ d : ℕ, 4 ≤ d → d + 1 ≤ Z d
  even : ∀ d : ℕ, Even d → Z d = 5 * d / 4
  oddUpper : ∀ d : ℕ, Odd d → Z d ≤ 5 * (d + 1) / 4

/-- The cited threshold construction makes the defining set for `Z` nonempty. -/
theorem zCandidateSet_nonempty_of_publishedZResults (hz : PublishedZResults)
    (d : ℕ) (hd : 2 ≤ d) :
    {i | IsZCandidate d i}.Nonempty :=
  zCandidateSet_nonempty d (hz.candidateExists d hd)

/-- Hence the `sInf` in `Z` is an actual least threshold candidate. -/
theorem Z_isZCandidate_of_publishedZResults (hz : PublishedZResults)
    (d : ℕ) (hd : 2 ≤ d) : IsZCandidate d (Z d) :=
  Z_isZCandidate d (hz.candidateExists d hd)

/-- The three fixed-order regimes follow arithmetically from the published
parity bounds on `Z`; they are not an additional graph-theoretic assumption. -/
theorem belowZRegime_of_publishedZResults (hz : PublishedZResults) :
    ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d →
      FixedOrderRegime (d : ℤ) ((i : ℤ) - (d : ℤ)) := by
  intro d i hd hdi hiZ
  have hmod : d % 4 = 0 ∨ d % 4 = 1 ∨ d % 4 = 2 ∨ d % 4 = 3 := by omega
  rcases hmod with h0 | h1 | h2 | h3
  · have heven : Even d := ⟨d / 2, by omega⟩
    have hZ := hz.even d heven
    left
    omega
  · have hodd : Odd d := ⟨d / 2, by omega⟩
    have hZ := hz.oddUpper d hodd
    by_cases hlarge : 4 * ((i : ℤ) - (d : ℤ)) + 3 ≤ (d : ℤ)
    · exact Or.inl hlarge
    · right; right
      refine ⟨(i : ℤ) - (d : ℤ), ?_, ?_, rfl⟩ <;> omega
  · have heven : Even d := ⟨d / 2, by omega⟩
    have hZ := hz.even d heven
    left
    omega
  · by_cases hd7 : d = 7
    · subst d
      have hseven := Z_seven_le_nine
      left
      norm_num at hdi hiZ ⊢
      omega
    · have hodd : Odd d := ⟨d / 2, by omega⟩
      have hZ := hz.oddUpper d hodd
      by_cases hlarge : 4 * ((i : ℤ) - (d : ℤ)) + 3 ≤ (d : ℤ)
      · exact Or.inl hlarge
      · right; left
        refine ⟨(i : ℤ) - (d : ℤ), ?_, ?_, rfl⟩ <;> omega

/-- The coarser threshold inequalities required by packing are arithmetic
consequences of the same published `Z` results. -/
theorem thresholdBounds_of_publishedZResults (hz : PublishedZResults) :
    ∀ d : ℕ, 7 ≤ d →
      1 ≤ Z d ∧ Z d < 2 * d ∧ Z d ≤ d + d / 2 - 1 := by
  intro d hd
  have hlower := hz.lower d (by omega)
  have hmod : d % 4 = 0 ∨ d % 4 = 1 ∨ d % 4 = 2 ∨ d % 4 = 3 := by omega
  rcases hmod with h0 | h1 | h2 | h3
  · have heven : Even d := ⟨d / 2, by omega⟩
    have hZ := hz.even d heven
    omega
  · have hodd : Odd d := ⟨d / 2, by omega⟩
    have hZ := hz.oddUpper d hodd
    omega
  · have heven : Even d := ⟨d / 2, by omega⟩
    have hZ := hz.even d heven
    omega
  · by_cases hd7 : d = 7
    · subst d
      have hseven := Z_seven_le_nine
      omega
    · have hodd : Odd d := ⟨d / 2, by omega⟩
      have hZ := hz.oddUpper d hodd
      omega

/-- Inputs whose fields correspond directly to identified results in the
two supplied papers. The standard odd-cycle argument and Proposition 4.1 are
proved locally. -/
structure AuditedLiteratureInputs : Prop where
  publishedZ : PublishedZResults
  leftEndpoint : ∀ d : ℕ, 2 ≤ d → fDelta d d = d * d + 1
  thresholdEndpoint : ∀ d : ℕ, 2 ≤ d →
    fDelta d (Z d) = d * Z d + d / 2
  structuralDecomposition : ∀ d m : ℕ, 2 ≤ d → d ≤ m →
    StructuralDecomposition d m
  publishedSolvedCases : ∀ d m : ℕ, 2 ≤ d → (m ≤ d ∨ d < 7) →
    fDelta d m = conjecturedValue d m
  construction22Assembly : ∀ d m : ℕ, 7 ≤ d →
    (∀ r : ℕ, d < r → r < Z d →
      d * r + r - d + 1 ≤ fDelta d r) →
    conjecturedValue d m ≤ fDelta d m

/-- Proposition 4.1 applies throughout the open interval below `Z(d)`.
The published upper bound on `Z(d)` supplies its numerical side condition. -/
theorem proposition41_below_Z (hz : PublishedZResults) (d i : ℕ)
    (hd : 7 ≤ d) (hdi : d < i) (hiz : i < Z d) :
    d * i + i - d + 1 ≤ fDelta d i := by
  obtain ⟨_hzpos, _hz2, hzUpper⟩ := thresholdBounds_of_publishedZResults hz d hd
  have ht : 2 * (i - d) + 2 ≤ d := by omega
  have h := IntermediateConstruction.proposition41_lower_bound
    (d := d) (t := i - d) (by omega) ht
  have hi : d + (i - d) = i := by omega
  simpa [hi] using h

/-- Conjecture 2.1 from component-level, theorem-by-theorem literature inputs. -/
theorem conjecture21_of_audited_literature (lit : AuditedLiteratureInputs) :
    Conjecture21 := by
  have hregimes := belowZRegime_of_publishedZResults lit.publishedZ
  intro d i hd hdi hiz
  have hfixed : ∀ j : ℕ, d < j → j < Z d → FixedOrderComponentUpper d j := by
    intro j hdj hjz
    exact fixed_order_component_upper d j hd hdj
      (hregimes d j hd hdj hjz)
  have hcert := packingCertificate_of_fixed_components d i
    (lit.leftEndpoint d (by omega)) hfixed
    (lit.structuralDecomposition d i (by omega) (by omega))
  dsimp [PackingCertificate] at hcert
  obtain ⟨q₀, xs, hvalid, hweight, hupper⟩ := hcert
  obtain ⟨hzposN, hz2N, hzN⟩ :=
    thresholdBounds_of_publishedZResults lit.publishedZ d hd
  have hpacking := packing_utility_bound
    (d : ℤ) (Z d : ℤ) (d / 2 : ℕ) 0 (i : ℤ) q₀ xs
    (by omega) (by exact_mod_cast hzposN) (by positivity)
    (by exact_mod_cast hz2N)
    (by
      have hzN' : Z d ≤ d + d / 2 - 1 := hzN
      have hsum : 1 ≤ d + d / 2 := by omega
      have hh : (Z d : ℤ) ≤ ((d + d / 2 - 1 : ℕ) : ℤ) := by exact_mod_cast hzN'
      rw [Nat.cast_sub hsum, Nat.cast_add] at hh
      exact hh)
    (by positivity) (by exact_mod_cast hiz) hvalid (by simpa using hweight)
  have hres : residualUtility (d : ℤ) (i : ℤ) =
      (i : ℤ) - (d : ℤ) + 1 := by
    rw [residualUtility, if_neg]
    exact_mod_cast (Nat.not_lt.mpr (Nat.le_of_lt hdi))
  have hupperZ : (fDelta d i : ℤ) ≤
      (d : ℤ) * (i : ℤ) + (i : ℤ) - (d : ℤ) + 1 := by
    rw [hres] at hpacking
    push_cast at hpacking
    nlinarith
  have hformula := openComponentValue_natCast d i (Nat.le_of_lt hdi)
  have hupperN : fDelta d i ≤ d * i + i - d + 1 := by
    exact_mod_cast (hformula ▸ hupperZ)
  exact le_antisymm hupperN
    (proposition41_below_Z lit.publishedZ d i hd hdi hiz)

/-- Conjecture 2.2 from the same audited inputs.  The old broad interface is
constructed here as a theorem, rather than assumed. -/
theorem conjecture22_of_audited_literature (lit : AuditedLiteratureInputs) :
    Conjecture22 := by
  have hc21 := conjecture21_of_audited_literature lit
  have hregimes := belowZRegime_of_publishedZResults lit.publishedZ
  have hthresholds := thresholdBounds_of_publishedZResults lit.publishedZ
  apply conjecture22_of_literature_inputs
  exact {
    belowZRegime := hregimes
    construction21 := proposition41_below_Z lit.publishedZ
    structuralUpper21 := by
      intro hfixed d i hd hdi hiz
      exact (hc21 d i hd hdi hiz).le
    publishedSolvedCases := lit.publishedSolvedCases
    thresholdBounds := hthresholds
    construction22 := by
      intro d m hd
      exact lit.construction22Assembly d m hd
        (fun r hdr hrz => proposition41_below_Z lit.publishedZ d r hd hdr hrz)
    leftEndpoint := lit.leftEndpoint
    thresholdEndpoint := lit.thresholdEndpoint
    knapsackDecomposition := by
      intro d m hd hdm
      exact knapsackDecomposition_of_structural d m
        (lit.thresholdEndpoint d hd) (lit.structuralDecomposition d m hd hdm) }

end TriangleProof
