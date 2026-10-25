import TriangleProof.FixedOrder
import TriangleProof.InitialCycle
import TriangleProof.ExtremalFoundations
import TriangleProof.Packing

/-!
# Published inputs and the final reduction

No field below is an axiom declaration.  A theorem using `LiteratureInputs`
is conditional on a value of this structure, so `#print axioms` continues to
show only Lean's foundational axioms.  The fields isolate results to be cited
from the papers:

* `construction21`: Banak--Ekim--Taskin (2023), Proposition 4.1;
* `structuralUpper21`: their Lemma 2.3 together with the published endpoint
  component bounds;
* `belowZRegime`: the established formula for even `Z(d)` and the published
  odd-degree construction/bound.

The initial-cycle principle is proved in `InitialCycle.lean`, using the
shortest-cycle extraction and incidence bounds in `ShortestOddCycle.lean`.
It is not an external field.
-/

namespace TriangleProof

/-- The fixed-order component upper bound, stated on canonical finite vertex
types and without factor-criticality, which the checked proof does not need. -/
def FixedOrderComponentUpper (d i : ℕ) : Prop :=
  ∀ (n : ℕ) (G : SimpleGraph (Fin n)),
    letI : DecidableRel G.Adj := Classical.decRel G.Adj
    (Fintype.card (Fin n) : ℤ) = 2 * (i : ℤ) + 1 →
    (∀ a b c, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)) →
    (∀ v, degreeZ G v ≤ (d : ℤ)) →
    edgeCountZ G ≤ (d : ℤ) * (i : ℤ) + (i : ℤ) - (d : ℤ) + 1

theorem fixed_order_component_upper (d i : ℕ)
    (hd : 7 ≤ d) (hdi : d < i)
    (hregime : FixedOrderRegime (d : ℤ) ((i : ℤ) - (d : ℤ))) :
    FixedOrderComponentUpper d i := by
  intro n G
  let _ : DecidableRel G.Adj := Classical.decRel G.Adj
  intro hn htriangle hdeg
  let t : ℤ := (i : ℤ) - (d : ℤ)
  have hdZ : (7 : ℤ) ≤ (d : ℤ) := by exact_mod_cast hd
  have ht : (1 : ℤ) ≤ t := by dsimp [t]; omega
  have htd : t < (d : ℤ) := by
    dsimp [t]
    rcases hregime with hlarge | hboundary
    · omega
    · rcases hboundary with ⟨a, ha, hda, hta⟩ | ⟨a, ha, hda, hta⟩ <;> omega
  have hn' : (Fintype.card (Fin n) : ℤ) =
      2 * (d : ℤ) + 2 * t + 1 := by
    dsimp [t]
    omega
  have hub := fixed_order_edge_bound G htriangle
    (d : ℤ) t hdZ ht htd hn' hdeg hregime
  dsimp [t] at hub
  nlinarith

/-- Precisely the external inputs used to derive Conjecture 2.1 from the
checked fixed-order theorem. -/
structure LiteratureInputs21 : Prop where
  belowZRegime : ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d →
    FixedOrderRegime (d : ℤ) ((i : ℤ) - (d : ℤ))
  construction21 : ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d →
    d * i + i - d + 1 ≤ fDelta d i
  structuralUpper21 :
    (∀ d i : ℕ, 7 ≤ d → d < i → i < Z d → FixedOrderComponentUpper d i) →
    ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d →
      fDelta d i ≤ d * i + i - d + 1

theorem conjecture21_of_literature_inputs (lit : LiteratureInputs21) :
    Conjecture21 := by
  intro d i hd hdi hiz
  have hall : ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d →
      FixedOrderComponentUpper d i := by
    intro d i hd hdi hiz
    exact fixed_order_component_upper d i hd hdi
      (lit.belowZRegime d i hd hdi hiz)
  exact le_antisymm
    (lit.structuralUpper21 hall d i hd hdi hiz)
    (lit.construction21 d i hd hdi hiz)

theorem residualUtility_natCast (d r : ℕ) :
    ((if r < d then 0 else r - d + 1 : ℕ) : ℤ) =
      residualUtility (d : ℤ) (r : ℤ) := by
  simp only [residualUtility]
  by_cases h : r < d
  · simp [h]
  · have hdr : d ≤ r := Nat.le_of_not_gt h
    simp [h, hdr, Nat.cast_sub]

theorem conjecturedValue_natCast (d m : ℕ) :
    (conjecturedValue d m : ℤ) =
      (d : ℤ) * (m : ℤ) + ((m / Z d : ℕ) : ℤ) * ((d / 2 : ℕ) : ℤ) +
        residualUtility (d : ℤ) ((m % Z d : ℕ) : ℤ) := by
  unfold conjecturedValue
  rw [Nat.cast_add, Nat.cast_add, Nat.cast_mul, Nat.cast_mul,
    residualUtility_natCast]

/-- Numerical certificate produced by the published component reduction once
the fixed-order component values are known. -/
noncomputable def PackingCertificate (d m : ℕ) : Prop :=
  let z : ℤ := (Z d : ℕ)
  let b : ℤ := (d / 2 : ℕ)
  ∃ (q₀ : ℕ) (xs : List ℤ),
    (∀ j ∈ xs, (d : ℤ) ≤ j ∧ j < z) ∧
    (q₀ : ℤ) * z + xs.sum ≤ (m : ℤ) ∧
    (fDelta d m : ℤ) ≤ (d : ℤ) * (m : ℤ) +
      ((q₀ : ℤ) * b + (xs.map fun j => j - (d : ℤ) + 1).sum)

/-- Direct list form of Corollary 3.5 and Lemmas 4.5/4.9 of
Ahanjideh--Ekim--Yildiz (2022). Threshold components are counted by `q₀`;
`xs` records the matching numbers of all other non-star components. -/
noncomputable def KnapsackDecomposition (d m : ℕ) : Prop :=
  let z : ℤ := (Z d : ℕ)
  ∃ (q₀ : ℕ) (xs : List ℤ),
    (∀ j ∈ xs, (d : ℤ) ≤ j ∧ j < z) ∧
    (q₀ : ℤ) * z + xs.sum ≤ (m : ℤ) ∧
    (fDelta d m : ℤ) ≤ (d : ℤ) * (m : ℤ) +
      (q₀ : ℤ) * ((fDelta d (Z d) : ℤ) - (d : ℤ) * z) +
      (xs.map fun j => (fDelta d j.toNat : ℤ) - (d : ℤ) * j).sum

theorem packingCertificate_of_component_values
    (hc21 : Conjecture21) (d m : ℕ) (hd : 7 ≤ d)
    (hleft : fDelta d d = d * d + 1)
    (hthreshold : fDelta d (Z d) = d * Z d + d / 2)
    (hdecomp : KnapsackDecomposition d m) : PackingCertificate d m := by
  dsimp [KnapsackDecomposition] at hdecomp
  obtain ⟨q₀, xs, hvalid, hweight, hupper⟩ := hdecomp
  refine ⟨q₀, xs, hvalid, hweight, ?_⟩
  have hthresholdZ : (fDelta d (Z d) : ℤ) -
      (d : ℤ) * (Z d : ℤ) = (d / 2 : ℕ) := by
    rw [hthreshold]
    push_cast
    ring
  have hcomponent : ∀ j ∈ xs,
      (fDelta d j.toNat : ℤ) - (d : ℤ) * j = j - (d : ℤ) + 1 := by
    intro j hj
    have hjv := hvalid j hj
    have hjnonneg : 0 ≤ j := by omega
    have hjcast : (j.toNat : ℤ) = j := Int.toNat_of_nonneg hjnonneg
    by_cases hjd : j = (d : ℤ)
    · subst j
      simp only [Int.toNat_natCast]
      rw [hleft]
      push_cast
      ring
    · have hdjZ : (d : ℤ) < j := by omega
      have hdjZ' : (d : ℤ) < (j.toNat : ℤ) := by simpa [hjcast] using hdjZ
      have hdj : d < j.toNat := by exact_mod_cast hdjZ'
      have hjz : j.toNat < Z d := by
        have hjzZ : (j.toNat : ℤ) < (Z d : ℤ) := by
          simpa [hjcast] using hjv.2
        exact_mod_cast hjzZ
      have hv := hc21 d j.toNat hd hdj hjz
      dsimp [Conjecture21At] at hv
      have hsub : d ≤ d * j.toNat + j.toNat := by nlinarith
      rw [hv, Nat.cast_add, Nat.cast_sub hsub, Nat.cast_add, Nat.cast_mul]
      rw [hjcast]
      ring
  have hlist :
      (xs.map fun j => (fDelta d j.toNat : ℤ) - (d : ℤ) * j).sum =
      (xs.map fun j => j - (d : ℤ) + 1).sum := by
    exact congrArg List.sum (List.map_congr_left hcomponent)
  rw [hthresholdZ, hlist] at hupper
  simpa only [add_assoc] using hupper

/-- External results used only in the final passage from component values to
the all-`m` formula.  `packingCertificate` is the formal interface to the
component-reduction lemma; `construction22` records realizability by the
published components and stars. -/
structure LiteratureInputs22 : Prop extends LiteratureInputs21 where
  publishedSolvedCases : ∀ d m : ℕ, 2 ≤ d → (m ≤ d ∨ d < 7) →
    fDelta d m = conjecturedValue d m
  thresholdBounds : ∀ d : ℕ, 7 ≤ d →
    1 ≤ Z d ∧ Z d < 2 * d ∧ Z d ≤ d + d / 2 - 1
  construction22 : ∀ d m : ℕ, 7 ≤ d →
    conjecturedValue d m ≤ fDelta d m
  leftEndpoint : ∀ d : ℕ, 2 ≤ d → fDelta d d = d * d + 1
  thresholdEndpoint : ∀ d : ℕ, 2 ≤ d →
    fDelta d (Z d) = d * Z d + d / 2
  knapsackDecomposition : ∀ d m : ℕ, 2 ≤ d → d ≤ m →
    KnapsackDecomposition d m

theorem conjecture22_of_literature_inputs (lit : LiteratureInputs22) :
    Conjecture22 := by
  have hc21 : Conjecture21 :=
    conjecture21_of_literature_inputs lit.toLiteratureInputs21
  intro d m hd
  by_cases hsolved : m ≤ d ∨ d < 7
  · exact lit.publishedSolvedCases d m hd hsolved
  have hd7 : 7 ≤ d := by omega
  have hdm : d ≤ m := by omega
  obtain ⟨hzposN, hz2N, hzN⟩ := lit.thresholdBounds d hd7
  let zN : ℕ := Z d
  let bN : ℕ := d / 2
  let kN : ℕ := m / zN
  let rN : ℕ := m % zN
  have hzpos : (1 : ℤ) ≤ (zN : ℤ) := by exact_mod_cast hzposN
  have hz2 : (zN : ℤ) < 2 * (d : ℤ) := by exact_mod_cast hz2N
  have hb : (0 : ℤ) ≤ (bN : ℤ) := by positivity
  have hzN' : zN ≤ d + bN - 1 := by simpa [zN, bN] using hzN
  have hsum : 1 ≤ d + bN := by omega
  have hz : (zN : ℤ) ≤ (d : ℤ) + (bN : ℤ) - 1 := by
    have hh : (zN : ℤ) ≤ ((d + bN - 1 : ℕ) : ℤ) := by exact_mod_cast hzN'
    rw [Nat.cast_sub hsum, Nat.cast_add] at hh
    exact hh
  have hr0 : (0 : ℤ) ≤ (rN : ℤ) := by positivity
  have hrzN : rN < zN := by
    dsimp [rN, zN]
    exact Nat.mod_lt m hzposN
  have hrz : (rN : ℤ) < (zN : ℤ) := by exact_mod_cast hrzN
  have hmN : kN * zN + rN = m := by
    simpa [kN, rN, mul_comm] using Nat.div_add_mod m zN
  have hm : (kN : ℤ) * (zN : ℤ) + (rN : ℤ) = (m : ℤ) := by
    exact_mod_cast hmN
  have hlower0 := lit.construction22 d m hd7
  have hlowerCast : (conjecturedValue d m : ℤ) ≤ (fDelta d m : ℤ) := by
    exact_mod_cast hlower0
  have hcv := conjecturedValue_natCast d m
  have hlower : (d : ℤ) * ((kN : ℤ) * (zN : ℤ) + (rN : ℤ)) +
      (kN : ℤ) * (bN : ℤ) + residualUtility (d : ℤ) (rN : ℤ) ≤
      (fDelta d m : ℤ) := by
    rw [hm]
    rw [← hcv]
    exact hlowerCast
  have hcert0 := packingCertificate_of_component_values hc21 d m hd7
    (lit.leftEndpoint d (by omega)) (lit.thresholdEndpoint d (by omega))
    (lit.knapsackDecomposition d m (by omega) hdm)
  dsimp [PackingCertificate, zN, bN] at hcert0
  obtain ⟨q₀, xs, hvalid, hweight0, hupper0⟩ := hcert0
  have hweight : (q₀ : ℤ) * (zN : ℤ) + xs.sum ≤
      (kN : ℤ) * (zN : ℤ) + (rN : ℤ) := by
    rw [hm]
    simpa [zN] using hweight0
  have hupper : (fDelta d m : ℤ) ≤
      (d : ℤ) * ((kN : ℤ) * (zN : ℤ) + (rN : ℤ)) +
        ((q₀ : ℤ) * (bN : ℤ) +
          (xs.map fun j => j - (d : ℤ) + 1).sum) := by
    rw [hm]
    simpa [bN] using hupper0
  have hF := extremal_value_from_packing
    (d : ℤ) (zN : ℤ) (bN : ℤ) (kN : ℤ) (rN : ℤ) (fDelta d m : ℤ)
    (by omega) hzpos hb hz2 hz hr0 hrz hlower
    ⟨q₀, xs, by simpa [zN] using hvalid, hweight, hupper⟩
  apply Int.ofNat_inj.mp
  rw [conjecturedValue_natCast]
  dsimp [zN, bN, kN, rN] at hm hF ⊢
  nlinarith

end TriangleProof
