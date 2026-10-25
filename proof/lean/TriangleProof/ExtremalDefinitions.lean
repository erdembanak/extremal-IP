import TriangleProof.GraphCounting
import Mathlib.Combinatorics.SimpleGraph.Matching

/-!
# Definitions corresponding to the two conjectures

All finite graphs can be represented on `Fin n`, so the extremal sets below
quantify over those canonical vertex types.  This avoids quotienting graphs
by isomorphism while retaining the paper's definitions.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V]
variable (G : SimpleGraph V)

def IsTriangleFree : Prop :=
  ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)

noncomputable def degreeN (v : V) : ℕ := (G.neighborSet v).ncard

noncomputable def HasMaximumDegreeAtMost (d : ℕ) : Prop := ∀ v, degreeN G v ≤ d

/-- Maximum cardinality in edges of a matching.  A matching has twice as
many covered vertices as edges, hence the division by two. -/
noncomputable def matchingNumber : ℕ :=
  sSup {k | ∃ M : G.Subgraph, M.IsMatching ∧ k = M.verts.ncard / 2}

noncomputable def IsAdmissible (d m : ℕ) : Prop :=
  IsTriangleFree G ∧ HasMaximumDegreeAtMost G d ∧ matchingNumber G ≤ m

/-- Edge counts attained by admissible finite graphs. -/
noncomputable def attainableEdgeCounts (d m : ℕ) : Set ℕ :=
  {e | ∃ n : ℕ, ∃ G : SimpleGraph (Fin n),
    IsAdmissible G d m ∧ G.edgeSet.ncard = e}

/-- The extremal function from the papers. `ExtremalFoundations` proves that
the set is bounded and nonempty, and that this supremum is attained. -/
noncomputable def fDelta (d m : ℕ) : ℕ := sSup (attainableEdgeCounts d m)

def HasPerfectMatching {W : Type*} (H : SimpleGraph W) : Prop :=
  ∃ M : H.Subgraph, M.IsPerfectMatching

/-- Deleting any one vertex leaves a graph with a perfect matching. -/
def IsFactorCritical : Prop :=
  ∀ v : V, HasPerfectMatching (G.induce {w | w ≠ v})

/-- Regular for even `d`, and total degree deficiency one for odd `d`.
Together with the maximum-degree condition, the latter says exactly one
vertex has degree `d-1` and all others have degree `d`. -/
noncomputable def HasThresholdDegreeProfile (d : ℕ) : Prop :=
  (Even d ∧ ∀ v, degreeN G v = d) ∨
  (Odd d ∧ HasMaximumDegreeAtMost G d ∧
    ∑ v, (d - degreeN G v) = 1)

noncomputable def IsZCandidate (d i : ℕ) : Prop :=
  ∃ G : SimpleGraph (Fin (2 * i + 1)),
    IsTriangleFree G ∧ IsFactorCritical G ∧
    matchingNumber G = i ∧ HasThresholdDegreeProfile G d

/-- The least matching number of a threshold component. `ExtremalFoundations`
turns the published existence construction into a proof that this infimum is
an actual candidate. -/
noncomputable def Z (d : ℕ) : ℕ := sInf {i | IsZCandidate d i}

/-- The pointwise statement of Conjecture 2.1. -/
noncomputable def Conjecture21At (d i : ℕ) : Prop :=
  fDelta d i = d * i + i - d + 1

noncomputable def Conjecture21 : Prop :=
  ∀ d i : ℕ, 7 ≤ d → d < i → i < Z d → Conjecture21At d i

noncomputable def conjecturedValue (d m : ℕ) : ℕ :=
  let z := Z d
  let k := m / z
  let r := m % z
  d * m + k * (d / 2) + if r < d then 0 else r - d + 1

/-- The pointwise statement of Conjecture 2.2. -/
noncomputable def Conjecture22 : Prop :=
  ∀ d m : ℕ, 2 ≤ d → fDelta d m = conjecturedValue d m

end TriangleProof
