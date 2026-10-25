import TriangleProof.DegreeCases
import TriangleProof.HardCase

/-!
# Explicit interface for the standard initial-cycle argument

The proof in the accompanying draft uses the standard characterization of
bipartite graphs by odd cycles, followed by a shortest-odd-cycle argument.
`ShortestOddCycle.lean` proves the extraction and neighbor bound from Mathlib's
closed-walk parity characterization. `InitialCycle.lean` proves the statements
defined here, and `FixedOrder.lean` uses them without an external hypothesis.
-/

namespace TriangleProof

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- The graph-level initial-cycle statement, proved in `InitialCycle.lean`. -/
def DenseTriangleFreeFiveCycle : Prop :=
  ∀ (d t : ℤ),
    7 ≤ d → 1 ≤ t → t < d →
    (Fintype.card V : ℤ) = 2 * d + 2 * t + 1 →
    (∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a)) →
    (∀ v, degreeZ G v ≤ d) →
    d * (d + t) + t + 1 < edgeCountZ G →
    ∃ φ, IsCycleCopy G φ

/-- Canonical form of the initial-five-cycle principle for all finite graphs. -/
def InitialCyclePrinciple : Prop :=
  ∀ (n : ℕ) (G : SimpleGraph (Fin n)),
    @DenseTriangleFreeFiveCycle (Fin n) inferInstance G (Classical.decRel G.Adj)

/-- A compact description of the three numerical regimes needed by the
fixed-order proof.  The first is the easy degree-counting range; the other
two are the boundary classes handled by the obstruction arguments. -/
def FixedOrderRegime (d t : ℤ) : Prop :=
  4 * t + 3 ≤ d ∨
    (∃ a : ℤ, 3 ≤ a ∧ d = 4 * a - 1 ∧ t = a) ∨
    (∃ a : ℤ, 2 ≤ a ∧ d = 4 * a + 1 ∧ t = a)

end TriangleProof
