import Mathlib.Combinatorics.SimpleGraph.Finite
import Mathlib.Data.Finset.Powerset
import Mathlib.Tactic

/-!
# The finite obstructions in Section 6 of the proof draft

All proofs below use Lean's ordinary kernel-checked tactics. The finite
decisions use `decide`, not native evaluation. These theorems concern the
specific eight-vertex graphs; they do not assert the extremal conjecture.
-/

namespace TriangleProof

abbrev Vertex := Fin 8

private def obstructionEdges : List (Vertex × Vertex) :=
  [(0, 2), (0, 3), (1, 2), (1, 3),
   (0, 4), (4, 5), (5, 1), (2, 6), (6, 7), (7, 3)]

/-- Labels `0,1,2,3,4,5,6,7` denote `a,b,u,v,p,q,x,y`. -/
def obstructionRel (i j : Vertex) : Prop :=
  (i, j) ∈ obstructionEdges ∨ (j, i) ∈ obstructionEdges

instance (i j : Vertex) : Decidable (obstructionRel i j) :=
  inferInstanceAs (Decidable ((i, j) ∈ obstructionEdges ∨ (j, i) ∈ obstructionEdges))

def obstruction : SimpleGraph Vertex where
  Adj := obstructionRel
  symm := ⟨by decide⟩
  loopless := ⟨by decide⟩

instance : DecidableRel obstruction.Adj :=
  fun i j => inferInstanceAs (Decidable (obstructionRel i j))

def mobiusRel (i j : Vertex) : Prop :=
  (i.val + 1) % 8 = j.val ∨
  (j.val + 1) % 8 = i.val ∨
  (i.val + 4) % 8 = j.val

instance (i j : Vertex) : Decidable (mobiusRel i j) :=
  inferInstanceAs (Decidable (
    (i.val + 1) % 8 = j.val ∨
    (j.val + 1) % 8 = i.val ∨
    (i.val + 4) % 8 = j.val))

/-- The eight-cycle with its four opposite chords. -/
def mobius : SimpleGraph Vertex where
  Adj := mobiusRel
  symm := ⟨by decide⟩
  loopless := ⟨by decide⟩

instance : DecidableRel mobius.Adj :=
  fun i j => inferInstanceAs (Decidable (mobiusRel i j))

def Independent (G : SimpleGraph Vertex) (s : Finset Vertex) : Prop :=
  ∀ u ∈ s, ∀ v ∈ s, ¬ G.Adj u v

instance (G : SimpleGraph Vertex) [DecidableRel G.Adj] (s : Finset Vertex) :
    Decidable (Independent G s) :=
  inferInstanceAs (Decidable (∀ u ∈ s, ∀ v ∈ s, ¬ G.Adj u v))

set_option maxRecDepth 100000
set_option maxHeartbeats 0

theorem obstruction_triangle_free :
    ∀ a b c : Vertex,
      ¬ (obstruction.Adj a b ∧ obstruction.Adj b c ∧ obstruction.Adj c a) := by
  decide

theorem mobius_triangle_free :
    ∀ a b c : Vertex,
      ¬ (mobius.Adj a b ∧ mobius.Adj b c ∧ mobius.Adj c a) := by
  decide

theorem obstruction_independence_bound :
    ∀ s : Finset Vertex, Independent obstruction s → s.card ≤ 3 := by
  decide

theorem mobius_independence_bound :
    ∀ s : Finset Vertex, Independent mobius s → s.card ≤ 3 := by
  decide

theorem mobius_independent_triples :
    ∀ s : Finset Vertex,
      (Independent mobius s ∧ s.card = 3) ↔
      ∃ v : Vertex, s = mobius.neighborFinset v := by
  decide

theorem mobius_degree_three :
    ∀ v : Vertex, (mobius.neighborFinset v).card = 3 := by
  decide

theorem mobius_common_neighbor :
    ∀ u v : Vertex, u ≠ v → ¬ mobius.Adj u v →
      ∃ w : Vertex, mobius.Adj u w ∧ mobius.Adj v w := by
  decide

/-- In any triangle-free graph, at most three labels of an `H` copy can
be adjacent to a given vertex. Injectivity is not needed for this bound. -/
theorem obstruction_copy_neighbor_bound {V : Type*} (G : SimpleGraph V)
    [DecidableRel G.Adj]
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V)
    (hedges : ∀ i j : Vertex, obstruction.Adj i j → G.Adj (φ i) (φ j))
    (v : V) :
    (Finset.univ.filter (fun i : Vertex => G.Adj v (φ i))).card ≤ 3 := by
  apply obstruction_independence_bound
  intro i hi j hj hij
  have hvi : G.Adj v (φ i) := (Finset.mem_filter.mp hi).2
  have hvj : G.Adj v (φ j) := (Finset.mem_filter.mp hj).2
  exact htriangle v (φ i) (φ j) ⟨hvi, hedges i j hij, hvj.symm⟩

/-- The analogous bound for a Möbius-ladder copy in an arbitrary graph. -/
theorem mobius_copy_neighbor_bound {V : Type*} (G : SimpleGraph V)
    [DecidableRel G.Adj]
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : Vertex → V)
    (hedges : ∀ i j : Vertex, mobius.Adj i j → G.Adj (φ i) (φ j))
    (v : V) :
    (Finset.univ.filter (fun i : Vertex => G.Adj v (φ i))).card ≤ 3 := by
  apply mobius_independence_bound
  intro i hi j hj hij
  have hvi : G.Adj v (φ i) := (Finset.mem_filter.mp hi).2
  have hvj : G.Adj v (φ j) := (Finset.mem_filter.mp hj).2
  exact htriangle v (φ i) (φ j) ⟨hvi, hedges i j hij, hvj.symm⟩

/-- Bag-neighbor sums in an eight-vertex Möbius-ladder blow-up. -/
def mobiusPotential (a : Vertex → ℤ) (j : Vertex) : ℤ :=
  a (j - 1) + a (j + 1) + a (j + 4)

theorem mobius_bag_identity (a : Vertex → ℤ) (j : Vertex) :
    mobiusPotential a (j - 1) + mobiusPotential a j + mobiusPotential a (j + 1)
      = (∑ k : Vertex, a k) + a j := by
  fin_cases j <;>
    simp only [mobiusPotential, Fin.sum_univ_succ] <;>
    dsimp [Fin.add, Fin.sub, Fin.succ] <;> ring

end TriangleProof
