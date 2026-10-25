import TriangleProof.SaturatedCycle
import TriangleProof.ObstructionExclusion

namespace TriangleProof

set_option maxRecDepth 100000
set_option maxHeartbeats 0

theorem cycle_independent_pair : ∀ s, cycleIndependent s → s.card = 2 →
    ∃ j, s = fiveCycle.neighborFinset j := by decide

theorem cycle_nonedge_common_neighbor : ∀ i j : CycleVertex,
    ¬ fiveCycle.Adj i j → ∃ k, fiveCycle.Adj i k ∧ fiveCycle.Adj j k := by decide

theorem cycle_add_one_adj : ∀ j : CycleVertex, fiveCycle.Adj j (j + 1) := by decide

theorem cycle_rotate_injective : ∀ j : CycleVertex, Function.Injective (fun i : CycleVertex => j + i) := by decide

theorem cycle_rotate_adj : ∀ j a b : CycleVertex,
    fiveCycle.Adj (j + a) (j + b) ↔ fiveCycle.Adj a b := by decide

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

def cycleNeighbors (φ : CycleVertex → V) (v : V) : Finset CycleVertex :=
  Finset.univ.filter (fun j => G.Adj v (φ j))

def cycleBag (φ : CycleVertex → V) (j : CycleVertex) : Finset V :=
  Finset.univ.filter (fun v => cycleNeighbors G φ v = fiveCycle.neighborFinset j)

def cycleRemainder (φ : CycleVertex → V) : Finset V :=
  Finset.univ.filter (fun v => (cycleNeighbors G φ v).card < 2)

omit [Fintype V] [DecidableEq V] in
theorem cycleNeighbors_independent
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (v : V) :
    cycleIndependent (cycleNeighbors G φ v) := by
  intro i hi j hj hij
  exact htriangle v (φ i) (φ j)
    ⟨(Finset.mem_filter.mp hi).2, hφ.2 i j hij, ((Finset.mem_filter.mp hj).2).symm⟩

omit [Fintype V] [DecidableEq V] in
theorem cycleNeighbors_on_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (j : CycleVertex) :
    cycleNeighbors G φ (φ j) = fiveCycle.neighborFinset j := by
  apply cycle_pattern_of_containing_neighborhood _ j (cycleNeighbors_independent G htriangle φ hφ _)
  intro k hk
  exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hφ.2 j k hk⟩

omit [DecidableEq V] in
theorem cycleBag_pattern (φ : CycleVertex → V) {v : V} {j : CycleVertex}
    (hv : v ∈ cycleBag G φ j) (k : CycleVertex) :
    G.Adj v (φ k) ↔ fiveCycle.Adj j k := by
  have hh := Finset.ext_iff.mp (Finset.mem_filter.mp hv).2 k
  simpa [cycleNeighbors] using hh

omit [DecidableEq V] in
theorem cycleBag_unique (φ : CycleVertex → V) {v : V} {i j : CycleVertex}
    (hi : v ∈ cycleBag G φ i) (hj : v ∈ cycleBag G φ j) : i = j := by
  exact cycle_neighborhood_injective ((Finset.mem_filter.mp hi).2.symm.trans (Finset.mem_filter.mp hj).2)

omit [DecidableEq V] in
theorem cycleBag_cycle_vertex
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (j : CycleVertex) :
    φ j ∈ cycleBag G φ j := by
  exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, cycleNeighbors_on_cycle G htriangle φ hφ j⟩

omit [DecidableEq V] in
theorem cycleRemainder_not_cycle
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) {v : V}
    (hv : v ∈ cycleRemainder G φ) : v ∉ Set.range φ := by
  rintro ⟨j, rfl⟩
  have hh := (Finset.mem_filter.mp hv).2
  rw [cycleNeighbors_on_cycle G htriangle φ hφ, cycle_neighborhood_card] at hh
  omega

omit [DecidableEq V] in
theorem cycleRemainder_not_bag (φ : CycleVertex → V) {v : V}
    (hv : v ∈ cycleRemainder G φ) (j : CycleVertex) : v ∉ cycleBag G φ j := by
  intro hb
  have hh := (Finset.mem_filter.mp hv).2
  rw [(Finset.mem_filter.mp hb).2, cycle_neighborhood_card] at hh
  omega

omit [DecidableEq V] in
theorem cycle_partition_covers
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (v : V)
    (hv : v ∉ cycleRemainder G φ) : ∃ j, v ∈ cycleBag G φ j := by
  have hind := cycleNeighbors_independent G htriangle φ hφ v
  have hcard := cycle_independent_bound _ hind
  have hnot : ¬ (cycleNeighbors G φ v).card < 2 := by simpa [cycleRemainder] using hv
  obtain ⟨j, hj⟩ := cycle_independent_pair _ hind (by omega)
  exact ⟨j, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hj⟩⟩

omit [DecidableEq V] in
theorem cycle_bag_edges
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) {v w : V} {i j : CycleVertex}
    (hv : v ∈ cycleBag G φ i) (hw : w ∈ cycleBag G φ j) (hvw : G.Adj v w) :
    fiveCycle.Adj i j := by
  by_contra h
  obtain ⟨k, hik, hjk⟩ := cycle_nonedge_common_neighbor i j h
  exact htriangle v w (φ k)
    ⟨hvw, (cycleBag_pattern G φ hw k).mpr hjk, ((cycleBag_pattern G φ hv k).mpr hik).symm⟩

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- Ten edge checks suffice to exhibit a copy of the obstruction. -/
theorem obstruction_map_edges (f : Vertex → V)
    (h02 : G.Adj (f 0) (f 2)) (h03 : G.Adj (f 0) (f 3))
    (h12 : G.Adj (f 1) (f 2)) (h13 : G.Adj (f 1) (f 3))
    (h04 : G.Adj (f 0) (f 4)) (h45 : G.Adj (f 4) (f 5)) (h51 : G.Adj (f 5) (f 1))
    (h26 : G.Adj (f 2) (f 6)) (h67 : G.Adj (f 6) (f 7)) (h73 : G.Adj (f 7) (f 3)) :
    ∀ i j, obstruction.Adj i j → G.Adj (f i) (f j) := by
  intro i j hij
  fin_cases i <;> fin_cases j <;>
    first | exact h02 | exact h03 | exact h12 | exact h13 | exact h04 | exact h45 | exact h51 |
      exact h26 | exact h67 | exact h73 | exact h02.symm | exact h03.symm | exact h12.symm |
      exact h13.symm | exact h04.symm | exact h45.symm | exact h51.symm | exact h26.symm |
      exact h67.symm | exact h73.symm | (exfalso; revert hij; decide)

omit [DecidableEq V] in
theorem cycle_no_adjacent_contacts
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (x u v : V) (j : CycleVertex) (hx : x ∈ cycleRemainder G φ)
    (hu : u ∈ cycleBag G φ j) (hv : v ∈ cycleBag G φ (j + 1))
    (hxu : G.Adj x u) (hxv : G.Adj x v) : False := by
  classical
  let C := fun k : CycleVertex => φ (j + k)
  have hC : Function.Injective C := hφ.1.comp (cycle_rotate_injective j)
  have hCE : ∀ a b, fiveCycle.Adj a b → G.Adj (C a) (C b) := by
    intro a b hab
    exact hφ.2 _ _ ((cycle_rotate_adj j a b).mpr hab)
  have huN : ∀ k, G.Adj u (C k) ↔ fiveCycle.Adj 0 k := by
    intro k
    rw [cycleBag_pattern G φ hu]
    simpa using cycle_rotate_adj j 0 k
  have hvN : ∀ k, G.Adj v (C k) ↔ fiveCycle.Adj 1 k := by
    intro k
    rw [cycleBag_pattern G φ hv]
    exact cycle_rotate_adj j 1 k
  have hxc : ∀ k, x ≠ C k := by
    intro k hk
    exact cycleRemainder_not_cycle G htriangle φ hφ hx ⟨j + k, hk.symm⟩
  have huc : ∀ k, u ≠ C k := by
    intro k hk
    have hkj := cycleBag_unique G φ hu (by rw [hk]; exact cycleBag_cycle_vertex G htriangle φ hφ (j + k))
    have hu0 : u = C 0 := by simpa only [C, ← hkj, add_zero] using hk
    have huv : G.Adj u v := by rw [hu0]; exact ((hvN 0).mpr (by decide)).symm
    exact htriangle x u v ⟨hxu, huv, hxv.symm⟩
  have hvc : ∀ k, v ≠ C k := by
    intro k hk
    have hkj := cycleBag_unique G φ hv (by rw [hk]; exact cycleBag_cycle_vertex G htriangle φ hφ (j + k))
    have hv1 : v = C 1 := by simpa only [C, ← hkj] using hk
    have huv : G.Adj u v := by rw [hv1]; exact (huN 1).mpr (by decide)
    exact htriangle x u v ⟨hxu, huv, hxv.symm⟩
  have huv : u ≠ v := by
    intro huv
    have hh := cycleBag_unique G φ hu (by rw [huv]; exact hv)
    have hn : ∀ j : CycleVertex, j ≠ j + 1 := by decide
    exact hn j hh
  have hxu' : x ≠ u := hxu.ne
  have hxv' : x ≠ v := hxv.ne
  let f : Vertex → V := ![C 0, C 2, v, C 1, C 4, C 3, x, u]
  have hf : Function.Injective f := by
    intro a b hab
    fin_cases a <;> fin_cases b <;>
      dsimp only [f] at hab <;>
      first | rfl | (exfalso; exact hxc _ hab) | (exfalso; exact hxc _ hab.symm) |
        (exfalso; exact huc _ hab) | (exfalso; exact huc _ hab.symm) |
        (exfalso; exact hvc _ hab) | (exfalso; exact hvc _ hab.symm) |
        (exfalso; exact huv hab) | (exfalso; exact huv hab.symm) |
        (exfalso; exact hxu' hab) | (exfalso; exact hxu' hab.symm) |
        (exfalso; exact hxv' hab) | (exfalso; exact hxv' hab.symm) |
        (have hh := congrArg Fin.val (hC hab); norm_num at hh)
  apply hnoH f hf
  apply obstruction_map_edges
  · exact ((hvN 0).mpr (by decide)).symm
  · exact hCE 0 1 (by decide)
  · exact ((hvN 2).mpr (by decide)).symm
  · exact hCE 2 1 (by decide)
  · exact hCE 0 4 (by decide)
  · exact hCE 4 3 (by decide)
  · exact hCE 3 2 (by decide)
  · exact hxv.symm
  · exact hxu
  · exact (huN 1).mpr (by decide)

omit [DecidableEq V] in
theorem cycle_no_shared_bag_on_remainder_edge
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (x y u v : V) (j : CycleVertex)
    (hx : x ∈ cycleRemainder G φ) (hy : y ∈ cycleRemainder G φ)
    (hu : u ∈ cycleBag G φ j) (hv : v ∈ cycleBag G φ j)
    (hxy : G.Adj x y) (hxu : G.Adj x u) (hyv : G.Adj y v) : False := by
  classical
  let C := fun k : CycleVertex => φ (j + k)
  have hC : Function.Injective C := hφ.1.comp (cycle_rotate_injective j)
  have hCE : ∀ a b, fiveCycle.Adj a b → G.Adj (C a) (C b) := by
    intro a b hab
    exact hφ.2 _ _ ((cycle_rotate_adj j a b).mpr hab)
  have huN : ∀ k, G.Adj u (C k) ↔ fiveCycle.Adj 0 k := by
    intro k
    rw [cycleBag_pattern G φ hu]
    simpa using cycle_rotate_adj j 0 k
  have hvN : ∀ k, G.Adj v (C k) ↔ fiveCycle.Adj 0 k := by
    intro k
    rw [cycleBag_pattern G φ hv]
    simpa using cycle_rotate_adj j 0 k
  have hxc : ∀ k, x ≠ C k := by
    intro k hk
    exact cycleRemainder_not_cycle G htriangle φ hφ hx ⟨j + k, hk.symm⟩
  have hyc : ∀ k, y ≠ C k := by
    intro k hk
    exact cycleRemainder_not_cycle G htriangle φ hφ hy ⟨j + k, hk.symm⟩
  have huc : ∀ k, k ≠ 0 → u ≠ C k := by
    intro k hk heq
    have hjk := cycleBag_unique G φ hu (by rw [heq]; exact cycleBag_cycle_vertex G htriangle φ hφ (j + k))
    have hh : j + k = j + 0 := by simpa using hjk.symm
    exact hk (cycle_rotate_injective j hh)
  have hvc : ∀ k, k ≠ 0 → v ≠ C k := by
    intro k hk heq
    have hjk := cycleBag_unique G φ hv (by rw [heq]; exact cycleBag_cycle_vertex G htriangle φ hφ (j + k))
    have hh : j + k = j + 0 := by simpa using hjk.symm
    exact hk (cycle_rotate_injective j hh)
  have huv : u ≠ v := by
    intro huv
    exact htriangle x y u ⟨hxy, huv ▸ hyv, hxu.symm⟩
  have hxu' : x ≠ u := hxu.ne
  have hyv' : y ≠ v := hyv.ne
  have hxv : x ≠ v := by
    intro h
    exact cycleRemainder_not_bag G φ hx j (by rw [h]; exact hv)
  have hyu : y ≠ u := by
    intro h
    exact cycleRemainder_not_bag G φ hy j (by rw [h]; exact hu)
  have hxy' : x ≠ y := hxy.ne
  let f : Vertex → V := ![C 4, C 1, u, v, C 3, C 2, x, y]
  have hf : Function.Injective f := by
    intro a b hab
    fin_cases a <;> fin_cases b <;> dsimp only [f] at hab <;>
      first | rfl | (exfalso; exact hxc _ hab) | (exfalso; exact hxc _ hab.symm) |
        (exfalso; exact hyc _ hab) | (exfalso; exact hyc _ hab.symm) |
        (exfalso; exact huc _ (by decide) hab) | (exfalso; exact huc _ (by decide) hab.symm) |
        (exfalso; exact hvc _ (by decide) hab) | (exfalso; exact hvc _ (by decide) hab.symm) |
        (exfalso; exact huv hab) | (exfalso; exact huv hab.symm) |
        (exfalso; exact hxu' hab) | (exfalso; exact hxu' hab.symm) |
        (exfalso; exact hyv' hab) | (exfalso; exact hyv' hab.symm) |
        (exfalso; exact hxv hab) | (exfalso; exact hxv hab.symm) |
        (exfalso; exact hyu hab) | (exfalso; exact hyu hab.symm) |
        (exfalso; exact hxy' hab) | (exfalso; exact hxy' hab.symm) |
        (have hh := congrArg Fin.val (hC hab); norm_num at hh)
  apply hnoH f hf
  apply obstruction_map_edges
  · exact ((huN 4).mpr (by decide)).symm
  · exact ((hvN 4).mpr (by decide)).symm
  · exact ((huN 1).mpr (by decide)).symm
  · exact ((hvN 1).mpr (by decide)).symm
  · exact hCE 4 3 (by decide)
  · exact hCE 3 2 (by decide)
  · exact hCE 2 1 (by decide)
  · exact hxu.symm
  · exact hxy
  · exact hyv

def cycleSupport (φ : CycleVertex → V) (v : V) : Finset CycleVertex :=
  Finset.univ.filter (fun j => ∃ u ∈ cycleBag G φ j, G.Adj v u)

omit [DecidableEq V] in
theorem cycle_support_independent
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (x : V)
    (hx : x ∈ cycleRemainder G φ) : cycleIndependent (cycleSupport G φ x) := by
  intro i hi j hj hij
  obtain ⟨u, hu, hxu⟩ := (Finset.mem_filter.mp hi).2
  obtain ⟨v, hv, hxv⟩ := (Finset.mem_filter.mp hj).2
  have hadj : ∀ i j : CycleVertex, fiveCycle.Adj i j → j = i + 1 ∨ i = j + 1 := by decide
  rcases hadj i j hij with rfl | rfl
  · exact cycle_no_adjacent_contacts G htriangle hnoH φ hφ x u v i hx hu hv hxu hxv
  · exact cycle_no_adjacent_contacts G htriangle hnoH φ hφ x v u j hx hv hu hxv hxu

/-- The two obstruction embeddings suffice to extend the partition.
In particular, no separate same-color degree argument is necessary. -/
theorem cycle_partition_extension
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (hnoH : ∀ f : Vertex → V, Function.Injective f →
      ¬ (∀ i j, obstruction.Adj i j → G.Adj (f i) (f j)))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ)
    (hmany : ∀ x ∈ cycleRemainder G φ, 2 ≤ (cycleSupport G φ x).card) :
    ∃ c : V → CycleVertex, (∀ j, c (φ j) = j) ∧
      (∀ v w, G.Adj v w → fiveCycle.Adj (c v) (c w)) := by
  classical
  have hassign : ∀ v, ∃ j,
      (v ∈ cycleRemainder G φ → cycleSupport G φ v = fiveCycle.neighborFinset j) ∧
      (v ∉ cycleRemainder G φ → v ∈ cycleBag G φ j) := by
    intro v
    by_cases hv : v ∈ cycleRemainder G φ
    · have hi := cycle_support_independent G htriangle hnoH φ hφ v hv
      have hb := cycle_independent_bound _ hi
      have hl := hmany v hv
      obtain ⟨j, hj⟩ := cycle_independent_pair _ hi (by omega)
      exact ⟨j, fun _ => hj, fun hn => (hn hv).elim⟩
    · obtain ⟨j, hj⟩ := cycle_partition_covers G htriangle φ hφ v hv
      exact ⟨j, fun hh => (hv hh).elim, fun _ => hj⟩
  choose c hc using hassign
  have hcycle : ∀ j, c (φ j) = j := by
    intro j
    have hb := cycleBag_cycle_vertex G htriangle φ hφ j
    have hn : φ j ∉ cycleRemainder G φ := by
      intro h
      exact cycleRemainder_not_bag G φ h j hb
    exact cycleBag_unique G φ ((hc (φ j)).2 hn) hb
  refine ⟨c, hcycle, fun v w hvw => ?_⟩
  by_cases hv : v ∈ cycleRemainder G φ <;> by_cases hw : w ∈ cycleRemainder G φ
  · by_contra hn
    obtain ⟨k, hvk, hwk⟩ := cycle_nonedge_common_neighbor (c v) (c w) hn
    have hvs : k ∈ cycleSupport G φ v := by rw [(hc v).1 hv]; simpa using hvk
    have hws : k ∈ cycleSupport G φ w := by rw [(hc w).1 hw]; simpa using hwk
    obtain ⟨u, hu, hvu⟩ := (Finset.mem_filter.mp hvs).2
    obtain ⟨z, hz, hwz⟩ := (Finset.mem_filter.mp hws).2
    exact cycle_no_shared_bag_on_remainder_edge G htriangle hnoH φ hφ v w u z k hv hw hu hz hvw hvu hwz
  · have hs : c w ∈ cycleSupport G φ v :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, w, (hc w).2 hw, hvw⟩
    rw [(hc v).1 hv] at hs
    simpa using hs
  · have hs : c v ∈ cycleSupport G φ w :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, v, (hc v).2 hv, hvw.symm⟩
    rw [(hc w).1 hw] at hs
    exact (show fiveCycle.Adj (c w) (c v) by simpa using hs).symm
  · exact cycle_bag_edges G htriangle φ ((hc v).2 hv) ((hc w).2 hw) hvw

theorem cycle_support_large_of_neighbor_bound
    (htriangle : ∀ a b c : V, ¬ (G.Adj a b ∧ G.Adj b c ∧ G.Adj c a))
    (φ : CycleVertex → V) (hφ : IsCycleCopy G φ) (x : V)
    (hlarge : ∀ j, (cycleBag G φ j).card <
      (Finset.univ.filter (fun v => G.Adj x v ∧ v ∉ cycleRemainder G φ)).card) :
    2 ≤ (cycleSupport G φ x).card := by
  classical
  by_contra h
  have hfinite : ∀ s : Finset CycleVertex, s.card < 2 → ∃ j, ∀ k ∈ s, k = j := by decide
  obtain ⟨j, hj⟩ := hfinite (cycleSupport G φ x) (by omega)
  have hsub : Finset.univ.filter (fun v => G.Adj x v ∧ v ∉ cycleRemainder G φ) ⊆ cycleBag G φ j := by
    intro v hv
    obtain ⟨hxv, hvr⟩ := (Finset.mem_filter.mp hv).2
    obtain ⟨k, hk⟩ := cycle_partition_covers G htriangle φ hφ v hvr
    have hkj := hj k (Finset.mem_filter.mpr ⟨Finset.mem_univ _, v, hk, hxv⟩)
    rwa [hkj] at hk
  have hh := Finset.card_le_card hsub
  have := hlarge j
  omega

end TriangleProof
