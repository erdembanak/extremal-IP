Initial validation, 19 September 2026

The implementation ran with OR-Tools 9.15.6755, NetworkX 3.4.2, and
pynauty 2.8.8.1 in the workspace-local virtual environment.

All eleven tests passed in 1.512 seconds. In addition to exhaustive small-graph
comparisons, the model was tested with all graph edges fixed to explicit
17-vertex controls:

- The C5 blowup (4,4,3,3,3) is accepted with 58 edges, α≤7, and
  factor-criticality required.
- B_(7,8) is rejected when α≤7 is required.
- B_(7,8) is accepted when that independence restriction is omitted.

An unseeded solver pilot used n=17, Δ≤7, e=58, α≤7, factor-criticality,
degree pattern 0 (fourteen 7s and three 6s), four workers, and at most
128 additional relabeling exclusions per graph. It found three edge
assignments in 13.784 seconds and stopped at the requested three-solution
limit. Nauty identified a single isomorphism class, and an independent
certificate comparison identified it as the known C5 blowup. The run
installed 387 distinct edge-assignment exclusion constraints.

Artifacts: [pattern-0 summary](runs/pilot-58-pattern0/summary.json),
[saved graphs](runs/pilot-58-pattern0/solutions.jsonl), and
[manifest](runs/pilot-58-pattern0/manifest.json).

A separate pilot on degree pattern 2 (sixteen 7s and one 4), with the same
α and factor-criticality restrictions, ended after approximately thirty
seconds with UNKNOWN and no graph. This establishes neither feasibility
nor infeasibility. See its [summary](runs/pilot-58-pattern2/summary.json).

Neither pilot is a completed enumeration. Isomorphic duplicates remain
possible despite the extra cuts; full enumeration performance is unmeasured.
No biplanarity computations were performed.

Subsequent full run: all ten degree cases are now complete, with thirteen
unit tests passing and a second complement-variable encoding agreeing.
The pilot limitations above describe the earlier runs. See
[RESULTS.md](RESULTS.md) for the completed classification.
