# CP-SAT formulation and completeness argument

Let V = {0,...,n−1}. For each unordered pair {i,j}, a Boolean x_ij indicates
an edge. The fixed-edge enumeration imposes:

1. **Triangle-free:** x_ij + x_ik + x_jk ≤ 2 for every triple.
2. **Degree:** d_i = Σ_(j≠i) x_ij ≤ Δ.
3. **Size:** Σ_(i<j) x_ij = e.
4. **Independence bound α ≤ a (optional):** for every (a+1)-subset S,
   Σ_({i,j}⊂S) x_ij ≥ 1.
5. **Factor-criticality (optional, odd n):** for each deleted vertex r,
   introduce Boolean matching variables y^r_ij for pairs avoiding r;
   impose y^r_ij ≤ x_ij and Σ_(j≠i,r) y^r_ij = 1 for every i≠r.

Matching variables are existential witnesses, not part of a graph's identity.
At n=17 the matching number is at most 8 automatically; factor-criticality
makes it exactly 8. This is a fixed-order graph model.

## Symmetry and exclusions

Vertices are ordered by nonincreasing degree. Within equal-degree classes,
optional row comparisons require the adjacency vector to be no larger than
its image under each vertex transposition, in lexicographic order. A globally
lexicographically least upper-triangular edge vector in each such orbit
satisfies every comparison, so some representative survives. Row comparisons
omit the two swapped vertices and list the other vertices in increasing order.

After finding an edge set E, impose Σ_(e∈E) x_e ≤ |E|−1. Because total size is
fixed, this excludes exactly that graph assignment, including every choice of
matching witnesses. Extra cuts exclude explicit isomorphic copies generated
within equal-degree classes. Nauty identifies repeated isomorphism classes.
These cuts remove only classes already represented.

## Covering all degree cases

For n=17 and Δ=7, total deficiency Σ_i(7−d_i)=119−2e is 3 at e=58 and 5 at
e=57. All integer partitions, padded with zeros, give three and seven degree
sequences. Every case is included and finishes with residual INFEASIBLE after
exclusions. An OPTIMAL feasibility response alone indicates one solution,
not completion.

## Second encoding

`audit_enumeration.py` uses nonedge variables, opposite lexicographic orientation,
no factor-critical matching constraints, and only individual assignment
exclusions. It independently generates degree sequences and compares classes
using NetworkX isomorphism. All ten cases terminate with INFEASIBLE and agree
with the four-class catalogue; no non-factor-critical graphs are found.

Both encodings use CP-SAT. Logs, source hashes, tests against small brute-force
catalogues, and the second encoding support a computational result; they are
not a formally checked infeasibility certificate.
