# Unfiltered 17-vertex catalogue

Two classes at 58 edges and 21 at 57 edges, all triangle-free with maximum
degree at most 7. No independence-number or factor-criticality filter was used.
See [`../../UNRESTRICTED.md`](../../UNRESTRICTED.md) for scope and validation.

Labels are 0,...,16. Each JSON has `num_vertices`, `edges`, `properties`,
`graph6`, and a nauty `certificate`. IDs use the edge count and a one-based
index sorted by certificate; they need not match IDs in a fresh run with
different canonicalization versions. Check isomorphism when comparing runs.

`index.json` lists every file and the edge-deletion relationships to the two
58-edge parents. Empty parent maps mean the graph cannot be obtained by one
edge deletion from either 58-edge class. The graph6 file order matches the
index order within each edge count. `SHA256SUMS` covers the JSON and graph6 data.
