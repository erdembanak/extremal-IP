# Individual graphs and cones

This directory's four graphs use the independence-number restriction.
The larger 23-graph fixed-order catalogue is in [`unrestricted/`](unrestricted/).

`index.json` lists four candidates and their cones. Labels are zero-based.
Independent bags: `[0,1,2,3]`, `[4,5,6,7]`, `[8,9,10]`, `[11,12,13]`,
`[14,15,16]`. The blowup joins consecutive bags cyclically.

| Candidate | Removed edge | Cone edges |
|---|---|---|
| F58 | none | 95 |
| F57_44 | (0,4) | 96 |
| F57_43 | (4,8) | 96 |
| F57_33 | (8,11) | 96 |

Candidate JSON uses `num_vertices` and `edges`. Cone JSON uses the downstream
solver's `num_vertices`, `edges_part1`, and `edges_part2` fields, with apex 17.
The two edge fields are **input storage, not a planar partition certificate**.
All cones are K1 joined to the complement of their candidate.

`SHA256SUMS` covers the JSON files. Portable graph6 files are in
[`../runs/full/`](../runs/full/): `candidates.g6` and `cones.g6`, in table order.

Regenerate with `python cp_sat/export_interchange.py`; check using
`python cp_sat/verify_catalogue.py`. The checker verifies the explicit data
and all 58 single-edge deletions, but does not replace exhaustive search.
