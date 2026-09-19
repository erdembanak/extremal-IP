# Reproducing the catalogue

The repository includes completed runs in `cp_sat/runs/full/`. Install
`cp_sat/requirements.txt` in a virtual environment (Python 3.10 was used).
Run these commands from the repository root:

```sh
.venv/bin/python -m unittest discover -s cp_sat -v
.venv/bin/python cp_sat/verify_catalogue.py
```

To repeat the search rather than read completed checkpoints, preserve the
supplied records first. In a fresh checkout, run once:

```sh
mv cp_sat/runs/full cp_sat/runs/published
.venv/bin/python cp_sat/full_enumeration.py --round-seconds 120
.venv/bin/python cp_sat/audit_enumeration.py
.venv/bin/python cp_sat/export_catalogue.py
.venv/bin/python cp_sat/export_interchange.py
.venv/bin/python cp_sat/verify_catalogue.py
```

The main run uses three jobs with three solver workers each and resumes
timed-out cases until complete. The audit has a 90-second per-case limit and
must report `complete_agreement: true`. A timeout is not infeasibility.
Representative labelings can differ between runs; compare isomorphism classes.

`runs/full/status.json` must report all ten cases complete and counts 1 at 58
edges and 3 at 57 edges. Every case must finish with residual INFEASIBLE.
`runs/full/audit.json` records the second encoding and source hash. Production
and audit sources are distributed unchanged to preserve their recorded hashes.

## Downstream solver input

In the separately installed `earth-moon-alpha2` repository, a cone file can be
used with:

```sh
python verifier/biplanar_sat_prop.py --json /path/to/extremal-IP/cp_sat/data/cone_F58.json --no-sym
```

The inspected downstream version imports a missing `c5_inflation_attack`
helper; that import must be repaired first, although the JSON route never
calls the helper. Cone inputs store all edges in `edges_part1` and leave
`edges_part2` empty as a storage convention, **not a planar partition claim**.

Only the F58 cone needs testing to exclude all four candidates: each other
cone contains it. This enumeration does not verify non-biplanarity.
