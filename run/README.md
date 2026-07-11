# Hardware run provenance

Receipts for the ibm_kingston run. IBM doesn't have public pages for
jobs (only the submitting account can query them), so the raw data is
committed here and the analysis can be rerun from it directly.

| File | What it is |
|---|---|
| `jobs.json` | the six job IDs, creation times, IBM's execution timestamps, billed usage (14 s each), archived runtime options, gate stats, and the measured physical qubits |
| `counts/` | raw bitstring counts for all 81 tomography circuits of each job, plus the measurement basis per circuit (`m_idx`, 0=Z 1=X 2=Y) |
| `circuits/` | one representative submitted circuit per job (pub 0 of 81), OpenQASM 3, on the 156-qubit layout. The other 80 differ only in measurement-basis rotations, which `counts/` records |
| `calibration_ibm_kingston.json` | the backend's calibration snapshot from run time |
| `requirements-lock.txt` | pip freeze of the environment that submitted everything |
| `states/` | density matrices reconstructed from the counts, with fidelity and purity per job |
| `fetch.py` | pulls all of the above from IBM (needs the submitting account) |
| `refit.py` | redoes the tomography from `counts/` alone, no account needed |
| `rebuild_husky.py` | rebuilds the final husky image from `counts/` and prints its sha256 |
| `plots.py` | makes the figures below, including the calibration-based prediction |

To check everything yourself:

```
python3 -m venv .venv
.venv/bin/pip install -r run/requirements-lock.txt
.venv/bin/python run/refit.py
.venv/bin/python run/rebuild_husky.py
```

`refit.py` prints, from my machine:

```
job                    label                            fidelity  purity
d98ra8af47jc73a896ng   single-component test (muzzle)      0.267   0.163
d98rg3if47jc73a89ct0   component 1: head |0>               0.309   0.169
d98rgcgtcv6s73dmgbhg   component 2: head |1>               0.322   0.171
d98rglgtcv6s73dmgbsg   component 3: left ear               0.306   0.165
d98rgugtcv6s73dmgc60   component 4: right ear              0.318   0.162
d98rh74qp3as739tajjg   component 5: muzzle                 0.315   0.164
```

`rebuild_husky.py` closes the chain from raw counts to the final image:
it reconstructs the five component states, applies the scene weights,
sums the Husimi functions and writes `husky_rebuilt.png`
(sha256 `ef5bf914...`). It is visually identical to the published
`output/husky_kingston.png` but not byte-identical, because the live run
fitted the states with qiskit-experiments' PSD-projected fitter while
the refit here uses plain linear inversion. The same difference shows in
the numbers, about 0.005 in fidelity (0.267 here vs 0.263 quoted for
the test job).

Gate counts, stated precisely because three different numbers float
around: each component's logical unitary decomposes to 94-95 entangling
gates in a cz basis before routing. The submitted circuits, after
routing onto the physical qubits and appending measurement rotations,
contain 230-245 physical CZ gates at Qiskit depths of 1028-1110. The
`physical_qubits` field in `jobs.json` lists which qubit each classical
bit measured in the archived circuit; routing settles each circuit on
its own final mapping within qubits 148-151.

The archived `runtime_options` show `dynamical_decoupling: {enable:
true}` on the five render jobs and no DD on the test job. The test job
measured 0.267 on the same muzzle circuit that scored 0.315 with DD.
That is one paired comparison, not a controlled benchmark, so treat the
+0.05 as indicative.

## What the errors look like

<p align="center"><img src="fidelity.png" width="75%"></p>

The dashed line is computed by `plots.py` from the archived calibration:
take the cz, sx and readout errors for qubits 148-151, raise (1 - error)
to the power of each submitted circuit's gate counts, and multiply. That
static per-gate model predicts 0.46. The machine delivered ~0.31. The
missing third is dominated by effects a per-gate model doesn't see,
mainly decoherence while the ~1000-deep circuit executes.

<p align="center"><img src="displacement_decay.png" width="65%"></p>

The same loss seen in phase space: every reconstructed state sits closer
to the vacuum than where it was prepared, and the ears, displaced
furthest, lose the most. The direction and ordering are consistent with
amplitude damping during the circuit, though readout bias, preparation
error and tomography artifacts will contribute too.
