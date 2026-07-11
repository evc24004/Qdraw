# Hardware run provenance

Receipts for the ibm_kingston run. IBM doesn't have public pages for
jobs (only the submitting account can query them), so the raw data is
committed here instead and the analysis can be rerun from it directly.

| File | What it is |
|---|---|
| `jobs.json` | the six job IDs, creation times, IBM's execution timestamps, billed usage (14 s each), and gate stats for the submitted circuits |
| `counts/` | raw bitstring counts for all 81 tomography circuits of each job, plus the measurement basis per circuit (`m_idx`, 0=Z 1=X 2=Y) |
| `circuits/` | one full transpiled circuit per job, OpenQASM 3, on the 156-qubit layout as actually submitted |
| `calibration_ibm_kingston.json` | the backend's calibration snapshot from run time |
| `requirements-lock.txt` | pip freeze of the environment that submitted everything |
| `states/` | density matrices reconstructed from the counts, with fidelity and purity per job |
| `fetch.py` | pulls all of the above from IBM (needs the submitting account) |
| `refit.py` | redoes the tomography from `counts/` alone, no account needed |
| `plots.py` | makes the figures below from the reconstructed states |

To check the numbers yourself:

```
cd ibm && python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && cd ..
ibm/.venv/bin/python run/refit.py
```

which prints, from my machine:

```
job                    label                            fidelity  purity
d98ra8af47jc73a896ng   single-component test (muzzle)      0.267   0.163
d98rg3if47jc73a89ct0   component 1: head |0>               0.309   0.169
d98rgcgtcv6s73dmgbhg   component 2: head |1>               0.322   0.171
d98rglgtcv6s73dmgbsg   component 3: left ear               0.306   0.165
d98rgugtcv6s73dmgc60   component 4: right ear              0.318   0.162
d98rh74qp3as739tajjg   component 5: muzzle                 0.315   0.164
```

The refit is plain linear inversion, so it lands ~0.005 away from the
PSD-projected fitter used during the live run (0.267 here vs the 0.263
quoted for the test job). The test job also ran without dynamical
decoupling while the five render jobs had it on, which is worth about
+0.05 fidelity on comparable circuits.

## What the errors look like

<p align="center"><img src="fidelity.png" width="75%"></p>

Kingston's calibration sheet (median gate errors, depolarizing model)
predicts ~0.77 fidelity for these circuits. The machine delivered ~0.31.
The gap is decoherence during the ~360-layer circuits, worst-case qubits,
and correlated errors, none of which show up in per-gate medians.

<p align="center"><img src="displacement_decay.png" width="65%"></p>

Same effect seen in phase space: energy relaxation during the circuit
drags every state back toward the vacuum. The ears, displaced furthest,
lose the most, which is why they nearly vanish in the hardware image
while the head barely moves.
