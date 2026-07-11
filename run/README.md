# Hardware run provenance

Everything needed to check the ibm_kingston run without taking my word
for it. IBM has no public per-job pages (jobs are only visible from the
account that submitted them), so the raw artifacts are published here
instead.

| File | What it is |
|---|---|
| `jobs.json` | The six job IDs with backend, creation time, IBM's own execution timestamps, billed usage (14 s each), shot count, and gate stats of the submitted circuits |
| `counts/<job>.json` | Raw measurement outcomes for all 81 tomography circuits per job: measurement-basis indices (`m_idx`, 0=Z 1=X 2=Y per qubit) and the bitstring counts exactly as returned by the backend |
| `circuits/<job>_pub0.qasm` | One full transpiled ISA circuit per job (OpenQASM 3, 156-qubit layout as actually submitted) |
| `calibration_ibm_kingston.json` | Backend calibration snapshot at run time |
| `requirements-lock.txt` | Exact package versions (`pip freeze`) of the environment that submitted the jobs |
| `states/<job>.json` | Density matrices reconstructed from the published counts, with fidelity and purity |
| `fetch.py` | The script that pulled all of the above from IBM (needs the submitting account) |
| `refit.py` | Independent verification: linear-inversion tomography from `counts/` alone, no IBM account needed |

To verify the fidelity numbers yourself:

```
cd ibm && python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && cd ..
ibm/.venv/bin/python run/refit.py
```

Output from my machine:

```
job                    label                            fidelity  purity
d98ra8af47jc73a896ng   single-component test (muzzle)      0.267   0.163
d98rg3if47jc73a89ct0   component 1: head |0>               0.309   0.169
d98rgcgtcv6s73dmgbhg   component 2: head |1>               0.322   0.171
d98rglgtcv6s73dmgbsg   component 3: left ear               0.306   0.165
d98rgugtcv6s73dmgc60   component 4: right ear              0.318   0.162
d98rh74qp3as739tajjg   component 5: muzzle                 0.315   0.164
```

Two notes on the numbers. The refit here is plain linear inversion, so
it differs by ~0.005 from the PSD-projected fitter qiskit-experiments
uses during a live run (0.267 vs the 0.263 quoted in the main README for
the test job). And the test job was submitted without dynamical
decoupling while the five render jobs had it enabled, which shows up as
the fidelity gap between 0.267 and 0.31-0.32 on otherwise comparable
circuits.
