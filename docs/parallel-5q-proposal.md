# Proposal: high-detail hardware husky via 5-qubit components on parallel tomography lines

Status: draft for review. Nothing in this document has been submitted to
hardware. All hardware numbers quoted as "measured" come from the two
archived Kingston runs in `run/`; all other numbers are estimates and are
labeled as such.

## 1. Objective

Render the detailed 8-component husky (the one shown in the simulated
poster) on ibm_kingston with recognizable quality, by combining two
changes to the July 11 approach:

1. 5 qubits per mode instead of 4, enlarging the truncated Fock space
   from 16 to 32 levels and the usable phase space from about +-4 to
   about +-5. This is what the simulated renders already use.
2. Parallel state tomography: place several components on disjoint
   qubit lines of the 156-qubit chip and measure them in a single job,
   instead of one job per component.

Neither change alters the rendering method: each component is still a
prepared state, reconstructed by full state tomography, imaged as a
Husimi function, and combined classically with fixed weights.

## 2. Measured baseline

From the archived runs (see `run/README.md`):

| Run | Circuits | CZ submitted | Depth | Component fidelity |
|---|---|---|---|---|
| July 10, unitary construction, 4 qubits | 5 + 1 test | 230-245 | 1028-1110 | 0.306-0.322 |
| July 11, state preparation, 4 qubits | 5 + 2 pilots | 23-29 | 122-151 | 0.872-0.905 |

Locally measured this session (transpile only, not submitted): the eight
5-qubit husky components cost 48-50 CZ at depth 165-169 when prepared as
states and routed on a 5-qubit line. For comparison, the same components
as full unitaries would cost roughly 500 two-qubit gates, which is why
5-qubit scenes were previously ruled out for hardware.

## 3. Why 5 qubits helps the image

The husky in the simulated poster uses 8 components at 5 qubits: head
disc (two Fock layers), two ears, muzzle, chin and two cheeks. The
4-qubit hardware husky dropped the chin and cheeks and compressed the
geometry to fit the smaller phase space, which is why it looks stubbier
than the simulation even at high fidelity. A 5-qubit run renders the
same scene as the poster.

## 4. Proposed design

### 4.1 Circuits

Identical to the July 11 methodology (`--preparation state`), at 5
qubits: `StatePreparation(U[:, fock])` per component, transpiled to the
backend target. Expected submitted size, by analogy with the 4-qubit
case where routing and tomography rotations added about 10 percent: 50-60
CZ, depth under 250. To be confirmed by `--dry-run` before any
submission.

### 4.2 Tomography

Full state tomography needs 3^5 = 243 basis settings per component
(versus 81 at 4 qubits). At 512 shots per setting this is 124,416 shots
per component, 3x the July runs.

### 4.3 Parallelism

qiskit-experiments provides `ParallelExperiment`, which composes k
experiments on disjoint physical qubits into one set of circuits; the
analysis marginalizes each subsystem independently. With 8 components on
8 disjoint 5-qubit lines (40 measured qubits), one parallel job replaces
8 sequential jobs. The number of circuits per job stays 243 (the basis
settings run in lockstep across lines).

Qubit selection: extend `select_low_error_line` to greedily pick k
disjoint lines: pick the best line, remove its qubits and their
neighbors from the graph, repeat. Kingston has 156 qubits; 8 lines of 5
with one-qubit buffers needs roughly 75 qubits of real estate, which the
heavy-hex lattice accommodates, but the 8th line will be measurably
worse than the 1st. The report from the selector (per-line error
estimates) is part of the pilot go/no-go data.

### 4.4 What parallelism costs

- Crosstalk: simultaneous CZ gates and readout on neighboring lines add
  error that sequential jobs do not see. Not modeled in any of our
  simulators; this is precisely what the hardware pilot must measure.
- Worse qubits: sequential jobs reuse the single best line; parallel
  jobs must use 8 different lines. Expected penalty grows with k.
- Drift immunity is the compensating benefit: all components see the
  same calibration window, which removes the day-apart caveat that
  applies to comparisons in the current archive.

## 5. Cost estimates (not commitments)

Per the July billing pattern (13-14 s per 81-setting job, roughly
proportional to total shots plus fixed overhead):

| Step | Jobs | Estimated QPU |
|---|---|---|
| Dry run (read-only transpile) | 0 | 0 s |
| Pilot A: one 5-qubit component, best line | 1 | 30-40 s |
| Pilot B: two components on two parallel lines | 1 | 30-40 s |
| Full parallel render: 8 components, 8 lines | 1 | 35-50 s |
| Total | 3 | 95-130 s |

Remaining monthly allowance is approximately 5.5 minutes, so the
sequence fits with about 3.5 minutes of reserve. If pilots force a
redesign, the reserve covers one more full attempt. Sequential fallback
(8 jobs at 30-40 s each) does not fit comfortably and is not planned.

## 6. Expected quality, stated carefully

Per-gate arithmetic with July 11 calibration values gives a survival
proxy around 0.80-0.86 for 50-60 CZ circuits. The July 11 experience
suggests such arithmetic can be approximately right when circuits are
shallow (measured 0.87-0.91 against a 0.92 proxy at 23-29 CZ), but the
July 10 experience shows it fails badly for deep circuits. Crosstalk
from parallel operation is an unknown on top. Therefore: the working
hypothesis is component fidelities in the 0.6-0.85 range, and no number
goes into the README until measured.

## 7. Validation ladder before any submission

### Phase 1: MATLAB, noiseless, full scale

Purpose: validate the 5-qubit scene, state-prep circuits and tomography
pipeline end to end in a second implementation, as was done for the
4-qubit method.

1. Add an `initGate`-based state-preparation path to the MATLAB glyph
   engine (initGate is MATLAB's equivalent of StatePreparation; a spot
   check this session measured 22 two-qubit gates for a 4-qubit
   component versus 115 via unitaryGate).
2. Run full 243-basis tomography of all 8 husky components at 4096
   shots in simulation, reconstruct, compose the image.
3. Acceptance: per-component fidelity of the reconstructed state against
   the exact state above 0.99 (shot noise only); composed image visually
   identical to the existing `output/husky_quantum_inferno.png`; total
   variation between measured distribution and exact distribution
   consistent with shot noise.

### Phase 2: Qiskit, ideal Aer

1. Same render through `ibm/render.py` at 5 qubits with
   `--preparation state` (exists already; needs the HUSKY scene wired to
   the ibm path, which currently guards against mixed qubit counts).
2. Parallel correctness: build `ParallelExperiment` with 2, then 4,
   then 8 StateTomography instances on disjoint lines of a 40-qubit ideal
   Aer target. Acceptance: every marginal density matrix equals the
   corresponding single-line reconstruction to within shot noise
   (fidelity between the two reconstructions above 0.99), and the
   composed image matches Phase 1.
3. Memory check: 243 circuits on 40 qubits is the largest simulation in
   the project so far; the run must stay under the 2.5 GB guard used
   throughout. If statevector simulation of 40 qubits is infeasible (it
   is, at 2^40 amplitudes), the parallel validation runs on 2x4-qubit
   and 2x5-qubit configurations (10-18 qubits total), which is
   sufficient to validate the marginalization logic; full 8-line
   validation is structural only (circuit construction, layout, qubit
   accounting).

### Phase 3: Qiskit, rehearsal noise

1. Per-line noisy simulation at 5 qubits with the archived calibration
   noise model (density-matrix method). This bounds the no-crosstalk
   expectation.
2. Acceptance to proceed: rehearsal fidelity above 0.8 per component.
   Below that, the geometry gets re-examined before spending hardware
   time.
3. Explicit limitation: no simulator in this project models crosstalk
   between parallel lines. That risk is retired only by hardware pilot B.

### Phase 4: gated hardware sequence

Each step requires the previous one to pass, and aborts if its result
is more than 0.15 below the rehearsal expectation:

1. `--dry-run` against the live backend: confirm CZ counts and depths
   are in the estimated range, record the 8 selected lines and their
   predicted errors.
2. Pilot A: single 5-qubit component (muzzle) on the best line, one
   job. Confirms the 5-qubit fidelity level in isolation.
3. Pilot B: the same component duplicated on two parallel lines in one
   job. The difference between the two marginals, and between them and
   pilot A, is the direct crosstalk-plus-line-quality measurement.
4. Full render only if pilot B's degradation versus pilot A is under
   0.1 in fidelity. Otherwise fall back to two sequential half-parallel
   jobs (4 components each) and re-budget.

### Phase 5: archive and reporting

Same discipline as the existing runs: extend `run/fetch.py` labels,
`refit.py` needs a 5-qubit path (its linear inversion is already
n-generic; the component map and Pauli dimension change), rebuild script
gains the new image, README numbers come only from the refit.

## 8. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| ParallelExperiment + SamplerV2 + physical_qubits interaction untested by us | medium | Phase 2.2 exercises exactly this path on Aer before hardware |
| Crosstalk between parallel lines | unknown | Pilot B measures it before the full spend |
| 8 disjoint good lines don't exist on that day's calibration | low | Selector reports predicted per-line error; abort threshold before submission |
| 243-setting job exceeds runtime payload or session limits | low | Dry-run constructs the exact job; if too large, split into 3 sub-jobs of 81 settings |
| Tomography refit assumptions (basis order) differ at 5 qubits | low | Phase 2 validates refit against ideal Aer counts at 5 qubits |
| Budget overrun | low | Hard cap: abort the sequence if cumulative new spend exceeds 150 s |

## 9. Open questions for review

1. Is pilot B's design (same component on two lines, one job) the right
   crosstalk probe, or should it duplicate the *busiest* expected
   configuration (8 lines) at reduced shots instead?
2. Is 512 shots per setting still right at 243 settings, or is 256
   acceptable to halve cost? (Shot-noise error bars in `run/plots.py`
   suggest 512 gives about +-0.012 at 95% on fidelity; 256 would give
   about +-0.017.)
3. Should the parallel run keep dynamical decoupling enabled given the
   shorter circuits, or is that worth its own pilot arm? (Proposed:
   keep it on; the July pilots were not conclusive either way and the
   budget favors fewer arms.)
4. The 5-qubit scene has 8 components with weights that sum to 1.0
   exactly (0.16+0.22+0.11+0.11+0.16+0.08+0.08+0.08). Any objection to
   keeping the simulated poster's weights unchanged so hardware and
   simulation stay directly comparable?
5. MATLAB cannot model noise (the support package simulator is
   noiseless), so Phase 1 validates method, not error behavior. Is a
   noiseless second implementation sufficient for the MATLAB phase, or
   should MATLAB also cross-check the *rehearsal* by importing Aer's
   noisy counts and refitting them independently? (Proposed: yes to the
   cross-check; it is cheap and reuses the archive tooling.)

## 10. What is explicitly out of scope

- 6 or more qubits per mode: 729 tomography settings and roughly 120 CZ
  put both cost and fidelity outside the free plan's envelope.
- Treating the whole chip as one canvas: tomography scales as 3^N
  settings; it is a measurement-cost wall, not an engineering gap.
- Readout error mitigation (M3): deliberately deferred; it adds
  calibration circuits and a correction stage whose benefit at these
  readout error levels (0.3-0.7 percent per qubit) is smaller than the
  crosstalk unknown this proposal is designed to measure.

## 11. Revisions after review (accepted)

The reviewer approved the architecture with mandatory changes. All are
adopted; the sections above stand as the original record.

1. The one-job plan is invalid: IBM caps low-level control instructions
   near 26.8M per qubit and the 243x512 job would hit roughly 68M on the
   busiest qubit. The render becomes four jobs of at most 64 settings,
   512 shots each, preserving the eight-way spatial parallelism.
2. Gate estimates were optimistic. Reviewer compilation against the live
   target measured 62-86 CZ per component and child depths of 256-403.
   These numbers replace the 50-60 CZ / depth-250 estimates.
3. Line selection moves from greedy to compiled-cost global selection:
   enumerate five-qubit paths, transpile the real circuit onto
   candidates, score with calibrated errors, and solve the eight-path
   set-packing lexicographically (worst line first, then total error,
   then separation).
4. The two-line pilot is replaced by a 1-2-4-8 concurrency ladder at 128
   shots per setting, same representative component throughout, with
   acceptance criteria: no line below 0.70, median eight-line
   degradation under 0.05 versus lower concurrency, single bad lines
   replaced rather than the design abandoned.
5. Dynamical decoupling is an experimental arm, not a default:
   instruction counts are checked on the scheduled circuits after DD
   insertion, and DD is dropped if it erodes the job-limit margin
   without measured benefit. Twirling stays off.
6. Reconstruction reports both linear inversion (transparent) and a
   PSD/trace-one constrained fit (physical), with bootstrap intervals;
   the constrained fit feeds the image, the raw fit stays in the
   archive.
7. MATLAB's role is scoped to what it can honestly do: independent
   convention and reconstruction checks from exported counts, not noise
   realism.
8. Classical shadows and compressed sensing are deferred until compared
   offline against full tomography on Husimi-pixel error per QPU
   execution.
