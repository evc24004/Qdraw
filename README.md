# Quantum CFS rendering

This repository reproduces part of the rendering method in Schloss and Usui's
paper, [General purpose graphical rendering on quantum devices with composable
function systems](https://www.leioslabs.com/publications/res/QCG_paper.pdf).
The paper's method is not mine. The additions here are an apple, a UConn wordmark
and a simplified husky assembled from the same Fock-state primitives and affine
operators.

<p align="center">
  <img src="output/uconn_quantum.png" width="68%" alt="UConn glyph render">
</p>
<p align="center">
  <img src="output/apple_qcg.png" width="32%" alt="Apple render">
  <img src="output/husky_quantum_inferno.png" width="32%" alt="Husky render">
</p>

## Method

An `N`-qubit register represents a harmonic oscillator truncated to `2^N`
Fock levels. Each component starts in a Fock state and is displaced, rotated or
squeezed. Full state tomography reconstructs its density matrix. The image is
the Husimi function of that matrix, and multi-part objects are weighted sums of
component images.

The MATLAB scripts use 4096 shots per tomography basis. The archived Kingston
run used four qubits, 81 bases and 512 shots per basis.

## MATLAB

The scripts require MATLAB's Support Package for Quantum Computing. They were
run with R2026a.

```sh
matlab -batch quantum_apple_qcg
matlab -batch qcg_glyphs
matlab -batch qcg_compose
matlab -batch qcg_husky_solo
```

The first two scripts perform the tomography. `qcg_compose.m` and
`qcg_husky_solo.m` only arrange and color the reconstructed images.

## Qiskit

```sh
python3 -m venv ibm/.venv
ibm/.venv/bin/pip install -r ibm/requirements.txt

# Noiseless Aer simulation; this is the default and needs no IBM account.
ibm/.venv/bin/python ibm/render.py

# Aer with a simple noise model derived from the archived Kingston calibration.
ibm/.venv/bin/python ibm/render.py --backend rehearse

# Submit to an IBM backend. This requires saved IBM Quantum credentials.
ibm/.venv/bin/python ibm/render.py --backend ibm --backend-name ibm_kingston
```

Run `ibm/render.py --help` to set the scene, shot count, resolution, output path
or physical qubits. The IBM mode chooses a connected line with a heuristic that
weights CZ and SX errors by the approximate gate mix of these circuits. It is a
placement heuristic, not an optimizer.

## Archived Kingston run

The hardware render consists of five tomography jobs submitted on July 10,
2026, plus one earlier muzzle test. The representative submitted circuits have
230-245 physical CZ gates and Qiskit depths of 1028-1110. The reconstructed
component fidelities are 0.309-0.322; the earlier test is 0.267.

<p align="center">
  <img src="output/husky_quantum_inferno.png" width="38%" alt="Simulated husky">
  <img src="output/husky_kingston.png" width="38%" alt="Kingston husky">
</p>

The hardware image loses most of the ear detail. The reconstructed component
centers also move toward the vacuum. That is consistent with amplitude damping,
but these data do not separate damping from preparation, correlated and
tomography errors.

The raw counts, one submitted circuit from each job, runtime options,
calibration snapshot and reconstruction scripts are in [`run/`](run/). They can
be checked without an IBM account:

```sh
ibm/.venv/bin/python run/refit.py
ibm/.venv/bin/python run/rebuild_husky.py
ibm/.venv/bin/python run/plots.py
```

The rehearsal backend and the calibration-only fidelity estimate are deliberately
simple. They model independent gate and readout errors but not idle time,
coherence, crosstalk or drift, so neither should be treated as a prediction of a
future hardware run.

## References

- J. Schloss and A. Usui, *General purpose graphical rendering on quantum
  devices with composable function systems*.
- J. Schloss, *Composable function systems as a general-purpose rendering
  framework*.
- [MATLAB Support Package for Quantum Computing](https://www.mathworks.com/products/quantum-computing.html)
