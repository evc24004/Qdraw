# Quantum CFS Rendering

Drawing images with quantum circuits, based on the method from
[Schloss & Usui, "General purpose graphical rendering on quantum devices
with composable function systems"](https://www.leioslabs.com/publications/res/QCG_paper.pdf).
I rendered an apple, then went a bit further and did the UConn husky logo
and the word UCONN. There are two implementations: one in MATLAB (Quantum
Computing support package) and one in Python/Qiskit that can run on real
IBM hardware.

<p align="center">
  <img src="output/uconn_quantum.png" width="70%" alt="UCONN">
</p>
<p align="center">
  <img src="output/apple_qcg.png" width="32%" alt="apple">
  <img src="output/husky_quantum_inferno.png" width="32%" alt="husky">
</p>

## The idea

N qubits encode one bosonic mode truncated to 2^N Fock levels. The mode's
phase space (position vs momentum) is the canvas. Fock states are the
primitives: the vacuum is a Gaussian blob, higher Fock states are rings,
and mixing them fills in a circle. The letter O is just a mixture of the
|2> and |3> Fock states, no drawing involved. Everything else is built by
displacing, squeezing and rotating Gaussians, which are all unitary
operations on the mode. Letters take about 3 bars each.

To get the image back out you can't look at the state directly. Each
circuit goes through full state tomography (all 3^N measurement bases,
4096 shots each), the density matrix is reconstructed from the counts,
and the Husimi function of that reconstructed state is the image. Objects
built from several pure states get one circuit per component and the
Husimi functions are summed with weights. This is the same readout the
paper used for their runs on IBM Kingston.

Each glyph gets its own mode because the usable phase space only spans
about 2N in each direction, which is nowhere near enough room for five
letters side by side. The panels are tiled into the final poster
afterwards. Pixel brightness comes from the measured probabilities, then
each panel gets normalized and a colormap (UConn navy or inferno) at the
end.

## MATLAB version

Needs the Quantum Computing support package (tested on R2026a; anything
recent enough to have `unitaryGate` should work).

```
matlab -batch quantum_apple_qcg
matlab -batch qcg_glyphs
matlab -batch qcg_compose
matlab -batch qcg_husky_solo
```

`qcg_glyphs.m` renders the letter and husky panels (edit the `glyphs`
list at the top if you only want some of them, a full run takes around
20 minutes). `qcg_compose.m` assembles the poster. `quantum_apple_qcg.m`
is the 5-qubit apple and `quantum_apple_qcg_hw.m` is a 4-qubit version
sized for real hardware. Outputs land in `output/`.

## IBM version

`ibm/render.py` is a standalone Qiskit implementation of the same thing.
It builds the same circuits, runs qiskit-experiments state tomography and
produces the same Husimi images. Settings are variables at the top of the
file.

```
cd ibm
python3 -m venv .venv          # needs python 3.11 or newer
.venv/bin/pip install -r requirements.txt
.venv/bin/python render.py
```

The default mode (`BACKEND = "rehearse"`) simulates with a rough noise
model made from the backend's reported median error rates. It is a
sanity check, not a digital twin. In practice it was optimistic by about
a factor of 2, see `run/` for the comparison. `"aer"` gives a noiseless
simulation and `"ibm"` submits to real hardware. For that you need an
IBM Quantum account (the free plan gives 10 minutes of QPU time a month)
and your API key saved once:

```
.venv/bin/python -c "from qiskit_ibm_runtime import QiskitRuntimeService; \
  QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token='YOUR_KEY')"
```

Stick to the 4-qubit scene on real devices. A generic 4-qubit unitary
decomposes to ~95 entangling gates, and 230-245 physical CZ gates once
routed onto hardware; the 5-qubit glyphs are several times worse, which
is past the point where anything comes back. The paper made the same
call, their hardware runs were 4 qubits per mode too.

## Run on ibm_kingston

I ran the husky on ibm_kingston, the machine from the paper. Five
tomography jobs, one per mixture component, 81 circuits x 512 shots each
with dynamical decoupling on, placed on the best connected line of four
qubits in that day's calibration data (qubits 151-148). IBM billed 14
seconds per job, 70 for the render plus 14 for an earlier test.

<p align="center">
  <img src="output/husky_quantum_inferno.png" width="38%" alt="matlab simulation">
  <img src="output/husky_kingston.png" width="38%" alt="measured on ibm_kingston">
</p>
<p align="center"><i>left: MATLAB simulation. right: measured on ibm_kingston.</i></p>

The head and muzzle survive, the ears mostly don't. Component fidelities
came out at 0.31-0.32 (0.27 for a test job that ran without dynamical
decoupling). The submitted circuits carry 230-245 CZ gates each, and the
qubits lose a lot of state over a circuit that deep; the loss also moves
displaced states back toward the center of the image, which is
consistent with amplitude damping though other errors contribute. The
paper's hardware figures show the same kind of damage.

The [`run/`](run/) folder has the receipts: job IDs, raw counts for all
486 circuits, representative submitted circuits, archived runtime
options, a calibration snapshot, locked package versions, plots of how
the errors behave, and scripts that recheck the fidelity numbers and
rebuild the husky image from the counts alone.

## References

- J. Schloss and A. Usui, General purpose graphical rendering on quantum
  devices with composable function systems.
  https://www.leioslabs.com/publications/res/QCG_paper.pdf
- J. Schloss, Composable function systems as a general-purpose rendering
  framework. https://www.leioslabs.com/publications/res/CFS_paper.pdf
- MATLAB Support Package for Quantum Computing.
  https://www.mathworks.com/products/quantum-computing.html
