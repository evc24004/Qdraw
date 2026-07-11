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
afterwards, and the colormaps (UConn navy or inferno) are applied at the
end. Brightness is measured probability throughout.

## MATLAB version

Needs the Quantum Computing support package (R2023b or newer).

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
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python render.py
```

By default it renders the 4-qubit husky with a noise model built from the
target backend's live calibration (`BACKEND = "rehearse"`). Set it to
`"aer"` for an ideal simulation or `"ibm"` to submit to real hardware
(needs saved IBM Quantum credentials, see their docs).

This actually ran on ibm_kingston, the same machine used in the paper.
Five tomography jobs (one per mixture component, 81 circuits x 512 shots
each, dynamical decoupling on) on its best-calibrated 4-qubit line:

<p align="center">
  <img src="output/husky_kingston.png" width="40%" alt="husky measured on ibm_kingston">
</p>

The head and muzzle survive; the ears mostly wash out. Measured component
fidelity was around 0.26 against the ideal states, dominated by
decoherence over the ~150 two-qubit gates each state preparation needs,
which also drags displaced states back toward the origin. That is the
honest state of NISQ hardware, and the same degradation is visible in the
paper's own hardware figures. A word of warning before spending
QPU time: a generic 4-qubit unitary transpiles to roughly 95 two-qubit
gates, which survives current devices reasonably well. The 5-qubit
glyphs are around 5x more expensive and mostly come back as noise, so
those are better left on the simulator. The paper made the same
compromise, their hardware runs were 4 qubits per mode too.

## References

- J. Schloss and A. Usui, General purpose graphical rendering on quantum
  devices with composable function systems.
  https://www.leioslabs.com/publications/res/QCG_paper.pdf
- J. Schloss, Composable function systems as a general-purpose rendering
  framework. https://www.leioslabs.com/publications/res/CFS_paper.pdf
- MATLAB Support Package for Quantum Computing.
  https://www.mathworks.com/products/quantum-computing.html
