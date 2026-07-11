import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here.parent / "ibm"))
import render

PAULI = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]),
         "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.array([[1, 0], [0, -1]])}

jobs = json.loads((here / "jobs.json").read_text())
states = {j: json.loads((here / "states" / f"{j}.json").read_text()) for j in jobs}

nq, lim, comps = render.scene("HUSKY4")
targets = {
    "single-component test (muzzle)": comps[4],
    "component 1: head |0>": comps[0],
    "component 2: head |1>": comps[1],
    "component 3: left ear": comps[2],
    "component 4: right ear": comps[3],
    "component 5: muzzle": comps[4],
}

def calibration_prediction():
    cal = json.loads((here / "calibration_ibm_kingston.json").read_text())
    line = [151, 150, 149, 148]
    edges = [{151, 150}, {150, 149}, {149, 148}]
    cz_errs, sx_errs = [], []
    for g in cal["gates"]:
        err = next((p["value"] for p in g["parameters"]
                    if p["name"] == "gate_error"), None)
        if err is None:
            continue
        if g["gate"] == "cz" and set(g["qubits"]) in edges:
            cz_errs.append(err)
        elif g["gate"] == "sx" and g["qubits"][0] in line:
            sx_errs.append(err)
    ro_errs = []
    for q in line:
        for p in cal["qubits"][q]:
            if p["name"] == "readout_error":
                ro_errs.append(p["value"])
    e_cz, e_sx = np.mean(cz_errs), np.mean(sx_errs)
    preds = []
    for meta in jobs.values():
        ops = meta["sample_circuit"]["ops"]
        n1 = ops.get("sx", 0) + ops.get("x", 0)
        preds.append((1 - e_cz) ** ops["cz"] * (1 - e_sx) ** n1
                     * np.prod([1 - r for r in ro_errs]))
    return float(np.mean(preds))

labels = ["test\n(no DD)", "head |0>", "head |1>", "ear L", "ear R", "muzzle"]
fids = [states[j]["fidelity_vs_ideal"] for j in jobs]

fig, ax = plt.subplots(figsize=(7, 4.2))
colors = ["#888888"] + ["#1f77b4"] * 5
ax.bar(labels, fids, color=colors)
pred = calibration_prediction()
print(f"calibration-based prediction: {pred:.3f}")
ax.axhline(pred, ls="--", c="green", lw=1)
ax.text(5.45, pred + 0.01,
        f"predicted from archived calibration ({pred:.2f}):\n"
        "per-gate errors on qubits 151-148, no decoherence",
        ha="right", va="bottom", fontsize=8, color="green")
ax.axhline(1 / 16, ls="--", c="red", lw=1)
ax.text(5.45, 0.075, "fully mixed state", ha="right", va="bottom",
        fontsize=8, color="red")
ax.set_ylabel("fidelity to target state")
ax.set_ylim(0, 1)
ax.set_title("ibm_kingston: measured fidelity vs calibration-sheet prediction")
fig.tight_layout()
fig.savefig(here / "fidelity.png", dpi=140)

d = 16
a = np.diag(np.sqrt(np.arange(1, d)), 1)
qop = (a.conj().T + a) / np.sqrt(2)
pop = 1j * (a.conj().T - a) / np.sqrt(2)

fig, ax = plt.subplots(figsize=(5.6, 5.6))
for j, meta in jobs.items():
    if "test" in meta["label"]:
        continue
    st = states[j]
    rho = np.array(st["rho_real"]) + 1j * np.array(st["rho_imag"])
    fock, U, w0 = targets[meta["label"]]
    psi = U[:, fock]
    qd = float(np.real(psi.conj() @ qop @ psi))
    pd_ = float(np.real(psi.conj() @ pop @ psi))
    qm = float(np.real(np.trace(rho @ qop)))
    pm = float(np.real(np.trace(rho @ pop)))
    ax.annotate("", xy=(qm, pm), xytext=(qd, pd_),
                arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.4))
    ax.plot(qd, pd_, "o", mfc="none", mec="black", ms=8)
    ax.plot(qm, pm, "o", color="#1f77b4", ms=6)
    ax.text(qd + 0.1, pd_ + 0.1, meta["label"].split(": ")[1], fontsize=8)
ax.plot(0, 0, "r+", ms=12)
ax.text(0.1, -0.35, "vacuum", fontsize=8, color="red")
ax.set_xlim(-2.6, 2.6)
ax.set_ylim(-2.6, 2.6)
ax.set_xlabel("q")
ax.set_ylabel("p")
ax.set_aspect("equal")
ax.set_title("where each state was placed (open) vs where it was\n"
             "measured (filled): decay drags everything toward vacuum")
fig.tight_layout()
fig.savefig(here / "displacement_decay.png", dpi=140)
print("wrote fidelity.png and displacement_decay.png")
