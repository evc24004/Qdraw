import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here.parent / "ibm"))
sys.path.insert(0, str(here))
import render
import refit


jobs = json.loads((here / "jobs.json").read_text())
states = {
    job_id: json.loads((here / "states" / f"{job_id}.json").read_text())
    for job_id in jobs
}

_, _, components = render.scene("HUSKY4")
targets = {
    "single-component test (muzzle)": components[4],
    "component 1: head |0>": components[0],
    "component 2: head |1>": components[1],
    "component 3: left ear": components[2],
    "component 4: right ear": components[3],
    "component 5: muzzle": components[4],
    "state-prep pilot: muzzle": components[4],
    "state-prep pilot: muzzle, gate twirl": components[4],
    "state-prep 1: head |0>": components[0],
    "state-prep 2: head |1>": components[1],
    "state-prep 3: left ear": components[2],
    "state-prep 4: right ear": components[3],
    "state-prep 5: muzzle": components[4],
}


def is_state_prep(label):
    return label.startswith("state-prep")


def is_pilot(label):
    return "test" in label or "pilot" in label


def short_label(label):
    if label == "single-component test (muzzle)":
        return "test\n(no DD)"
    if label == "state-prep pilot: muzzle, gate twirl":
        return "sp pilot\n(twirl)"
    if label == "state-prep pilot: muzzle":
        return "sp pilot"
    part = label.split(": ")[1]
    prefix = "sp\n" if is_state_prep(label) else ""
    return prefix + part


def pauli_matrix(key):
    matrix = np.array([[1.0]])
    for q in reversed(range(len(key))):
        matrix = np.kron(matrix, refit.PAULI[key[q]])
    return matrix


def fidelity_with_shot_error(records, psi, n=4):
    """Linear-inversion fidelity and propagated multinomial shot error."""
    key_counts = {}
    record_terms = []
    for record in records:
        basis = [refit.BASIS[index] for index in record["m_idx"]]
        terms = []
        for mask in range(1, 2**n):
            key = tuple(
                basis[q] if (mask >> q) & 1 else "I" for q in range(n)
            )
            key_counts[key] = key_counts.get(key, 0) + 1
            terms.append((mask, key))
        record_terms.append(terms)

    d = 2**n
    coefficients = {
        key: float(np.real(psi.conj() @ pauli_matrix(key) @ psi))
        / (d * count)
        for key, count in key_counts.items()
    }

    fidelity = 1 / d
    variance = 0.0
    for record, terms in zip(records, record_terms):
        counts = record["counts"]
        shots = sum(counts.values())
        values = []
        weights = []
        for bits, count in counts.items():
            value = 0.0
            for mask, key in terms:
                parity = sum(
                    int(bits[len(bits) - 1 - q])
                    for q in range(n)
                    if (mask >> q) & 1
                )
                value += coefficients[key] * (-1) ** (parity % 2)
            values.append(value)
            weights.append(count)
        values = np.asarray(values)
        weights = np.asarray(weights)
        mean = float(np.dot(values, weights) / shots)
        fidelity += mean
        if shots > 1:
            variance += float(
                np.dot(weights, (values - mean) ** 2) / (shots * (shots - 1))
            )
    return fidelity, np.sqrt(variance)


labels = [short_label(metadata["label"]) for metadata in jobs.values()]
fidelities = []
errors = []
colors = []
for job_id, metadata in jobs.items():
    fock, unitary, _ = targets[metadata["label"]]
    psi = unitary[:, fock]
    records = json.loads((here / "counts" / f"{job_id}.json").read_text())
    fidelity, standard_error = fidelity_with_shot_error(records, psi)
    fidelities.append(fidelity)
    errors.append(1.96 * standard_error)
    if is_pilot(metadata["label"]):
        colors.append("#777777")
    elif is_state_prep(metadata["label"]):
        colors.append("#ff7f0e")
    else:
        colors.append("#1f77b4")

figure, axis = plt.subplots(figsize=(11, 4.8))
bars = axis.bar(
    labels,
    fidelities,
    yerr=errors,
    capsize=3,
    color=colors,
    edgecolor="none",
)
axis.axhline(1 / 16, color="#b22222", linestyle="--", linewidth=1)
axis.text(len(labels) - 0.55, 1 / 16 + 0.012,
          "fully mixed state (1/16)", ha="right", color="#b22222")
axis.set_ylabel("fidelity to the prepared target state")
axis.set_ylim(0, 1.0)
axis.legend(
    handles=[
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#1f77b4", label="unitary circuits (230-245 CZ)"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#ff7f0e", label="state-preparation circuits (17 CZ)"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#777777", label="single-component pilots"),
    ],
    loc="upper left",
    frameon=False,
    fontsize=8,
)
axis.set_title("Fidelity reconstructed from the archived tomography counts")
for bar, value in zip(bars, fidelities):
    axis.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.014,
        f"{value:.3f}",
        ha="center",
        va="bottom",
        fontsize=8,
    )
figure.text(
    0.5,
    0.01,
    "Linear inversion; error bars are 95% shot-noise intervals. Systematic hardware errors are not included.",
    ha="center",
    fontsize=8,
)
figure.tight_layout(rect=(0, 0.04, 1, 1))
figure.savefig(here / "fidelity.png", dpi=160)


d = 16
annihilation = np.diag(np.sqrt(np.arange(1, d)), 1)
q_operator = (annihilation.conj().T + annihilation) / np.sqrt(2)
p_operator = 1j * (annihilation.conj().T - annihilation) / np.sqrt(2)

figure, axis = plt.subplots(figsize=(6, 6))
label_offsets = {
    "head |0>": (10, 9),
    "head |1>": (-48, -15),
    "left ear": (8, 8),
    "right ear": (8, 8),
    "muzzle": (8, -22),
}
for job_id, metadata in jobs.items():
    if is_pilot(metadata["label"]):
        continue
    color = "#ff7f0e" if is_state_prep(metadata["label"]) else "#1f77b4"
    state = states[job_id]
    rho = np.array(state["rho_real"]) + 1j * np.array(state["rho_imag"])
    fock, unitary, _ = targets[metadata["label"]]
    psi = unitary[:, fock]
    target_q = float(np.real(psi.conj() @ q_operator @ psi))
    target_p = float(np.real(psi.conj() @ p_operator @ psi))
    measured_q = float(np.real(np.trace(rho @ q_operator)))
    measured_p = float(np.real(np.trace(rho @ p_operator)))
    axis.annotate(
        "",
        xy=(measured_q, measured_p),
        xytext=(target_q, target_p),
        arrowprops={"arrowstyle": "->", "color": color, "lw": 1.4},
    )
    axis.plot(target_q, target_p, "o", mfc="none", mec="#222222", ms=8)
    axis.plot(measured_q, measured_p, "o", color=color, ms=6)
    if not is_state_prep(metadata["label"]):
        label = metadata["label"].split(": ")[1]
        axis.annotate(
            label,
            xy=(measured_q, measured_p),
            xytext=label_offsets[label],
            textcoords="offset points",
            fontsize=8,
        )

axis.plot(0, 0, "+", color="#b22222", ms=12)
axis.text(0.1, -0.12, "vacuum", fontsize=8, color="#b22222")
axis.set_xlim(-2.6, 2.6)
axis.set_ylim(-2.6, 2.6)
axis.set_xlabel("q expectation value")
axis.set_ylabel("p expectation value")
axis.set_aspect("equal")
axis.grid(color="#dddddd", linewidth=0.6)
axis.set_title("Target and reconstructed component centers")
axis.legend(
    handles=[
        Line2D([0], [0], marker="o", color="none", markeredgecolor="#222222", markerfacecolor="none", label="target"),
        Line2D([0], [0], marker="o", color="none", markeredgecolor="#1f77b4", markerfacecolor="#1f77b4", label="unitary circuits"),
        Line2D([0], [0], marker="o", color="none", markeredgecolor="#ff7f0e", markerfacecolor="#ff7f0e", label="state-preparation circuits"),
    ],
    loc="lower right",
    frameon=False,
)
figure.text(
    0.5,
    0.01,
    "Arrows show displacement in expectation value; they do not identify a specific noise mechanism.",
    ha="center",
    fontsize=8,
)
figure.tight_layout(rect=(0, 0.04, 1, 1))
figure.savefig(here / "displacement_decay.png", dpi=160)

print("wrote fidelity.png and displacement_decay.png")
