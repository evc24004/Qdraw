import numpy as np
from scipy.linalg import expm

GLYPHS = ["HUSKY4"]
BACKEND = "aer"
SHOTS = 4096
RES = 128


def mode_ops(d):
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    n = np.diag(np.arange(d))
    D = lambda al: expm(al * a.conj().T - np.conj(al) * a)
    R = lambda th: expm(1j * th * n)
    S = lambda z: expm((np.conj(z) * a @ a - z * a.conj().T @ a.conj().T) / 2)
    return D, R, S


def scene(name):
    if name == "HUSKY4":
        nq = 4
        D, R, S = mode_ops(2 ** nq)
        al = lambda q, p: (q + 1j * p) / np.sqrt(2)
        vbar = lambda q, p, s: D(al(q, p)) @ S(np.log(s))
        hbar = lambda q, p, s: D(al(q, p)) @ S(-np.log(s))
        ear = lambda q, p, deg: D(al(q, p)) @ R(np.radians(deg)) @ S(np.log(1.5))
        comps = [
            (0, D(al(0, 0.25)), 0.16),
            (1, D(al(0, 0.25)), 0.22),
            (0, ear(-1.2, 2.0, 15), 0.11),
            (0, ear(1.2, 2.0, -15), 0.11),
            (0, vbar(0, -1.0, 1.25), 0.16),
            (0, hbar(0, -1.85, 1.25), 0.08),
            (0, D(al(-1.35, -0.5)), 0.08),
            (0, D(al(1.35, -0.5)), 0.08),
        ]
        return nq, 3.6, comps

    nq = 5
    D, R, S = mode_ops(2 ** nq)
    al = lambda q, p: (q + 1j * p) / np.sqrt(2)
    vbar = lambda q, p, s: D(al(q, p)) @ S(np.log(s))
    hbar = lambda q, p, s: D(al(q, p)) @ S(-np.log(s))
    dbar = lambda q, p, s, deg: D(al(q, p)) @ R(np.radians(deg)) @ S(np.log(s))
    I = np.eye(2 ** nq)
    scenes = {
        "U": [(0, vbar(-1.4, 0.7, 2.1), 1), (0, vbar(1.4, 0.7, 2.1), 1),
              (0, hbar(0, -1.9, 2.0), 1)],
        "C": [(0, hbar(0.3, 2.1, 1.8), 1), (0, vbar(-1.4, 0, 2.1), 1),
              (0, hbar(0.3, -2.1, 1.8), 1)],
        "O": [(2, I, 0.45), (3, I, 0.55)],
        "N": [(0, vbar(-1.5, 0, 2.2), 1), (0, vbar(1.5, 0, 2.2), 1),
              (0, dbar(0, 0, 2.9, 44), 1)],
        "HUSKY": [
            (0, D(al(0, 0.3)), 0.16), (1, D(al(0, 0.3)), 0.22),
            (0, dbar(-1.5, 2.5, 1.6, 15), 0.11), (0, dbar(1.5, 2.5, 1.6, -15), 0.11),
            (0, vbar(0, -1.2, 1.3), 0.16), (0, hbar(0, -2.3, 1.3), 0.08),
            (0, D(al(-1.7, -0.6)), 0.08), (0, D(al(1.7, -0.6)), 0.08)],
    }
    scenes["N2"] = scenes["N"]
    return nq, 4.5, scenes[name]


def component_circuit(nq, fock, unitary):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import UnitaryGate

    qc = QuantumCircuit(nq)
    for q in range(nq):
        if (fock >> q) & 1:
            qc.x(q)
    qc.append(UnitaryGate(unitary), range(nq))
    return qc


def tomography(qc, backend, shots):
    from qiskit import transpile
    from qiskit_experiments.library import StateTomography

    kwargs = {"optimization_level": 3}
    if backend.name == "aer_simulator":
        kwargs["basis_gates"] = ["rz", "sx", "x", "cx"]
    tqc = transpile(qc, backend, **kwargs)
    data = StateTomography(tqc).run(backend, shots=shots)
    data.block_for_results()
    return np.asarray(data.analysis_results("state").value.data)


def husimi(rho, lim):
    d = rho.shape[0]
    qv = np.linspace(-lim, lim, RES)
    Q, P = np.meshgrid(qv, qv)
    alpha = (Q.ravel() + 1j * P.ravel()) / np.sqrt(2)
    ns = np.arange(d)
    logfact = np.cumsum(np.log(np.maximum(ns, 1)))
    r = np.abs(alpha)[:, None]
    with np.errstate(divide="ignore"):
        logmag = np.where(r > 0, ns * np.log(np.where(r > 0, r, 1)), 0.0)
    logmag = np.where((r == 0) & (ns > 0), -np.inf, logmag)
    C = np.exp(-r ** 2 / 2 + logmag - logfact / 2)
    C = C * np.exp(1j * ns * np.angle(alpha)[:, None])
    H = np.einsum("pi,ij,pj->p", C.conj(), rho, C).real
    return np.clip(H, 0, None).reshape(RES, RES)


def get_backend():
    if BACKEND == "aer":
        from qiskit_aer import AerSimulator
        return AerSimulator()
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService()
    return service.least_busy(operational=True)


def main():
    backend = get_backend()
    panels = {}
    for name in GLYPHS:
        nq, lim, comps = scene(name)
        H = np.zeros((RES, RES))
        for i, (fock, U, w) in enumerate(comps):
            print(f"{name}: component {i + 1}/{len(comps)}")
            rho = tomography(component_circuit(nq, fock, U), backend, SHOTS)
            H += w * husimi(rho, lim)
        panels[name] = H
        np.save(f"{name}.npy", H)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4))
    for ax, (name, H) in zip(np.atleast_1d(axes), panels.items()):
        ax.imshow(H, origin="lower", cmap="inferno")
        ax.set_title(name)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig("render.png", dpi=150, facecolor="black")
    print("wrote render.png")


if __name__ == "__main__":
    main()
