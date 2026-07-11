import numpy as np
from scipy.linalg import expm

GLYPHS = ["HUSKY4"]
BACKEND = "rehearse"
SHOTS = 512
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
        ear = lambda q, p, deg: D(al(q, p)) @ R(np.radians(deg)) @ S(np.log(1.4))
        comps = [
            (0, D(al(0, 0.1)), 0.16),
            (1, D(al(0, 0.1)), 0.22),
            (0, ear(-1.25, 1.8, 18), 0.11),
            (0, ear(1.25, 1.8, -18), 0.11),
            (0, vbar(0, -1.15, 1.25), 0.16),
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


def best_line(backend, n):
    t = backend.target
    cz = {tuple(sorted(k)): p.error for k, p in t["cz"].items()
          if p is not None and p.error}
    ro = {k[0]: p.error for k, p in t["measure"].items()
          if p is not None and p.error}
    adj = {}
    for a, b in cz:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    best, best_cost = None, None

    def walk(path, cost):
        nonlocal best, best_cost
        if len(path) == n:
            c = cost + sum(ro.get(q, 0.05) for q in path)
            if best_cost is None or c < best_cost:
                best, best_cost = list(path), c
            return
        for nxt in adj.get(path[-1], []):
            if nxt not in path:
                edge = cz.get(tuple(sorted((path[-1], nxt))), 0.05)
                walk(path + [nxt], cost + edge)

    for start in adj:
        walk([start], 0.0)
    return best


def tomography(qc, backend, basis, shots, qubits=None):
    from qiskit import transpile
    from qiskit_experiments.library import StateTomography

    if basis:
        qc = transpile(qc, basis_gates=basis, optimization_level=3)
        exp = StateTomography(qc, physical_qubits=qubits)
        data = exp.run(backend, shots=shots)
    else:
        from qiskit_ibm_runtime import SamplerV2
        sampler = SamplerV2(backend)
        sampler.options.default_shots = shots
        sampler.options.dynamical_decoupling.enable = True
        exp = StateTomography(qc, physical_qubits=qubits)
        data = exp.run(backend=backend, sampler=sampler)
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


def median_error(target, name):
    import statistics
    errs = [p.error for p in target[name].values() if p is not None and p.error]
    return statistics.median(errs)


def get_backend():
    if BACKEND == "aer":
        from qiskit_aer import AerSimulator
        return AerSimulator(), ["rz", "sx", "x", "cx"]
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService()
    names = [b.name for b in service.backends(operational=True)]
    if "ibm_kingston" in names:
        real = service.backend("ibm_kingston")
    else:
        real = service.least_busy(operational=True)
    print(f"backend: {real.name}")
    if BACKEND == "ibm":
        return real, None
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    e1 = median_error(real.target, "sx")
    e2 = median_error(real.target, "cz")
    ro = median_error(real.target, "measure")
    print(f"noise: sx={e1:.2e} cz={e2:.2e} readout={ro:.2e}")
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(e1, 1), ["sx", "x"])
    nm.add_all_qubit_quantum_error(depolarizing_error(e2, 2), ["cz"])
    nm.add_all_qubit_readout_error(ReadoutError([[1 - ro, ro], [ro, 1 - ro]]))
    sim = AerSimulator(noise_model=nm, method="density_matrix")
    return sim, ["rz", "sx", "x", "cz"]


def main():
    backend, basis = get_backend()
    panels = {}
    for name in GLYPHS:
        nq, lim, comps = scene(name)
        qubits = None
        if basis is None:
            qubits = best_line(backend, nq)
            print(f"physical qubits: {qubits}")
        H = np.zeros((RES, RES))
        for i, (fock, U, w) in enumerate(comps):
            print(f"{name}: component {i + 1}/{len(comps)}")
            rho = tomography(component_circuit(nq, fock, U), backend, basis,
                             SHOTS, qubits)
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
