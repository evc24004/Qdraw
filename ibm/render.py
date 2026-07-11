import argparse
import json
import math
import pathlib
import statistics

import numpy as np
from scipy.linalg import expm


RES = 128
ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_CALIBRATION = ROOT / "run" / "calibration_ibm_kingston.json"


def mode_ops(d):
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    n = np.diag(np.arange(d))
    displacement = lambda alpha: expm(alpha * a.conj().T - np.conj(alpha) * a)
    rotation = lambda theta: expm(1j * theta * n)
    squeeze = lambda z: expm(
        (np.conj(z) * a @ a - z * a.conj().T @ a.conj().T) / 2
    )
    return displacement, rotation, squeeze


def scene(name):
    if name == "HUSKY4":
        nq = 4
        D, R, S = mode_ops(2**nq)
        alpha = lambda q, p: (q + 1j * p) / np.sqrt(2)
        vertical = lambda q, p, s: D(alpha(q, p)) @ S(np.log(s))
        ear = lambda q, p, deg: (
            D(alpha(q, p)) @ R(np.radians(deg)) @ S(np.log(1.4))
        )
        components = [
            (0, D(alpha(0, 0.1)), 0.16),
            (1, D(alpha(0, 0.1)), 0.22),
            (0, ear(-1.25, 1.8, 18), 0.11),
            (0, ear(1.25, 1.8, -18), 0.11),
            (0, vertical(0, -1.15, 1.25), 0.16),
        ]
        return nq, 3.6, components

    nq = 5
    D, R, S = mode_ops(2**nq)
    alpha = lambda q, p: (q + 1j * p) / np.sqrt(2)
    vertical = lambda q, p, s: D(alpha(q, p)) @ S(np.log(s))
    horizontal = lambda q, p, s: D(alpha(q, p)) @ S(-np.log(s))
    diagonal = lambda q, p, s, deg: (
        D(alpha(q, p)) @ R(np.radians(deg)) @ S(np.log(s))
    )
    identity = np.eye(2**nq)
    scenes = {
        "U": [
            (0, vertical(-1.4, 0.7, 2.1), 1),
            (0, vertical(1.4, 0.7, 2.1), 1),
            (0, horizontal(0, -1.9, 2.0), 1),
        ],
        "C": [
            (0, horizontal(0.3, 2.1, 1.8), 1),
            (0, vertical(-1.4, 0, 2.1), 1),
            (0, horizontal(0.3, -2.1, 1.8), 1),
        ],
        "O": [(2, identity, 0.45), (3, identity, 0.55)],
        "N": [
            (0, vertical(-1.5, 0, 2.2), 1),
            (0, vertical(1.5, 0, 2.2), 1),
            (0, diagonal(0, 0, 2.9, 44), 1),
        ],
        "HUSKY": [
            (0, D(alpha(0, 0.3)), 0.16),
            (1, D(alpha(0, 0.3)), 0.22),
            (0, diagonal(-1.5, 2.5, 1.6, 15), 0.11),
            (0, diagonal(1.5, 2.5, 1.6, -15), 0.11),
            (0, vertical(0, -1.2, 1.3), 0.16),
            (0, horizontal(0, -2.3, 1.3), 0.08),
            (0, D(alpha(-1.7, -0.6)), 0.08),
            (0, D(alpha(1.7, -0.6)), 0.08),
        ],
    }
    scenes["N2"] = scenes["N"]
    if name not in scenes:
        raise ValueError(f"unknown scene: {name}")
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


def _safe_error(value, fallback):
    if value is None or not 0 <= value < 1:
        return fallback
    return value


def select_low_error_line(backend, n, expected_cz=240, expected_sx=680):
    """Choose a connected line using a circuit-weighted error heuristic."""
    target = backend.target
    cz = {
        tuple(sorted(qubits)): properties.error
        for qubits, properties in target["cz"].items()
        if properties is not None
    }
    sx = {
        qubits[0]: properties.error
        for qubits, properties in target["sx"].items()
        if properties is not None
    }
    readout = {
        qubits[0]: properties.error
        for qubits, properties in target["measure"].items()
        if properties is not None
    }
    adjacency = {}
    for a, b in cz:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    candidates = []

    def visit(path):
        if len(path) == n:
            edge_uses = expected_cz / (n - 1)
            one_qubit_uses = expected_sx / n
            cost = 0.0
            for a, b in zip(path, path[1:]):
                error = _safe_error(cz.get(tuple(sorted((a, b)))), 0.05)
                cost -= edge_uses * math.log1p(-error)
            for qubit in path:
                cost -= one_qubit_uses * math.log1p(
                    -_safe_error(sx.get(qubit), 0.01)
                )
                cost -= math.log1p(-_safe_error(readout.get(qubit), 0.05))
            candidates.append((cost, tuple(path)))
            return
        for next_qubit in adjacency.get(path[-1], ()):
            if next_qubit not in path:
                visit(path + [next_qubit])

    for start in adjacency:
        visit([start])
    if not candidates:
        raise RuntimeError(f"backend has no connected line of {n} qubits")
    return list(min(candidates)[1])


def tomography(qc, backend, basis_gates, shots, qubits=None, enable_dd=True):
    from qiskit import transpile
    from qiskit_experiments.library import StateTomography

    if basis_gates:
        qc = transpile(qc, basis_gates=basis_gates, optimization_level=3)
        experiment = StateTomography(qc)
        data = experiment.run(backend, shots=shots)
    else:
        from qiskit_ibm_runtime import SamplerV2

        sampler = SamplerV2(backend)
        sampler.options.default_shots = shots
        sampler.options.dynamical_decoupling.enable = enable_dd
        experiment = StateTomography(qc, physical_qubits=qubits)
        data = experiment.run(backend=backend, sampler=sampler)
    data.block_for_results()
    results = data.analysis_results("state", dataframe=True)
    return np.asarray(results["value"].iloc[0].data)


def husimi(rho, limit, resolution=RES):
    d = rho.shape[0]
    axis = np.linspace(-limit, limit, resolution)
    q, p = np.meshgrid(axis, axis)
    alpha = (q.ravel() + 1j * p.ravel()) / np.sqrt(2)
    levels = np.arange(d)
    log_factorial = np.cumsum(np.log(np.maximum(levels, 1)))
    radius = np.abs(alpha)[:, None]
    with np.errstate(divide="ignore"):
        log_magnitude = np.where(
            radius > 0,
            levels * np.log(np.where(radius > 0, radius, 1)),
            0.0,
        )
    log_magnitude = np.where(
        (radius == 0) & (levels > 0), -np.inf, log_magnitude
    )
    coherent = np.exp(-radius**2 / 2 + log_magnitude - log_factorial / 2)
    coherent = coherent * np.exp(1j * levels * np.angle(alpha)[:, None])
    image = np.einsum(
        "pi,ij,pj->p", coherent.conj(), rho, coherent
    ).real
    return np.clip(image, 0, None).reshape(resolution, resolution)


def _calibration_errors(path, qubits):
    calibration = json.loads(path.read_text())
    selected = set(qubits)
    sx_errors = []
    cz_errors = []
    readout_errors = []
    for gate in calibration["gates"]:
        error = next(
            (
                value["value"]
                for value in gate["parameters"]
                if value["name"] == "gate_error"
            ),
            None,
        )
        gate_qubits = gate["qubits"]
        if error is None:
            continue
        if gate["gate"] == "sx" and gate_qubits[0] in selected:
            sx_errors.append(error)
        elif gate["gate"] == "cz" and set(gate_qubits).issubset(selected):
            cz_errors.append(error)
    for qubit in qubits:
        error = next(
            value["value"]
            for value in calibration["qubits"][qubit]
            if value["name"] == "readout_error"
        )
        readout_errors.append(error)
    return (
        statistics.median(sx_errors),
        statistics.median(cz_errors),
        statistics.median(readout_errors),
    )


def rehearsal_backend(calibration_path, qubits):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error

    sx_error, cz_error, readout_error = _calibration_errors(
        calibration_path, qubits
    )
    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(
        depolarizing_error(sx_error, 1), ["sx", "x"]
    )
    noise.add_all_qubit_quantum_error(depolarizing_error(cz_error, 2), ["cz"])
    noise.add_all_qubit_readout_error(
        ReadoutError(
            [
                [1 - readout_error, readout_error],
                [readout_error, 1 - readout_error],
            ]
        )
    )
    print(
        "archived noise medians: "
        f"sx={sx_error:.2e} cz={cz_error:.2e} readout={readout_error:.2e}"
    )
    return AerSimulator(noise_model=noise, method="density_matrix")


def parse_qubits(value):
    qubits = [int(item) for item in value.split(",") if item.strip()]
    if len(qubits) != len(set(qubits)):
        raise argparse.ArgumentTypeError("physical qubits must be unique")
    return qubits


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Render the quantum-CFS scenes with Qiskit tomography."
    )
    parser.add_argument(
        "--backend",
        choices=("aer", "rehearse", "ibm"),
        default="aer",
        help="aer is noiseless and offline; rehearse uses the archived calibration",
    )
    parser.add_argument(
        "--glyph",
        action="append",
        choices=("HUSKY4", "U", "C", "O", "N", "N2", "HUSKY"),
        dest="glyphs",
    )
    parser.add_argument("--shots", type=int, default=512)
    parser.add_argument("--resolution", type=int, default=RES)
    parser.add_argument(
        "--output", type=pathlib.Path, default=ROOT / "ibm" / "render.png"
    )
    parser.add_argument("--backend-name", default="ibm_kingston")
    parser.add_argument("--physical-qubits", type=parse_qubits)
    parser.add_argument("--calibration", type=pathlib.Path, default=DEFAULT_CALIBRATION)
    parser.add_argument("--no-dd", action="store_true")
    args = parser.parse_args(argv)
    if args.shots < 1 or args.resolution < 2:
        parser.error("shots must be positive and resolution must be at least 2")
    args.glyphs = args.glyphs or ["HUSKY4"]
    return args


def get_backend(args, n_qubits):
    if args.backend == "aer":
        from qiskit_aer import AerSimulator

        return AerSimulator(), ["rz", "sx", "x", "cx"], None
    if args.backend == "rehearse":
        qubits = args.physical_qubits or [148, 149, 150, 151]
        if len(qubits) != n_qubits:
            raise ValueError(f"{n_qubits} physical qubits are required")
        return (
            rehearsal_backend(args.calibration, qubits),
            ["rz", "sx", "x", "cz"],
            None,
        )

    from qiskit_ibm_runtime import QiskitRuntimeService

    backend = QiskitRuntimeService().backend(args.backend_name)
    qubits = args.physical_qubits or select_low_error_line(backend, n_qubits)
    if len(qubits) != n_qubits:
        raise ValueError(f"{n_qubits} physical qubits are required")
    print(f"backend: {backend.name}")
    print(f"physical qubits: {qubits}")
    return backend, None, qubits


def main(argv=None):
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    qubit_counts = {scene(name)[0] for name in args.glyphs}
    if args.backend == "ibm" and len(qubit_counts) != 1:
        raise ValueError("an IBM run must use scenes with the same qubit count")
    backend, basis_gates, physical_qubits = get_backend(args, max(qubit_counts))

    panels = {}
    for name in args.glyphs:
        nq, limit, components = scene(name)
        if nq != max(qubit_counts) and args.backend != "aer":
            raise ValueError("mixed qubit counts are supported only by the Aer backend")
        image = np.zeros((args.resolution, args.resolution))
        for index, (fock, unitary, weight) in enumerate(components, start=1):
            print(f"{name}: component {index}/{len(components)}")
            rho = tomography(
                component_circuit(nq, fock, unitary),
                backend,
                basis_gates,
                args.shots,
                physical_qubits,
                enable_dd=not args.no_dd,
            )
            image += weight * husimi(rho, limit, args.resolution)
        panels[name] = image
        np.save(args.output.with_name(f"{name}.npy"), image)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4))
    for axis, (name, image) in zip(np.atleast_1d(axes), panels.items()):
        axis.imshow(image, origin="lower", cmap="inferno")
        axis.set_title(name)
        axis.axis("off")
    figure.tight_layout()
    figure.savefig(args.output, dpi=150, facecolor="black")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
