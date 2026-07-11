import argparse
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import render

INSTRUCTION_LIMIT = 24_000_000


def coupling_adjacency(backend):
    adjacency = {}
    for a, b in backend.target.build_coupling_map().get_edges():
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def enumerate_lines(adjacency, n):
    lines = set()

    def walk(path):
        if len(path) == n:
            if path[0] < path[-1]:
                lines.add(tuple(path))
            return
        for nxt in adjacency.get(path[-1], ()):
            if nxt not in path:
                walk(path + [nxt])

    for start in adjacency:
        walk([start])
    return sorted(lines)


def usage_profile(nq, components):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    line = CouplingMap.from_line(nq)
    edge_cz = np.zeros(nq - 1)
    qubit_1q = np.zeros(nq)
    for fock, unitary, _ in components:
        qc = transpile(render.component_state_circuit(nq, fock, unitary),
                       basis_gates=["rz", "sx", "x", "cz"],
                       coupling_map=line, optimization_level=3)
        for inst in qc.data:
            qubits = [qc.find_bit(q).index for q in inst.qubits]
            if inst.operation.name == "cz":
                edge_cz[min(qubits)] += 1
            elif inst.operation.name in ("sx", "x"):
                qubit_1q[qubits[0]] += 1
    count = len(components)
    return edge_cz / count, qubit_1q / count


def target_errors(backend):
    t = backend.target
    cz = {tuple(sorted(k)): p.error for k, p in t["cz"].items()
          if p is not None and p.error}
    sx = {k[0]: p.error for k, p in t["sx"].items()
          if p is not None and p.error}
    ro = {k[0]: p.error for k, p in t["measure"].items()
          if p is not None and p.error}
    return cz, sx, ro


def line_score(line, edge_cz, qubit_1q, cz, sx, ro):
    log_survival = 0.0
    for position, (a, b) in enumerate(zip(line, line[1:])):
        error = cz.get(tuple(sorted((a, b))), 0.05)
        log_survival += edge_cz[position] * math.log1p(-min(error, 0.5))
    for position, qubit in enumerate(line):
        log_survival += qubit_1q[position] * math.log1p(
            -min(sx.get(qubit, 0.01), 0.5))
        log_survival += math.log1p(-min(ro.get(qubit, 0.05), 0.5))
    return math.exp(log_survival)


def compatible(line, chosen, adjacency):
    occupied = set()
    for other in chosen:
        occupied.update(other)
        for qubit in other:
            occupied.update(adjacency.get(qubit, ()))
    return not (set(line) & occupied)


def pick_lines(backend, nq, k, components):
    adjacency = coupling_adjacency(backend)
    edge_cz, qubit_1q = usage_profile(nq, components)
    cz, sx, ro = target_errors(backend)
    candidates = [
        (max(line_score(line, edge_cz, qubit_1q, cz, sx, ro),
             line_score(line[::-1], edge_cz, qubit_1q, cz, sx, ro)),
         line if line_score(line, edge_cz, qubit_1q, cz, sx, ro)
         >= line_score(line[::-1], edge_cz, qubit_1q, cz, sx, ro)
         else line[::-1])
        for line in enumerate_lines(adjacency, nq)
    ]
    candidates.sort(reverse=True)

    best = None
    for start in range(min(40, len(candidates))):
        chosen = [candidates[start][1]]
        scores = [candidates[start][0]]
        for score, line in candidates:
            if len(chosen) == k:
                break
            if compatible(line, chosen, adjacency):
                chosen.append(line)
                scores.append(score)
        if len(chosen) < k:
            continue
        key = (min(scores), sum(scores))
        if best is None or key > best[0]:
            best = (key, chosen, scores)
    if best is None:
        raise RuntimeError(f"no {k} buffered disjoint lines found")
    return best[1], best[2]


def instruction_estimate(circuits, shots):
    per_qubit = {}
    for qc in circuits:
        for inst in qc.data:
            for qubit in inst.qubits:
                index = qc.find_bit(qubit).index
                per_qubit[index] = per_qubit.get(index, 0) + 1
    if not per_qubit:
        return 0
    return 2 * max(per_qubit.values()) * shots


def build_parallel(scene_name, component_indices, lines, backend):
    from qiskit_experiments.framework import ParallelExperiment
    from qiskit_experiments.library import StateTomography

    nq, _, components = render.scene(scene_name)
    children = []
    for component_index, line in zip(component_indices, lines):
        fock, unitary, _ = components[component_index]
        qc = render.component_state_circuit(nq, fock, unitary)
        children.append(StateTomography(qc, physical_qubits=line,
                                        backend=backend))
    parallel = ParallelExperiment(children, flatten_results=False,
                                  backend=backend)
    parallel.set_experiment_options(max_circuits=64)
    return parallel


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Parallel 5-qubit tomography on disjoint lines."
    )
    parser.add_argument("--scene", default="HUSKY")
    parser.add_argument("--mode", choices=("dry", "screen", "full"),
                        default="dry")
    parser.add_argument("--concurrency", type=int, default=8,
                        help="lines used in screen mode (1, 2, 4 or 8)")
    parser.add_argument("--shots", type=int)
    parser.add_argument("--backend-name", default="ibm_kingston")
    parser.add_argument("--output", type=pathlib.Path,
                        default=pathlib.Path(__file__).parent / "parallel.png")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    shots = args.shots or (128 if args.mode == "screen" else 512)

    from qiskit_ibm_runtime import QiskitRuntimeService

    backend = QiskitRuntimeService().backend(args.backend_name)
    nq, limit, components = render.scene(args.scene)
    k = args.concurrency if args.mode == "screen" else len(components)
    if args.mode == "screen":
        component_indices = [4] * k
    else:
        component_indices = list(range(len(components)))

    lines, scores = pick_lines(backend, nq, k, components)
    print(f"backend: {backend.name}")
    for line, score in zip(lines, scores):
        print(f"  line {list(line)}: predicted survival {score:.3f}")

    parallel = build_parallel(args.scene, component_indices, lines, backend)
    circuits = parallel._transpiled_circuits()
    per_job = math.ceil(len(circuits) / math.ceil(len(circuits) / 64))
    estimate = instruction_estimate(circuits[:per_job], shots)
    print(f"{len(circuits)} circuits, <= {per_job} per job, "
          f"busiest-qubit instruction estimate per job {estimate:,}")
    if estimate > INSTRUCTION_LIMIT:
        raise SystemExit("instruction estimate exceeds the guard; "
                         "reduce shots or split further")
    if args.mode == "dry":
        depths = [qc.depth() for qc in circuits]
        two_qubit = [sum(1 for inst in qc.data
                         if inst.operation.num_qubits == 2)
                     for qc in circuits]
        print(f"depth {min(depths)}-{max(depths)}, "
              f"2q gates {min(two_qubit)}-{max(two_qubit)} per circuit "
              "(not submitted)")
        return

    from qiskit_ibm_runtime import SamplerV2

    sampler = SamplerV2(backend)
    sampler.options.default_shots = shots
    sampler.options.dynamical_decoupling.enable = True
    data = parallel.run(backend=backend, sampler=sampler)
    data.block_for_results()

    image = np.zeros((render.RES, render.RES))
    for child, component_index, line in zip(
            data.child_data(), component_indices, lines):
        results = child.analysis_results("state", dataframe=True)
        rho = np.asarray(results["value"].iloc[0].data)
        fock, unitary, weight = components[component_index]
        target = unitary[:, fock]
        fidelity = float(np.real(target.conj() @ rho @ target))
        print(f"  line {list(line)}: fidelity vs target {fidelity:.3f}")
        if args.mode == "full":
            image += weight * render.husimi(rho, limit)

    if args.mode == "full":
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots(figsize=(4, 4))
        axis.imshow(image, origin="lower", cmap="inferno")
        axis.set_title(args.scene)
        axis.axis("off")
        figure.tight_layout()
        figure.savefig(args.output, dpi=150, facecolor="black")
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
