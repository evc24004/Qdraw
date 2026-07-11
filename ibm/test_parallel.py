import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import render


def reconstruct(data):
    results = data.analysis_results("state", dataframe=True)
    return np.asarray(results["value"].iloc[0].data)


def run_parallel(scene_name, component_indices, lines, shots=2048):
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from qiskit_experiments.framework import ParallelExperiment
    from qiskit_experiments.library import StateTomography

    backend = AerSimulator()
    nq, _, components = render.scene(scene_name)
    children = []
    targets = []
    for component_index, line in zip(component_indices, lines):
        fock, unitary, _ = components[component_index]
        qc = transpile(render.component_state_circuit(nq, fock, unitary),
                       basis_gates=["rz", "sx", "x", "cx"],
                       optimization_level=3)
        children.append(StateTomography(qc, physical_qubits=line))
        targets.append(unitary[:, fock])
    parallel = ParallelExperiment(children, flatten_results=False)
    data = parallel.run(backend, shots=shots)
    data.block_for_results()

    singles = []
    for child in children:
        single = child.copy()
        single_data = single.run(backend, shots=shots)
        single_data.block_for_results()
        singles.append(reconstruct(single_data))

    for index, (child_data, target, single_rho) in enumerate(
            zip(data.child_data(), targets, singles)):
        rho = reconstruct(child_data)
        fid_target = float(np.real(target.conj() @ rho @ target))
        overlap = float(np.real(np.trace(rho @ single_rho)))
        norm = float(np.sqrt(np.real(np.trace(rho @ rho))
                             * np.real(np.trace(single_rho @ single_rho))))
        agreement = overlap / norm
        print(f"  line {lines[index]}: fidelity vs target {fid_target:.4f}, "
              f"agreement with sequential {agreement:.4f}")
        assert fid_target > 0.97, fid_target
        assert agreement > 0.98, agreement


if __name__ == "__main__":
    print("2 x 4-qubit lines (HUSKY4 head, muzzle):")
    run_parallel("HUSKY4", [0, 4], [(0, 1, 2, 3), (5, 6, 7, 8)])
    print("2 x 5-qubit lines (HUSKY ear, muzzle):")
    run_parallel("HUSKY", [2, 4], [(0, 1, 2, 3, 4), (6, 7, 8, 9, 10)])
    print("4 x 5-qubit lines (HUSKY head, ears, muzzle):")
    run_parallel("HUSKY", [0, 2, 3, 4],
                 [(0, 1, 2, 3, 4), (6, 7, 8, 9, 10),
                  (12, 13, 14, 15, 16), (18, 19, 20, 21, 22)])
    print("parallel marginalization validated")
