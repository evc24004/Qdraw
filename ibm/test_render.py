import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import render


def state_of(qc):
    from qiskit.quantum_info import Statevector

    return Statevector(qc).data


def test_state_preparation_matches_target():
    for name in ("HUSKY4", "U", "O", "HUSKY"):
        nq, _, components = render.scene(name)
        for fock, unitary, _ in components:
            target = unitary[:, fock]
            for build in (render.component_circuit,
                          render.component_state_circuit):
                vec = state_of(build(nq, fock, unitary))
                fid = abs(np.vdot(target, vec)) ** 2
                assert fid > 1 - 1e-9, (name, build.__name__, fid)
    print("ok: both circuit forms prepare the target states (fid > 1-1e-9)")


def test_gate_counts():
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    nq, _, components = render.scene("HUSKY4")
    line = CouplingMap.from_line(nq)
    print("component  unitary(logical/line)  state(logical/line)")
    for index, (fock, unitary, _) in enumerate(components, start=1):
        row = []
        for build in (render.component_circuit,
                      render.component_state_circuit):
            qc = build(nq, fock, unitary)
            logical = transpile(qc, basis_gates=["rz", "sx", "x", "cz"],
                                optimization_level=3)
            routed = transpile(qc, basis_gates=["rz", "sx", "x", "cz"],
                               coupling_map=line, optimization_level=3)
            two = lambda c: sum(1 for i in c.data
                                if i.operation.num_qubits == 2)
            row.append((two(logical), two(routed)))
        (ul, ur), (sl, sr) = row
        print(f"{index:>9}  {ul:>7}/{ur:<12} {sl:>6}/{sr}")
        assert sl < ul / 5, "state preparation should be far cheaper"


def test_husimi_analytic():
    d = 16
    vacuum = np.zeros((d, d))
    vacuum[0, 0] = 1
    img = render.husimi(vacuum, 3.6, 101)
    peak = np.unravel_index(np.argmax(img), img.shape)
    assert peak == (50, 50), peak
    assert abs(img[50, 50] - 1.0) < 1e-6

    ring = np.zeros((d, d))
    ring[3, 3] = 1
    img = render.husimi(ring, 3.6, 201)
    axis = np.linspace(-3.6, 3.6, 201)
    q, p = np.meshgrid(axis, axis)
    radius = np.sqrt(q**2 + p**2)
    peak_radius = radius.ravel()[np.argmax(img.ravel())]
    assert abs(peak_radius - np.sqrt(2 * 3)) < 0.15, peak_radius
    print("ok: Husimi matches vacuum peak and |3> ring radius sqrt(6)")


def test_sampler_options_schema():
    from qiskit_ibm_runtime.options import SamplerOptions

    options = SamplerOptions()
    options.dynamical_decoupling.enable = True
    options.twirling.enable_gates = True
    options.twirling.enable_measure = True
    options.twirling.num_randomizations = 8
    options.twirling.shots_per_randomization = 64
    print("ok: runtime sampler accepts the DD and twirling options")


def test_rehearsal_ab():
    backend = render.rehearsal_backend(render.DEFAULT_CALIBRATION,
                                       [148, 149, 150, 151])
    basis = ["rz", "sx", "x", "cz"]
    nq, _, components = render.scene("HUSKY4")
    print("component  unitary fid/purity   state fid/purity")
    for index, (fock, unitary, _) in enumerate(components, start=1):
        target = unitary[:, fock]
        row = []
        for build in (render.component_circuit,
                      render.component_state_circuit):
            rho = render.tomography(build(nq, fock, unitary), backend,
                                    basis, 512)
            fid = float(np.real(target.conj() @ rho @ target))
            pur = float(np.real(np.trace(rho @ rho)))
            row.append((fid, pur))
        (uf, up), (sf, sp) = row
        print(f"{index:>9}  {uf:.3f}/{up:.3f}        {sf:.3f}/{sp:.3f}")
        assert sf > uf, "state preparation should win under the noise model"


if __name__ == "__main__":
    test_state_preparation_matches_target()
    test_husimi_analytic()
    test_sampler_options_schema()
    test_gate_counts()
    test_rehearsal_ab()
    print("all tests passed")
