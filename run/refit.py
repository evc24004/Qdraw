import json
import pathlib
import sys

import numpy as np

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here.parent / "ibm"))
import render

PAULI = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]),
         "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.array([[1, 0], [0, -1]])}
BASIS = "ZXY"

COMPONENT_INDEX = {
    "single-component test (muzzle)": 4,
    "component 1: head |0>": 0,
    "component 2: head |1>": 1,
    "component 3: left ear": 2,
    "component 4: right ear": 3,
    "component 5: muzzle": 4,
}


def linear_inversion(records, n=4):
    est, cnt = {}, {}
    for rec in records:
        basis = [BASIS[i] for i in rec["m_idx"]]
        counts = rec["counts"]
        shots = sum(counts.values())
        for mask in range(1, 2 ** n):
            key = tuple(basis[q] if (mask >> q) & 1 else "I" for q in range(n))
            val = 0.0
            for bits, c in counts.items():
                par = sum(int(bits[len(bits) - 1 - q]) for q in range(n)
                          if (mask >> q) & 1)
                val += ((-1) ** (par % 2)) * c / shots
            est[key] = est.get(key, 0.0) + val
            cnt[key] = cnt.get(key, 0) + 1
    rho = np.eye(2 ** n, dtype=complex) / 2 ** n
    for key, v in est.items():
        P = np.array([[1.0]])
        for q in reversed(range(n)):
            P = np.kron(P, PAULI[key[q]])
        rho += (v / cnt[key]) * P / 2 ** n
    return rho


def load_jobs():
    return json.loads((here / "jobs.json").read_text())


def refit_job(job_id):
    records = json.loads((here / "counts" / f"{job_id}.json").read_text())
    return linear_inversion(records)


if __name__ == "__main__":
    jobs = load_jobs()
    nq, lim, comps = render.scene("HUSKY4")
    (here / "states").mkdir(exist_ok=True)
    print(f"{'job':22s} {'label':32s} {'fidelity':>8s} {'purity':>7s}")
    for job_id, meta in jobs.items():
        rho = refit_job(job_id)
        fock, U, w = comps[COMPONENT_INDEX[meta["label"]]]
        psi = U[:, fock]
        fid = float(np.real(psi.conj() @ rho @ psi))
        pur = float(np.real(np.trace(rho @ rho)))
        print(f"{job_id:22s} {meta['label']:32s} {fid:8.3f} {pur:7.3f}")
        with open(here / "states" / f"{job_id}.json", "w") as f:
            json.dump({"label": meta["label"], "fidelity_vs_ideal": fid,
                       "purity": pur, "rho_real": np.real(rho).tolist(),
                       "rho_imag": np.imag(rho).tolist()}, f)
