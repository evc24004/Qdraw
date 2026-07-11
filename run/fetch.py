import json
import pathlib
import subprocess
import sys

import numpy as np
from qiskit import qasm3
from qiskit_ibm_runtime import QiskitRuntimeService

JOBS = {
    "d98ra8af47jc73a896ng": "single-component test (muzzle)",
    "d98rg3if47jc73a89ct0": "component 1: head |0>",
    "d98rgcgtcv6s73dmgbhg": "component 2: head |1>",
    "d98rglgtcv6s73dmgbsg": "component 3: left ear",
    "d98rgugtcv6s73dmgc60": "component 4: right ear",
    "d98rh74qp3as739tajjg": "component 5: muzzle",
}

here = pathlib.Path(__file__).parent
(here / "counts").mkdir(exist_ok=True)
(here / "circuits").mkdir(exist_ok=True)

service = QiskitRuntimeService()
manifest = {}
run_time = None
for job_id, label in JOBS.items():
    job = service.job(job_id)
    m = job.metrics()
    if run_time is None:
        run_time = job.creation_date
    manifest[job_id] = {
        "label": label,
        "backend": job.backend().name,
        "created": str(job.creation_date),
        "timestamps": m.get("timestamps"),
        "usage": m.get("usage"),
    }

    res = job.result()
    manifest[job_id]["num_circuits"] = len(res)
    recs = []
    for pub in res:
        recs.append({
            "m_idx": pub.metadata["circuit_metadata"]["m_idx"],
            "counts": pub.data.c_tomo.get_counts(),
        })
    with open(here / "counts" / f"{job_id}.json", "w") as f:
        json.dump(recs, f)
    manifest[job_id]["shots_per_circuit"] = sum(recs[0]["counts"].values())

    inputs = job.inputs
    manifest[job_id]["runtime_options"] = json.loads(
        json.dumps(inputs.get("options", {}), default=str))
    qc = inputs["pubs"][0][0]
    measured = sorted(
        (qc.find_bit(inst.clbits[0]).index, qc.find_bit(inst.qubits[0]).index)
        for inst in qc.data if inst.operation.name == "measure")
    manifest[job_id]["physical_qubits"] = [q for _, q in measured]
    ops = qc.count_ops()
    manifest[job_id]["sample_circuit"] = {
        "depth": qc.depth(),
        "ops": {k: int(v) for k, v in ops.items()},
    }
    with open(here / "circuits" / f"{job_id}_pub0.qasm", "w") as f:
        f.write(qasm3.dumps(qc))
    print(f"{job_id}  {label}: {len(res)} circuits saved")

with open(here / "jobs.json", "w") as f:
    json.dump(manifest, f, indent=2)

backend = service.backend("ibm_kingston")
try:
    props = backend.properties(datetime=run_time)
except Exception:
    props = backend.properties()
with open(here / "calibration_ibm_kingston.json", "w") as f:
    json.dump(props.to_dict(), f, default=str)

freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                        capture_output=True, text=True).stdout
(here / "requirements-lock.txt").write_text(freeze)
print("wrote jobs.json, calibration, requirements-lock.txt")
