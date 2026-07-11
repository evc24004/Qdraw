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
    "d98sl4gtcv6s73dmhi60": "state-prep pilot: muzzle",
    "d98smot2su3c739jmab0": "state-prep pilot: muzzle, gate twirl",
    "d98sngotcv6s73dmhkjg": "state-prep 1: head |0>",
    "d98sno52su3c739jmb8g": "state-prep 2: head |1>",
    "d98snvkqp3as739tbrr0": "state-prep 3: left ear",
    "d98so6t2su3c739jmbog": "state-prep 4: right ear",
    "d98sohcqp3as739tbsgg": "state-prep 5: muzzle",
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


PARALLEL_JOBS = {
    "d98tf34qp3as739tcl10": "v3 screen 1-line part 1",
    "d98tf3kqp3as739tcl20": "v3 screen 1-line part 2",
    "d98tf3qf47jc73a8bel0": "v3 screen 1-line part 3",
    "d98tf4cqp3as739tcl4g": "v3 screen 1-line part 4",
    "d98tgrqf47jc73a8bhp0": "v3 screen 8-line part 1",
    "d98tgs0tcv6s73dmihog": "v3 screen 8-line part 2",
    "d98tgsif47jc73a8bhq0": "v3 screen 8-line part 3",
    "d98tgsotcv6s73dmihqg": "v3 screen 8-line part 4",
    "d98ti0cqp3as739tcqeg": "v3 parallel render part 1",
    "d98ti0if47jc73a8bjlg": "v3 parallel render part 2",
    "d98ti10tcv6s73dmik00": "v3 parallel render part 3",
    "d98ti1d2su3c739jn900": "v3 parallel render part 4",
}

parallel_manifest = {}
for job_id, label in PARALLEL_JOBS.items():
    job = service.job(job_id)
    m = job.metrics()
    entry = {
        "label": label,
        "backend": job.backend().name,
        "created": str(job.creation_date),
        "timestamps": m.get("timestamps"),
        "usage": m.get("usage"),
        "runtime_options": json.loads(
            json.dumps(job.inputs.get("options", {}), default=str)),
    }
    res = job.result()
    entry["num_circuits"] = len(res)
    pubs = []
    lines = None
    for pub in res:
        md = pub.metadata["circuit_metadata"]
        lines = md["composite_qubits"]
        registers = sorted(pub.data.keys(), key=lambda n: int(n[1:]))
        children = []
        for child_index, register in enumerate(registers):
            children.append({
                "m_idx": md["composite_metadata"][child_index]["m_idx"],
                "counts": getattr(pub.data, register).get_counts(),
            })
        pubs.append(children)
    entry["lines"] = lines
    with open(here / "counts" / f"{job_id}.json", "w") as f:
        json.dump(pubs, f)
    parallel_manifest[job_id] = entry
    print(f"{job_id}  {label}: {len(res)} circuits saved")

with open(here / "jobs_parallel.json", "w") as f:
    json.dump(parallel_manifest, f, indent=2)
print("wrote jobs_parallel.json")
