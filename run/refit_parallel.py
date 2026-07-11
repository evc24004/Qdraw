import hashlib
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here.parent / "ibm"))
sys.path.insert(0, str(here))
import render
import refit

RENDER_PARTS = [
    "d98ti0cqp3as739tcqeg",
    "d98ti0if47jc73a8bjlg",
    "d98ti10tcv6s73dmik00",
    "d98ti1d2su3c739jn900",
]
COMPONENT_NAMES = ["head |0>", "head |1>", "left ear", "right ear",
                   "muzzle", "chin", "left cheek", "right cheek"]

nq, limit, components = render.scene("HUSKY")
manifest = json.loads((here / "jobs_parallel.json").read_text())
lines = manifest[RENDER_PARTS[0]]["lines"]

records_per_component = [[] for _ in components]
for job_id in RENDER_PARTS:
    pubs = json.loads((here / "counts" / f"{job_id}.json").read_text())
    for children in pubs:
        for component_index, child in enumerate(children):
            records_per_component[component_index].append(child)

(here / "states").mkdir(exist_ok=True)
image = np.zeros((render.RES, render.RES))
print(f"{'component':12s} {'line':22s} {'fidelity':>8s} {'purity':>7s}")
for component_index, records in enumerate(records_per_component):
    rho = refit.linear_inversion(records, n=nq)
    fock, unitary, weight = components[component_index]
    target = unitary[:, fock]
    fidelity = float(np.real(target.conj() @ rho @ target))
    purity = float(np.real(np.trace(rho @ rho)))
    line = lines[component_index]
    print(f"{COMPONENT_NAMES[component_index]:12s} {str(line):22s} "
          f"{fidelity:8.3f} {purity:7.3f}")
    with open(here / "states" / f"v3_component_{component_index}.json",
              "w") as f:
        json.dump({"component": COMPONENT_NAMES[component_index],
                   "line": line, "fidelity_vs_ideal": fidelity,
                   "purity": purity, "rho_real": np.real(rho).tolist(),
                   "rho_imag": np.imag(rho).tolist()}, f)
    image += weight * render.husimi(rho, limit)

figure, axis = plt.subplots(figsize=(4, 4))
axis.imshow(image, origin="lower", cmap="inferno")
axis.set_title("HUSKY")
axis.axis("off")
figure.tight_layout()
out = here / "husky_rebuilt_v3.png"
figure.savefig(out, dpi=150, facecolor="black")
print(f"wrote {out.name}")
print(f"sha256 {hashlib.sha256(out.read_bytes()).hexdigest()}")
