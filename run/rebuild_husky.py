import hashlib
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here.parent / "ibm"))
import render
import refit

RUNS = {
    "husky_rebuilt.png": {
        "d98rg3if47jc73a89ct0": 0,
        "d98rgcgtcv6s73dmgbhg": 1,
        "d98rglgtcv6s73dmgbsg": 2,
        "d98rgugtcv6s73dmgc60": 3,
        "d98rh74qp3as739tajjg": 4,
    },
    "husky_rebuilt_v2.png": {
        "d98sngotcv6s73dmhkjg": 0,
        "d98sno52su3c739jmb8g": 1,
        "d98snvkqp3as739tbrr0": 2,
        "d98so6t2su3c739jmbog": 3,
        "d98sohcqp3as739tbsgg": 4,
    },
}

nq, lim, comps = render.scene("HUSKY4")
for filename, render_jobs in RUNS.items():
    H = np.zeros((render.RES, render.RES))
    for job_id, idx in render_jobs.items():
        rho = refit.refit_job(job_id)
        fock, U, w = comps[idx]
        H += w * render.husimi(rho, lim)
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.imshow(H, origin="lower", cmap="inferno")
    ax.set_title("HUSKY4")
    ax.axis("off")
    fig.tight_layout()
    out = here / filename
    fig.savefig(out, dpi=150, facecolor="black")
    print(f"wrote {out.name}")
    print(f"sha256 {hashlib.sha256(out.read_bytes()).hexdigest()}")

for name in ("husky_kingston.png", "husky_kingston_v2.png"):
    published = here.parent / "output" / name
    if published.exists():
        digest = hashlib.sha256(published.read_bytes()).hexdigest()
        print(f"published {name} sha256 {digest}")
