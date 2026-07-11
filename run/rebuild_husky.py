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

RENDER_JOBS = {
    "d98rg3if47jc73a89ct0": 0,
    "d98rgcgtcv6s73dmgbhg": 1,
    "d98rglgtcv6s73dmgbsg": 2,
    "d98rgugtcv6s73dmgc60": 3,
    "d98rh74qp3as739tajjg": 4,
}

nq, lim, comps = render.scene("HUSKY4")
H = np.zeros((render.RES, render.RES))
for job_id, idx in RENDER_JOBS.items():
    rho = refit.refit_job(job_id)
    fock, U, w = comps[idx]
    H += w * render.husimi(rho, lim)
    print(f"{job_id}: component {idx + 1} weight {w}")

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
ax.imshow(H, origin="lower", cmap="inferno")
ax.set_title("HUSKY4")
ax.axis("off")
fig.tight_layout()
out = here / "husky_rebuilt.png"
fig.savefig(out, dpi=150, facecolor="black")

sha = hashlib.sha256(out.read_bytes()).hexdigest()
print(f"wrote {out.name}")
print(f"sha256 {sha}")
orig = here.parent / "output" / "husky_kingston.png"
if orig.exists():
    print(f"published image sha256 {hashlib.sha256(orig.read_bytes()).hexdigest()}")
