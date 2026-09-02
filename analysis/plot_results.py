"""Plot HelixTracking CSV products without depending on a C++ plotting library."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("helixtracking/outputs"))
    args = parser.parse_args()
    hits = read_csv(args.input / "hits.csv")
    curve = read_csv(args.input / "fitted_curve.csv")
    inlier = hits["reconstructed_inlier"].astype(bool)

    fig = plt.figure(figsize=(11, 4.8))
    xy = fig.add_subplot(1, 2, 1)
    xy.scatter(hits["x_m"][~inlier], hits["y_m"][~inlier], s=20, color="#b7b7b7", label="rejected hit")
    xy.scatter(hits["x_m"][inlier], hits["y_m"][inlier], s=20, color="#315f8c", label="accepted hit")
    xy.plot(curve["x_m"], curve["y_m"], color="#d14b40", lw=2, label="fitted helix")
    xy.set_aspect("equal")
    xy.set_xlabel("x [m]")
    xy.set_ylabel("y [m]")
    xy.set_title("Transverse circle reconstruction")
    xy.legend(frameon=False, fontsize=8)

    rz = fig.add_subplot(1, 2, 2, projection="3d")
    rz.scatter(hits["x_m"][~inlier], hits["y_m"][~inlier], hits["z_m"][~inlier], s=10, color="#b7b7b7")
    rz.scatter(hits["x_m"][inlier], hits["y_m"][inlier], hits["z_m"][inlier], s=15, color="#315f8c")
    rz.plot(curve["x_m"], curve["y_m"], curve["z_m"], color="#d14b40", lw=2)
    rz.set_xlabel("x [m]")
    rz.set_ylabel("y [m]")
    rz.set_zlabel("z [m]")
    rz.set_title("Three-dimensional track")
    fig.tight_layout()
    fig.savefig(args.input / "track_reconstruction.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()

