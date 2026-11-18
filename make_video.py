#!/usr/bin/env python3
"""Make a video from simulation.bin (written by isithasubasinghe_bhut.cpp).

Usage: python3 make_video.py [simulation.bin] [-o simulation.mp4] [--fps 10] [--fraction 0.8]

The camera is fixed. Its view is a cube that contains FRACTION of all the
planet positions in all the frames. Planets outside the cube are not drawn.
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    """Return (mass[n], pos[frames][n][3]). The layout is in the .cpp file."""
    n, frames = np.fromfile(path, dtype=np.int32, count=2)
    mass = np.fromfile(path, dtype=np.float64, count=n, offset=8)
    pos = np.fromfile(path, dtype=np.float64, offset=8 + 8 * n)
    # A run that is not complete has fewer frames than the header says.
    frames = min(frames, pos.size // (3 * n))
    return mass, pos[: frames * n * 3].reshape(frames, n, 3)


def fixed_view(pos, fraction):
    """Return (centre[3], half_width) of the cube that contains `fraction` of the positions."""
    points = pos.reshape(-1, 3)
    centre = np.median(points, axis=0)
    # Largest offset on an axis: a point is in the cube when this is <= half_width.
    offset = np.abs(points - centre).max(axis=1)
    return centre, np.quantile(offset, fraction)


def main():
    parser = argparse.ArgumentParser(description="Make a video from simulation.bin.")
    parser.add_argument("input", nargs="?", default="simulation.bin")
    parser.add_argument("-o", "--output", default="simulation.mp4")
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--fraction", type=float, default=0.8,
                        help="fraction of the positions that the view must contain")
    args = parser.parse_args()

    _, pos = load(args.input)
    frames, n, _ = pos.shape
    centre, half = fixed_view(pos, args.fraction)
    print(f"{n} planets, {frames} frames, centre {centre}, half width {half:.4g}")

    fig = plt.figure(figsize=(8, 8), dpi=150, facecolor="black")
    ax = fig.add_subplot(projection="3d", facecolor="black")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    with tempfile.TemporaryDirectory() as tmp:
        for i in range(frames):
            p = pos[i]
            inside = (np.abs(p - centre) <= half).all(axis=1)
            p = p[inside]

            ax.clear()
            ax.set_axis_off()
            ax.set_xlim(centre[0] - half, centre[0] + half)
            ax.set_ylim(centre[1] - half, centre[1] + half)
            ax.set_zlim(centre[2] - half, centre[2] + half)
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(elev=25, azim=45)
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=0.5, c="white", alpha=0.6,
                       linewidths=0, depthshade=False)
            ax.text2D(0.02, 0.97, f"step {i}", transform=ax.transAxes, color="white")
            fig.savefig(Path(tmp) / f"frame_{i:05d}.png", facecolor="black")
            print(f"frame {i + 1}/{frames}", end="\r", flush=True)
        print()

        subprocess.run(
            ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
             "-i", str(Path(tmp) / "frame_%05d.png"),
             "-c:v", "libx264", "-pix_fmt", "yuv420p", args.output],
            check=True,
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
