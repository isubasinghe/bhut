#!/usr/bin/env python3
"""Make a video from simulation.bin (written by isithasubasinghe_bhut.cpp).

Usage: python3 make_video.py [simulation.bin] [-o simulation.mp4] [--fps 10]
                           [--width 3840] [--height 2160]

Render smoothed logarithmic density at the output resolution with a fixed camera;
optionally draw individual particles. Colours track groups from the initial frame;
these are approximate galaxy memberships, since the file has no galaxy labels.
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


COLOURS = np.array([[0.25, 0.70, 1.0], [1.0, 0.55, 0.20]])


def project(pos, elev, azim):
    """Orthographic screen coordinates, using the same orientation as mplot3d."""
    elev, azim = np.deg2rad([elev, azim])
    right = [-np.sin(azim), np.cos(azim), 0]
    up = [-np.sin(elev) * np.cos(azim),
          -np.sin(elev) * np.sin(azim), np.cos(elev)]
    return pos @ np.array([right, up]).T


def initial_groups(pos):
    """Deterministic two-cluster fit; labels remain attached to particle IDs."""
    centre = np.mean(pos, axis=0)
    _, axes = np.linalg.eigh(np.cov(pos.T))
    along = (pos - centre) @ axes[:, -1]
    groups = along > np.median(along)
    for _ in range(50):
        if groups.all() or not groups.any():
            break
        centres = np.array([pos[~groups].mean(axis=0), pos[groups].mean(axis=0)])
        distances = ((pos[:, None, :] - centres) ** 2).sum(axis=2)
        updated = distances[:, 1] < distances[:, 0]
        if np.array_equal(groups, updated):
            break
        groups = updated
    return groups


def fixed_view(pos, fraction):
    """Square screen bounds containing at least `fraction` of all positions."""
    points = pos.reshape(-1, 2)
    centre = np.median(points, axis=0)
    offset = np.abs(points - centre).max(axis=1)
    return centre, max(float(np.quantile(offset, fraction)) * 1.05, 1e-12)


def density_layers(points, groups, centre, half, size):
    """Smoothed counts on a (width, height) grid with per-axis half extents."""
    bounds = [(c - h, c + h) for c, h in zip(centre, half)]
    layers = np.array([
        np.histogram2d(*points[groups == group].T, bins=size, range=bounds)[0].T
        for group in (False, True)
    ])
    offsets = np.arange(-4, 5)
    kernel = np.exp(-0.5 * (offsets / 1.4) ** 2)
    kernel /= kernel.sum()
    for axis in (1, 2):
        padding = [(0, 0)] * 3
        padding[axis] = (4, 4)
        padded = np.pad(layers, padding)
        layers = np.zeros_like(layers)
        for j, weight in enumerate(kernel):
            window = [slice(None)] * 3
            window[axis] = slice(j, j + layers.shape[axis])
            layers += weight * padded[tuple(window)]
    return layers


def colour_density(layers, exposure):
    total = layers.sum(axis=0)
    tint = np.einsum("ghw,gc->hwc", layers, COLOURS)
    tint /= np.maximum(total[..., None], 1e-12)
    # Boost faint outskirts while retaining structure in the dense cores.
    brightness = np.clip(np.log1p(10 * total) / np.log1p(10 * exposure), 0, 1) ** 0.6
    return tint * brightness[..., None]


def main():
    parser = argparse.ArgumentParser(description="Make a video from simulation.bin.")
    parser.add_argument("input", nargs="?", default="simulation.bin")
    parser.add_argument("-o", "--output", default="simulation.mp4")
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--width", type=int, default=3840,
                        help="output width in pixels, a positive even number (default: 3840)")
    parser.add_argument("--height", type=int, default=2160,
                        help="output height in pixels, a positive even number (default: 2160)")
    parser.add_argument("--render-mode", choices=("particles", "density"),
                        default="density", help="rendering style (default: density)")
    parser.add_argument("--particle-size", type=float, default=2.0,
                        help="particle marker diameter in output pixels (default: 2)")
    parser.add_argument("--fraction", type=float, default=0.95,
                        help="fraction of projected positions to include (default: 0.95)")
    parser.add_argument("--elev", type=float, default=25, help="camera elevation in degrees")
    parser.add_argument("--azim", type=float, default=45, help="camera azimuth in degrees")
    args = parser.parse_args()
    if args.fps <= 0 or not 0 < args.fraction <= 1:
        parser.error("fps must be positive and fraction must be in (0, 1]")
    if args.width <= 0 or args.height <= 0 or args.width % 2 or args.height % 2:
        parser.error("width and height must be positive even numbers for H.264 output")
    if not np.isfinite(args.particle_size) or args.particle_size <= 0:
        parser.error("particle size must be finite and positive")

    _, pos = load(args.input)
    frames, n, _ = pos.shape
    if frames == 0 or n < 2 or not np.isfinite(pos).all():
        parser.error("input must contain at least one finite frame with two particles")
    groups = initial_groups(pos[0])
    screen = project(pos, args.elev, args.azim)
    del pos
    centre, half = fixed_view(screen, args.fraction)
    print(f"{n} planets, {frames} frames, centre {centre}, half width {half:.4g}")
    print(f"Initial spatial groups: {(~groups).sum()} blue, {groups.sum()} orange")
    dpi = 150
    fig = plt.figure(figsize=(args.width / dpi, args.height / dpi),
                     dpi=dpi, facecolor="black")
    ax = fig.add_axes((0.02, 0.04, 0.96, 0.88), facecolor="black")
    ax.set_axis_off()
    if args.render_mode == "particles":
        # Scatter sizes are in points squared; keep marker diameter in pixels.
        picture = ax.scatter(*screen[0].T, s=(args.particle_size * 72 / dpi) ** 2,
                             c=COLOURS[groups.astype(int)], marker="o", linewidths=0)
        ax.set_xlim(centre[0] - half, centre[0] + half)
        ax.set_ylim(centre[1] - half, centre[1] + half)
        ax.set_aspect("equal")
        caption = "Individual particles · fixed camera"
    else:
        density_size = (args.width, args.height)
        # Expand the longer view axis so world distances stay equal in pixels.
        density_half = half * np.array(density_size) / min(density_size)
        first = density_layers(screen[0], groups, centre, density_half, density_size)
        exposure = max(float(first.sum(axis=0).max()), 1.0)
        # Fill the canvas: one density cell per output pixel, with text overlaid.
        ax.set_position((0, 0, 1, 1))
        picture = ax.imshow(colour_density(first, exposure), origin="lower",
                            interpolation="nearest")
        caption = "Projected particle density · logarithmic brightness · fixed camera"
    title = fig.text(0.04, 0.965, "", color="white", fontsize=13, va="top")
    fig.text(0.04, 0.925, "Initial group 1", color=COLOURS[0], fontsize=10)
    fig.text(0.24, 0.925, "Initial group 2", color=COLOURS[1], fontsize=10)
    fig.text(0.04, 0.025, caption, color="0.65", fontsize=9)

    with tempfile.TemporaryDirectory() as tmp:
        for i in range(frames):
            if args.render_mode == "particles":
                picture.set_offsets(screen[i])
            else:
                layers = first if i == 0 else density_layers(
                    screen[i], groups, centre, density_half, density_size)
                picture.set_data(colour_density(layers, exposure))
            title.set_text(f"Galaxy collision  |  step {i}")
            fig.savefig(Path(tmp) / f"frame_{i:05d}.png", facecolor="black")
            print(f"frame {i + 1}/{frames}", end="\r", flush=True)
        print()

        subprocess.run(
            ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
             "-i", str(Path(tmp) / "frame_%05d.png"),
             "-c:v", "libx264", "-crf", "18", "-pix_fmt", "yuv420p",
             "-vf", "setsar=1", "-movflags", "+faststart", args.output],
            check=True,
        )
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
