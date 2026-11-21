#!/usr/bin/env python3
"""CUDA regression checks; uses only Python's standard library and a real GPU."""
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest


BINARY = Path(sys.argv.pop(1) if len(sys.argv) > 1 else "build/bhut-gpu").resolve()


def direct_frames(rows, steps, dt, softening=0.0):
    """Independent FP64 all-pairs reference for the original kick/drift scheme."""
    position = [list(row[1:4]) for row in rows]
    velocity = [list(row[4:7]) for row in rows]
    frames = [[x for p in position for x in p]]
    for step in range(steps):
        acceleration = [[0.0] * 3 for _ in rows]
        for i in range(len(rows)):
            for j in range(len(rows)):
                if i == j:
                    continue
                delta = [position[j][k] - position[i][k] for k in range(3)]
                r2 = sum(x*x for x in delta) + softening*softening
                if r2:
                    scale = rows[j][0] / (r2 * math.sqrt(r2))
                    for k in range(3):
                        acceleration[i][k] += delta[k] * scale
        for i in range(len(rows)):
            for k in range(3):
                velocity[i][k] += acceleration[i][k] * dt * (0.5 if step == 0 else 1)
                position[i][k] += velocity[i][k] * dt
        frames.append([x for p in position for x in p])
    return frames


class GPUChecks(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="bhut-gpu-test-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "planets.txt"
        self.output = self.root / "simulation.bin"
        # Intentionally not in spatial order; masses identify original particle IDs.
        self.rows = [
            (0.75, 2.0, -1.0, 0.5, 0.01, 0.02, -0.01),
            (1.25, -2.0, 1.0, -0.5, -0.02, 0.01, 0.02),
            (0.5, 0.5, 2.0, -1.0, 0.03, -0.01, 0.01),
            (2.0, -0.5, -2.0, 1.0, -0.01, -0.02, -0.03),
        ]
        self.write_input(self.rows)

    def write_input(self, rows):
        self.source.write_text("".join(" ".join(map(str, row)) + "\n" for row in rows))

    def run_solver(self, *args, success=True):
        result = subprocess.run(
            [str(BINARY), "--input", str(self.source), "--output", str(self.output), *map(str, args)],
            text=True, capture_output=True, timeout=120,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def read_output(self):
        data = self.output.read_bytes()
        n, count = struct.unpack_from("=ii", data)
        self.assertEqual(len(data), 8 + 8*n + count*n*3*8)
        masses = struct.unpack_from(f"={n}d", data, 8)
        frames = [struct.unpack_from(f"={3*n}d", data, 8+8*n+frame*24*n) for frame in range(count)]
        self.assertTrue(all(math.isfinite(x) for frame in frames for x in frame))
        return masses, frames

    def compare_frames(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for frame, reference in zip(actual, expected):
            for value, target in zip(frame, reference):
                self.assertAlmostEqual(value, target, delta=3e-6 * max(1.0, abs(target)))

    def test_gpu_tree_and_force_cases(self):
        result = self.run_solver("--self-test")
        self.assertIn("Self-test 8/8 passed", result.stdout)
        self.assertFalse(self.output.exists())

    def test_file_layout_ids_and_integration(self):
        self.run_solver("--steps", 5, "--dt", 0.001, "--theta", 0, "--validate", 4,
                        "--max-force-error", 0.00001)
        masses, frames = self.read_output()
        self.assertEqual(masses, tuple(row[0] for row in self.rows))
        self.compare_frames(frames, direct_frames(self.rows, 5, 0.001))

    def test_softening_and_double_accumulation(self):
        self.run_solver("--steps", 5, "--dt", 0.001, "--theta", 0, "--softening", 0.2,
                        "--double-sum", "--validate", 4, "--max-force-error", 0.00001)
        _, frames = self.read_output()
        self.compare_frames(frames, direct_frames(self.rows, 5, 0.001, 0.2))

    def test_stride_includes_final_frame(self):
        self.run_solver("--steps", 5, "--dt", 0.001, "--theta", 0, "--output-stride", 2)
        _, frames = self.read_output()
        reference = direct_frames(self.rows, 5, 0.001)
        self.compare_frames(frames, [reference[i] for i in [0, 2, 4, 5]])

    def test_zero_steps_preserves_input_frame(self):
        self.rows = [(1.23456789123, 0.123456789123, -2.0, 0.0, 0.1, 0.2, 0.3)]
        self.write_input(self.rows)
        self.run_solver("--steps", 0)
        masses, frames = self.read_output()
        self.assertEqual(masses, (self.rows[0][0],))
        self.assertEqual(frames, [self.rows[0][1:4]])

    def test_single_particle_drifts_without_self_force(self):
        self.rows = [(1.0, 0.25, -0.5, 0.75, 0.1, 0.2, 0.3)]
        self.write_input(self.rows)
        self.run_solver("--steps", 5, "--dt", 0.01, "--validate", 1)
        _, frames = self.read_output()
        self.compare_frames(frames, direct_frames(self.rows, 5, 0.01))

    def test_existing_output_is_untouched(self):
        self.output.write_bytes(b"existing simulation")
        self.run_solver("--steps", 0, success=False)
        self.assertEqual(self.output.read_bytes(), b"existing simulation")

    def test_no_output(self):
        self.run_solver("--steps", 2, "--no-output")
        self.assertFalse(self.output.exists())

    def test_invalid_options(self):
        for args in [("--steps", "-1"), ("--steps", "2junk"), ("--dt", "nan"),
                     ("--dt", "0"), ("--theta", "1e30"), ("--softening", "1e-30"),
                     ("--output-stride", "0"), ("--unknown", "1")]:
            with self.subTest(args=args):
                self.run_solver(*args, success=False)
                self.assertFalse(self.output.exists())

    def test_invalid_input(self):
        for content in ["", "0 0 0 0 0 0 0\n", "-1 0 0 0 0 0 0\n", "1 0 0\n",
                        "1 nan 0 0 0 0 0\n", "1 0 0 0 0 0 0 extra\n", "1e-99 0 0 0 0 0 0\n"]:
            with self.subTest(content=content):
                self.source.write_text(content)
                self.run_solver("--steps", 0, success=False)
                self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
