# benchmarks/bench_parallel_mace.py
import time
import numpy as np
import torch
from ase.build import molecule
from mace.calculators import mace_off
from nebwalk import idpp_interpolate
from nebwalk.optimize import _eval_all

MODEL = "EGRET_1T.model"
N_IMAGES = 7
N_REPEATS = 10
N_WARMUP = 3

def make_images():
    start = molecule('C2H6')
    end = molecule('C2H6')
    end.rotate(60, 'z', center=end.positions[0])
    images = idpp_interpolate(start, end, N_IMAGES)
    for img in images:
        img.calc = mace_off(model=MODEL, device='cpu')
    return images

def bench(images, n_workers):
    for _ in range(N_WARMUP):
        _eval_all(images, n_workers)
    times = []
    for _ in range(N_REPEATS):
        t0 = time.perf_counter()
        _eval_all(images, n_workers)
        times.append(time.perf_counter() - t0)
    return np.mean(times) * 1000, np.std(times) * 1000

print(f"CPU cores visible: {torch.get_num_threads()} PyTorch threads\n")

# Build images once — same calc instances reused across both conditions
images = make_images()

for label, n_threads in [("default threads", None), ("1 thread/worker", 1)]:
    if n_threads is not None:
        torch.set_num_threads(n_threads)
    print(f"=== {label} ===")
    seq_ms, seq_std = bench(images, n_workers=1)
    par_ms, par_std = bench(images, n_workers=N_IMAGES)
    print(f"  Sequential: {seq_ms:.1f} ± {seq_std:.1f} ms")
    print(f"  Parallel:   {par_ms:.1f} ± {par_std:.1f} ms")
    print(f"  Speedup:    {seq_ms/par_ms:.2f}×\n")