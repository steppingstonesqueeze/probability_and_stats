#!/usr/bin/env python3
"""
stick_break_dirichlet.py
Sample from Dirichlet(alpha) via *finite stick-breaking*.

Usage examples:
  python stick_break_dirichlet.py --alphas 1,2,3 --n 5 --seed 42
  python stick_break_dirichlet.py --alphas 0.5,0.5,0.5,0.5 --n 10000 --summary --seed 1 --output samples.csv
"""
import argparse
import csv
import sys
from typing import List
import numpy as np

def parse_alphas(s: str) -> List[float]:
    try:
        vals = [float(x) for x in s.split(",")]
        if len(vals) < 2:
            raise ValueError("Need at least 2 alphas for a Dirichlet.")
        if any(a <= 0 for a in vals):
            raise ValueError("All alphas must be > 0.")
        return vals
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))

def stick_break_once(alpha: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One sample from Dirichlet(alpha) using finite stick-breaking.
    alpha: shape (K,), all > 0
    """
    K = alpha.shape[0]
    pies = np.zeros(K, dtype=float)
    remaining = 1.0
    # For i=1..K-1: V_i ~ Beta(alpha_i, sum_{j>i} alpha_j)
    # Then pi_i = remaining * V_i, update remaining *= (1 - V_i)
    for i in range(K - 1):
        a_i = alpha[i]
        b_i = alpha[i+1:].sum()
        V_i = rng.beta(a_i, b_i)
        pies[i] = remaining * V_i
        remaining *= (1.0 - V_i)
    pies[K - 1] = remaining  # last piece
    return pies

def main():
    p = argparse.ArgumentParser(description="Sample from Dirichlet(alpha) via finite stick-breaking.")
    p.add_argument("--alphas", type=parse_alphas, required=True,
                   help="Comma-separated positive alphas, e.g. 1,2,3")
    p.add_argument("--n", type=int, default=1, help="Number of samples (rows) to draw.")
    p.add_argument("--seed", type=int, default=None, help="Random seed.")
    p.add_argument("--output", type=str, default=None, help="Optional CSV output path.")
    p.add_argument("--summary", action="store_true",
                   help="Print summary stats (empirical mean vs. theoretical).")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    alpha = np.array(args.alphas, dtype=float)
    K = alpha.shape[0]

    samples = np.zeros((args.n, K), dtype=float)
    for i in range(args.n):
        samples[i, :] = stick_break_once(alpha, rng)

    # Sanity: rows sum to 1 (within numerical tolerance)
    row_sums = samples.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-10):
        print("Warning: some rows do not sum to 1 within tolerance.", file=sys.stderr)

    # Output
    if args.output:
        with open(args.output, "w", newline="") as f:
            writer = csv.writer(f)
            header = [f"pi_{j+1}" for j in range(K)]
            writer.writerow(header)
            writer.writerows(samples.tolist())
        print(f"Wrote {args.n} samples to {args.output}")
    else:
        # Print to stdout (up to a reasonable limit)
        max_show = min(args.n, 10)
        for i in range(max_show):
            print(", ".join(f"{x:.6f}" for x in samples[i]))
        if args.n > max_show:
            print(f"... ({args.n - max_show} more rows)")

    if args.summary:
        emp_mean = samples.mean(axis=0)
        th_mean = alpha / alpha.sum()
        print("\nSummary (empirical vs theoretical means):")
        for j in range(K):
            print(f"pi_{j+1}: emp={emp_mean[j]:.6f}  th={th_mean[j]:.6f}")
        print(f"Total alpha: {alpha.sum():.6f}")

if __name__ == "__main__":
    main()
