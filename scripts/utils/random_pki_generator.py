import argparse
import csv
import os
import random
import sys
from math import log10, floor

def round_sig(x: float, sig: int) -> float:
    if x == 0:
        return 0.0
    return round(x, sig - 1 - floor(log10(abs(x))))

def main():
    parser = argparse.ArgumentParser(
        description="Generate random pKi values and save to a CSV."
    )
    parser.add_argument(
        "--upper", type=float, required=True,
        help="Upper limit for pKi (exclusive of nothing; sampled uniformly in [lower, upper])."
    )
    parser.add_argument(
        "--num", type=int, required=True,
        help="Number of pKi values to generate."
    )
    parser.add_argument(
        "--out", required=True,
        help="Output CSV file path."
    )
    args = parser.parse_args()

    values = [random.uniform(1, args.upper) for _ in range(args.num)]
    rows = [(i, v) for i, v in enumerate(values, start=1)]

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "pKi"])
        writer.writerows(rows)

if __name__ == "__main__":
    main()