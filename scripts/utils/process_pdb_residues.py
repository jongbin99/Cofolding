"""Utilities for trimming and renumbering PDB files.

This script replicates the residue filtering logic from the notebook
`handling_cif_pdb_files.ipynb` and exposes it as a command-line tool.

Usage
-----

    python process_pdb_residues.py --base-dir /path/to/pdbs

By default each PDB file in the provided directory is rewritten in-place.
Non-`ATOM`/`HETATM` records are removed, residues 1-3 and 171-173 are
discarded (unless they contain ``LIG``/``LIG2``) and the remaining residues
are renumbered so that the first residue starts at index 3.

Optional flags allow recursion into sub-directories or changing the starting
residue number.
"""

import argparse
import os
from typing import Iterable, Sequence, Tuple

EXCLUDE_RANGES: Tuple[Tuple[int, int], ...] = ((1, 3), (171, 173))
KEEP_TAGS: Tuple[str, ...] = ("LIG2", "LIG")

def _should_exclude(residue_number: int,
                    ranges: Sequence[Tuple[int, int]]) -> bool:
    return any(start <= residue_number <= end for start, end in ranges)


def remove_residues(file_path: str,
                    keep_tags: Iterable[str] = KEEP_TAGS,
                    exclude_ranges: Sequence[Tuple[int, int]] = EXCLUDE_RANGES
                    ) -> None:
    """Keep only ATOM/HETATM records and drop selected residue ranges."""

    with open(file_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]

    filtered: list[str] = []
    for line in atom_lines:
        if any(tag in line for tag in keep_tags):
            filtered.append(line)
            continue

        try:
            residue_number = int(line[22:26].strip())
        except ValueError:
            filtered.append(line)
            continue

        if _should_exclude(residue_number, exclude_ranges):
            continue

        filtered.append(line)

    with open(file_path, "w", encoding="utf-8") as handle:
        handle.writelines(filtered)

#reset your new_start to 3 if you want to renumber the residues starting from 1
def renumber_residues(file_path: str, new_start: int = 3) -> None:
    """Renumber residues so the first residue starts at ``new_start``."""

    resi_offset = None
    with open(file_path, "r", encoding="utf-8") as infile:
        lines = infile.readlines()

    with open(file_path, "w", encoding="utf-8") as outfile:
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            try:
                residue_number = int(line[22:26].strip())
            except ValueError:
                outfile.write(line)
                continue

            if resi_offset is None:
                resi_offset = new_start - residue_number

            new_number = residue_number + resi_offset
            outfile.write(f"{line[:22]}{new_number:4}{line[26:]}")


def process_directory(base_dir: str, recursive: bool, new_start: int) -> None:
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Base directory does not exist: {base_dir}")

    paths: list[str] = []
    if recursive:
        for root, _, files in os.walk(base_dir):
            for name in files:
                if name.lower().endswith(".pdb"):
                    paths.append(os.path.join(root, name))
    else:
        for name in os.listdir(base_dir):
            if name.lower().endswith(".pdb"):
                paths.append(os.path.join(base_dir, name))

    paths.sort()

    for path in paths:
        remove_residues(path)
        renumber_residues(path, new_start=new_start)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter and renumber PDB files in a directory."
    )
    parser.add_argument(
        "--base-dir",
        required=True,
        help="Directory containing PDB files to process."
    )
    parser.add_argument(
        "--new-start",
        type=int,
        default=3,
        help="Residue number to assign to the first remaining residue (default: 3)."
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process PDB files in subdirectories as well."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    process_directory(args.base_dir, recursive=args.recursive, new_start=args.new_start)


if __name__ == "__main__":
    main()

