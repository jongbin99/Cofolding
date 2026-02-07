import argparse
import os
import re
from typing import Iterable, Sequence, Tuple
import gemmi

EXCLUDE_RANGES: Tuple[Tuple[int, int], ...] = ((1, 3), (171, 173))
KEEP_TAGS: Tuple[str, ...] = ("LIG2", "LIG")


def convert_cif_to_pdb(base_dir: str) -> None:
    """Convert all .cif files under base_dir to .pdb (mirrors subdirs)."""
    for root, _, files in os.walk(base_dir):
        for fname in files:
            if not fname.lower().endswith(".cif"):
                continue

            m = re.search(r"mac-x(\d+)", fname)
            if not m:
                continue

            in_path = os.path.join(root, fname)
            out_name = f"mac-x{m.group(1)}_pred_chain.pdb"
            out_path = os.path.join(root, out_name)

            structure = gemmi.read_structure(in_path)
            structure.write_pdb(out_path)


def _should_exclude(residue_number: int, ranges: Sequence[Tuple[int, int]]) -> bool:
    return any(start <= residue_number <= end for start, end in ranges)


def remove_residues(
    file_path: str,
    keep_tags: Iterable[str] = KEEP_TAGS,
    exclude_ranges: Sequence[Tuple[int, int]] = EXCLUDE_RANGES,
) -> None:
    with open(file_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]

    filtered: list[str] = []
    for line in atom_lines:
        if any(tag in line for tag in keep_tags):
            filtered.append(line)
            continue

        residue_number = int(line[22:26].strip())
        if _should_exclude(residue_number, exclude_ranges):
            continue

        filtered.append(line)

    with open(file_path, "w", encoding="utf-8") as handle:
        handle.writelines(filtered)


def renumber_residues(file_path: str, new_start: int = 3) -> None:
    resi_offset = None
    with open(file_path, "r", encoding="utf-8") as infile:
        lines = infile.readlines()

    with open(file_path, "w", encoding="utf-8") as outfile:
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            residue_number = int(line[22:26].strip())
            if resi_offset is None:
                resi_offset = new_start - residue_number

            new_number = residue_number + resi_offset
            outfile.write(f"{line[:22]}{new_number:4}{line[26:]}")


def process_directory(base_dir: str, recursive: bool, new_start: int) -> None:
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
        description="Convert Boltz CIFs to PDBs (optional), then filter + renumber PDBs in-place."
    )
    parser.add_argument(
        "--base-dir",
        required=True,
        help="Directory containing CIF/PDB files (conversion + processing happens here).",
    )
    parser.add_argument(
        "--convert-cif",
        action="store_true",
        help="Convert all .cif files under base-dir to .pdb before processing.",
    )
    parser.add_argument(
        "--new-start",
        type=int,
        default=1,
        help="Residue number to assign to the first remaining residue (default: 1).",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process files in subdirectories as well.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.convert_cif:
        convert_cif_to_pdb(args.base_dir)

    process_directory(
        args.base_dir,
        recursive=args.recursive,
        new_start=args.new_start,
    )


if __name__ == "__main__":
    main()

