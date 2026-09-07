#!/usr/bin/env python3
"""Extract confidence matrices and pLDDT values from Protenix full-data JSON."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path


PAIRWISE_FIELDS = ("token_pair_pae", "token_pair_pde", "contact_probs")


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def load_json(path: Path) -> dict:
    with open_text(path, "rt") as stream:
        return json.load(stream)


def write_matrix(path: Path, matrix) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(matrix)


def write_vector(path: Path, values) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("atom", "pLDDT"))
        writer.writerows(enumerate(values, start=1))


def sample_number(path: Path) -> int:
    match = re.search(r"sample_(\d+)", path.name)
    if not match:
        raise ValueError(f"Cannot determine sample number from {path.name}")
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    inputs = sorted(
        args.input_dir.glob("*_full_data_sample_*.json*"),
        key=sample_number,
    )
    if not inputs:
        raise FileNotFoundError(f"No Protenix full-data JSON files in {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for input_path in inputs:
        model = sample_number(input_path)
        print(f"Reading {input_path}")
        data = load_json(input_path)

        for field in PAIRWISE_FIELDS:
            if field not in data:
                raise KeyError(f"{field} is missing from {input_path}")
            short_name = {
                "token_pair_pae": "pae",
                "token_pair_pde": "pde",
                "contact_probs": "contact_probs",
            }[field]
            output = args.output_dir / f"{args.prefix}_{short_name}_m{model}.csv.gz"
            write_matrix(output, data[field])

        plddt_field = next(
            (field for field in ("atom_plddt", "atom_plddts") if field in data),
            None,
        )
        if plddt_field:
            output = args.output_dir / f"{args.prefix}_atom_plddt_m{model}.tsv"
            write_vector(output, data[plddt_field])


if __name__ == "__main__":
    main()
