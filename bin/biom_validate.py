#!/usr/bin/env python3
"""
Validate BIOM file

Usage:
    python biom_validate.py --input biom_table.biom
"""

import argparse
import biom


def validate_biom_file(filepath):
    """Validate BIOM file and print information."""
    print(f"Loading BIOM file: {filepath}")

    table = biom.load_table(filepath)

    print("Successfully loaded BIOM table")
    print(f"Format version: {table.format_version}")
    print(f"Table type: {table.type}")
    print(f"Shape: {table.shape[0]} observations × {table.shape[1]} samples")
    print(f"Non-zero entries: {table.nnz}")

    # Table.format_version is a tuple like (2, 1) coming from biom.util.__format_version__
    if table.format_version == (2, 1):
        print("BIOM 2.1 (HDF5) format confirmed")
    else:
        print(f"Warning: Format version is {table.format_version}, expected (2, 1) for BIOM 2.1")

    obs_md = table.metadata(axis='observation')
    sample_md = table.metadata(axis='sample')

    if obs_md and any(obs_md):
        print(f"\nObservation metadata present for {sum(1 for md in obs_md if md)} observations")
        for i, md in enumerate(obs_md):
            if md:
                print(f"  Example observation metadata (OTU {table.ids(axis='observation')[i]}): {md}")
                break

    if sample_md and any(sample_md):
        print(f"\nSample metadata present for {sum(1 for md in sample_md if md)} samples")
        for i, md in enumerate(sample_md):
            if md:
                print(f"  Example sample metadata (Sample {table.ids(axis='sample')[i]}): {md}")
                break

    if table.format_version == (2, 1):
        print("\nBIOM format: 2.1")
    else:
        print(f"\nBIOM format: {table.format_version}")
    
    if table.type in ('OTU table', 'FeatureTable'):
        print("...QIIME2-compatible table type")
    else:
        print(f"Table type '{table.type}' - may need to be 'OTU table' or 'FeatureTable' for QIIME2")

    print(
        "\nExample import command:\n"
        f"  qiime tools import --type 'FeatureTable[Frequency]' "
        f"--input-format BIOMV210Format --input-path {filepath}"
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Validate BIOM file")
    parser.add_argument('--input', required=True, help='BIOM file to validate')
    args = parser.parse_args()
    validate_biom_file(args.input)
