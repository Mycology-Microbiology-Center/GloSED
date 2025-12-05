#!/usr/bin/env python3
"""
Convert data to BIOM 2.1 HDF5 format with metadata

This script takes:
1. Long-format OTU table (columns: OTU, SampleID, Abundance); TSV format, can be gzip-compressed
2. OTU taxonomic annotations (columns: OTU, 7 taxonomic levels); TSV format, can be gzip-compressed
3. Sample metadata (SampleID + metadata columns); TSV format, can be gzip-compressed
4. OTU sequences (optional, FASTA format, can be gzip-compressed)

And creates a BIOM 2.1 HDF5 file (should be compatible with QIIME2)

Usage:
    python biom_from_longtab.py \
        --otu-table       otu_table.tsv.gz \
        --taxonomy        taxonomy.tsv.gz \
        --sample-metadata sample_metadata.tsv \
        --otu-fasta       otu_sequences.fasta.gz \
        --output          biom_table.biom
"""

import argparse
import pandas as pd
import numpy as np
from scipy import sparse
import biom
import gzip
import os
import sys
import re


## ANSI color codes
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    RESET = '\033[0m'


def _print_success(message):
    """Print success message in green."""
    print(f"{Colors.GREEN}{message}{Colors.RESET}")


def _print_error(message):
    """Print error message in red."""
    print(f"{Colors.RED}{message}{Colors.RESET}", file=sys.stderr)


def _is_gzipped(filepath):
    """Check if file is gzip compressed by reading magic bytes."""
    try:
        with open(filepath, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except:
        return False


def _find_column(columns, possible_names):
    """Find column name from list of possible names (case-insensitive)."""
    columns_lower = [col.lower() for col in columns]
    for name in possible_names:
        if name.lower() in columns_lower:
            return columns[columns_lower.index(name.lower())]
    return None


def _normalize_metadata_key(name):
    """Normalize metadata keys to be QIIME2-compatible.

    - Strip leading '#' and whitespace
    - Replace whitespace with underscores
    - Replace non-alphanumeric/underscore chars with underscores
    - Lowercase
    """
    name = str(name)
    ## drop leading comment characters and surrounding whitespace
    name = re.sub(r'^#+', '', name).strip()
    ## collapse whitespace to single underscore
    name = re.sub(r'\s+', '_', name)
    ## replace anything non-alphanumeric/underscore with underscore
    name = re.sub(r'[^\w]', '_', name)
    return name.lower()


def load_otu_sequences(filepath):
    """Load OTU sequences from FASTA (optionally gzipped).

    Returns
    -------
    dict
        Mapping {OTU_id: sequence}
    """
    ## Detect if file is gzip compressed
    if filepath.endswith('.gz') or (os.path.exists(filepath) and _is_gzipped(filepath)):
        opener = lambda p: gzip.open(p, 'rt')
    else:
        opener = lambda p: open(p, 'r')

    sequences = {}
    current_id = None
    current_seq_chunks = []

    with opener(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                ## flush previous record
                if current_id is not None:
                    seq = ''.join(current_seq_chunks).upper()
                    if current_id in sequences and sequences[current_id] != seq:
                        raise ValueError(f"Duplicate OTU ID '{current_id}' in FASTA with conflicting sequences.")
                    sequences[current_id] = seq
                header = line[1:].strip()
                current_id = header.split()[0]
                current_seq_chunks = []
            else:
                current_seq_chunks.append(line)

    ## flush last record
    if current_id is not None:
        seq = ''.join(current_seq_chunks).upper()
        if current_id in sequences and sequences[current_id] != seq:
            raise ValueError(f"Duplicate OTU ID '{current_id}' in FASTA with conflicting sequences.")
        sequences[current_id] = seq

    return sequences


def natural_sort_key(s):
    """Generate a sort key for natural sorting (handles numbers in strings)."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]


def load_long_format_otu_table(filepath):
    """Load long-format OTU table (SampleID, OTU, Abundance)."""
    ## Detect if file is gzip compressed
    if filepath.endswith('.gz') or (os.path.exists(filepath) and _is_gzipped(filepath)):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t')

    if len(df.columns) < 3:
        raise ValueError("OTU table must have at least 3 columns")

    ## Detect columns by name (case-insensitive)
    sample_col = _find_column(df.columns, ['sample', 'sampleid', 'sample_id', 'SampleID', 'Sample'])
    otu_col = _find_column(df.columns, ['otu', 'otuid', 'otu_id', 'OTU', 'feature', 'sequence'])
    abundance_col = _find_column(df.columns, ['abundance', 'count', 'reads', 'Abundance', 'Count'])

    if not all([sample_col, otu_col, abundance_col]):
        raise ValueError("Could not identify SampleID, OTU, and Abundance columns in OTU table. "
                        f"Available columns: {list(df.columns)}")

    return df, sample_col, otu_col, abundance_col


def load_taxonomy(filepath):
    """Load OTU taxonomy annotations."""
    ## Detect if file is gzip compressed
    if filepath.endswith('.gz') or (os.path.exists(filepath) and _is_gzipped(filepath)):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t')

    ## Find OTU column by name
    otu_col = _find_column(df.columns, ['otu', 'otuid', 'otu_id', 'OTU', 'feature', 'sequence'])
    if not otu_col:
        raise ValueError("Could not identify OTU column in taxonomy file. "
                        f"Available columns: {list(df.columns)}")

    return df.set_index(otu_col)


def load_sample_metadata(filepath):
    """Load sample metadata."""
    ## Detect if file is gzip compressed
    if filepath.endswith('.gz') or (os.path.exists(filepath) and _is_gzipped(filepath)):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t')

    ## Find SampleID column by name
    sample_col = _find_column(df.columns, ['sample', 'sampleid', 'sample_id', 'SampleID', 'Sample'])
    if not sample_col:
        raise ValueError("Could not identify SampleID column in sample metadata file. "
                        f"Available columns: {list(df.columns)}")

    return df.set_index(sample_col)


def create_sparse_matrix(long_df, sample_col, otu_col, abundance_col):
    """Convert long format to sparse matrix."""
    ## Get unique samples and OTUs
    samples = sorted(long_df[sample_col].unique(), key=natural_sort_key)
    otus = sorted(long_df[otu_col].unique())

    ## Create sample and OTU index mappings
    sample_to_idx = {sample: i for i, sample in enumerate(samples)}
    otu_to_idx = {otu: i for i, otu in enumerate(otus)}

    ## Create sparse matrix
    n_otus, n_samples = len(otus), len(samples)
    matrix = sparse.lil_matrix((n_otus, n_samples), dtype=float)

    ## Fill matrix
    for _, row in long_df.iterrows():
        otu_idx = otu_to_idx[row[otu_col]]
        sample_idx = sample_to_idx[row[sample_col]]
        matrix[otu_idx, sample_idx] = row[abundance_col]

    return matrix.tocsr(), otus, samples


def create_observation_metadata(taxonomy_df, otu_sequences=None):
    """Create observation metadata dictionary."""
    ## Treat all non-index columns as ordered taxonomy ranks
    tax_columns = list(taxonomy_df.columns)

    obs_metadata = {}
    for otu_id, row in taxonomy_df.iterrows():
        ## Create taxonomy list (only include non-empty levels)
        ## by BIOM/QIIME: list of ranks in order, e.g. ["k__Fungi", "p__Ascomycota", ...]
        taxonomy = []
        for col in tax_columns:
            if pd.notna(row[col]):
                value = str(row[col]).strip()
                if value:
                    taxonomy.append(value)

        ## Only include a single canonical 'taxonomy' field, which is what
        ## QIIME2-style tooling expects from BIOM observation metadata.
        metadata = {}
        if taxonomy:
            metadata['taxonomy'] = taxonomy

        ## Optionally attach underlying OTU sequence as additional observation
        ## metadata. This is not required by QIIME2, but can be useful as
        ## feature metadata after import.
        if otu_sequences is not None:
            seq = otu_sequences.get(otu_id)
            if seq:
                metadata['sequence'] = seq

        obs_metadata[otu_id] = metadata

    return obs_metadata


def create_sample_metadata(sample_df):
    """Create sample metadata dictionary."""
    ## Get all column names to ensure consistent metadata across samples
    all_columns = list(sample_df.columns)

    ## Normalize metadata keys to be QIIME2-style (lowercase, no spaces/#, etc.)
    normalized_columns = {
        col: _normalize_metadata_key(col) for col in all_columns
    }

    sample_metadata = {}
    for sample_id, row in sample_df.iterrows():
        metadata = {}
        for col in all_columns:
            key = normalized_columns[col]
            if pd.notna(row[col]):
                metadata[key] = row[col]
            else:
                ## Include None for missing values to maintain consistent metadata structure
                metadata[key] = None
        sample_metadata[sample_id] = metadata

    return sample_metadata


def validate_inputs(otu_df, taxonomy_df, sample_col, otu_col, abundance_col):
    """Validate OTU table and taxonomy for BIOM/QIIME2 compatibility.

    Checks:
    - abundances are numeric, integer, and non-negative
    - no duplicated OTU-sample combinations
    - all OTUs present in the OTU table have taxonomy
    """
    df = otu_df.copy()

    ## Ensure abundance is numeric
    try:
        df[abundance_col] = pd.to_numeric(df[abundance_col], errors='raise')
    except Exception as e:
        raise ValueError(f"Abundance column '{abundance_col}' must be numeric: {e}")

    ## Ensure abundance values are integers
    # Convert to float, then check integer equality to avoid dtype issues
    values = df[abundance_col].astype(float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Abundance column contains non-finite values (NaN/inf), which are invalid.")

    if not np.all(values == values.astype(int)):
        raise ValueError("Abundance values must be integers (no fractional counts allowed).")

    df[abundance_col] = values.astype(int)

    ## Ensure no negative abundances
    if (df[abundance_col] < 0).any():
        bad = df[df[abundance_col] < 0].iloc[0]
        raise ValueError(
            "Negative abundances detected. Example offending row: "
            f"Sample={bad[sample_col]}, OTU={bad[otu_col]}, {abundance_col}={bad[abundance_col]}"
        )

    ## Ensure no duplicated OTU–sample combinations
    dup_mask = df.duplicated(subset=[sample_col, otu_col], keep=False)
    if dup_mask.any():
        dup_rows = df.loc[dup_mask, [sample_col, otu_col, abundance_col]]
        examples = dup_rows.head().to_dict(orient='records')
        raise ValueError(
            "Found duplicated OTU–Sample combinations in the OTU table. "
            "Each (SampleID, OTU) pair must appear at most once.\n"
            f"Example duplicates (showing up to 5): {examples}"
        )

    ## Ensure all OTUs in the OTU table have taxonomy
    otu_ids = set(df[otu_col].unique())
    tax_ids = set(taxonomy_df.index)
    missing_in_tax = sorted(otu_ids - tax_ids, key=str)
    if missing_in_tax:
        example = missing_in_tax[:5]
        raise ValueError(
            "Some OTUs present in the OTU table are missing from the taxonomy file. "
            "This would result in incomplete BIOM observation metadata and is not "
            "allowed for QIIME2-compatible output.\n"
            f"Missing OTUs (showing up to 5 of {len(missing_in_tax)}): {example}"
        )

    return df


def validate_otu_sequences(otu_df, otu_col, otu_sequences):
    """Validate that OTU sequences cover all OTUs present in the OTU table."""
    otu_ids = set(otu_df[otu_col].unique())
    seq_ids = set(otu_sequences.keys())
    missing_seqs = sorted(otu_ids - seq_ids, key=str)
    if missing_seqs:
        example = missing_seqs[:5]
        raise ValueError(
            "Some OTUs present in the OTU table are missing from the OTU FASTA file. "
            "When providing sequences, every OTU in the BIOM table must have a sequence "
            "for QIIME2-compatible output.\n"
            f"OTUs without sequences (showing up to 5 of {len(missing_seqs)}): {example}"
        )


def main():
    parser = argparse.ArgumentParser(description="Convert long-format OTU table to BIOM 2.1 HDF5")
    parser.add_argument('--otu-table', required=True, help='Long-format OTU table (SampleID, OTU, Abundance)')
    parser.add_argument('--taxonomy', required=True, help='OTU taxonomy annotations')
    parser.add_argument('--sample-metadata', required=True, help='Sample metadata')
    parser.add_argument('--output', required=True, help='Output BIOM file path')
    parser.add_argument('--table-type', default='OTU table', help='Table type (default: OTU table)')
    parser.add_argument('--otu-fasta', required=False,
                        help='Optional OTU sequences in FASTA format (optionally gzipped) to store as observation metadata')

    args = parser.parse_args()

    try:
        ## Load data
        print("Loading OTU table...")
        otu_df, sample_col, otu_col, abundance_col = load_long_format_otu_table(args.otu_table)
        print("Loading taxonomy...")
        taxonomy_df = load_taxonomy(args.taxonomy)
        print("Loading sample metadata...")
        sample_df = load_sample_metadata(args.sample_metadata)

        otu_sequences = None
        if args.otu_fasta is not None:
            print("Loading OTU sequences from FASTA...")
            otu_sequences = load_otu_sequences(args.otu_fasta)

        ## Validate and normalize inputs before proceeding
        print("Validating input tables...")
        otu_df = validate_inputs(otu_df, taxonomy_df, sample_col, otu_col, abundance_col)
        if otu_sequences is not None:
            print("Validating OTU sequences against OTU table...")
            validate_otu_sequences(otu_df, otu_col, otu_sequences)
        _print_success(f"Loaded OTU table with {len(otu_df)} rows, columns: {sample_col}, {otu_col}, {abundance_col}")
        _print_success(f"  Unique samples: {otu_df[sample_col].nunique()}")
        _print_success(f"  Unique OTUs: {otu_df[otu_col].nunique()}")
        _print_success(f"  Total reads: {int(otu_df[abundance_col].sum())}")
        _print_success(f"Loaded taxonomy for {len(taxonomy_df)} OTUs")
        _print_success(f"Loaded metadata for {len(sample_df)} samples")

        ## Create sparse matrix
        print("Creating sparse matrix...")
        matrix, otus, samples = create_sparse_matrix(otu_df, sample_col, otu_col, abundance_col)
        _print_success(f"Created sparse matrix: {len(otus)} OTUs x {len(samples)} samples")

        ## Create metadata
        print("Creating observation metadata...")
        obs_metadata = create_observation_metadata(taxonomy_df, otu_sequences=otu_sequences)
        _print_success("Created observation metadata")

        print("Creating sample metadata...")
        sample_metadata = create_sample_metadata(sample_df)
        _print_success("Created sample metadata")

        ## Create observation metadata list (must be in same order as otus)
        obs_metadata_list = [obs_metadata.get(otu, {}) for otu in otus]

        ## Create sample metadata list (must be in same order as samples)
        sample_metadata_list = [sample_metadata.get(sample, {}) for sample in samples]

        ## Create BIOM table
        print("Creating BIOM table...")
        table = biom.Table(matrix, otus, samples,
                          observation_metadata=obs_metadata_list,
                          sample_metadata=sample_metadata_list,
                          type=args.table_type)
        _print_success("Created BIOM table")

        ## Save as HDF5 (BIOM 2.1 format)
        print(f"Saving to {args.output}...")
        biom.save_table(table, args.output, format_='2.1.0')
        _print_success(f"Saved BIOM file to {args.output}")

        _print_success("BIOM 2.1 HDF5 file created successfully!")
        print(f"  Table shape: {len(otus)} OTUs x {len(samples)} samples")
        print(f"  Non-zero entries: {matrix.nnz}")

    except Exception as e:
        _print_error(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
