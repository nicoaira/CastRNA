#!/usr/bin/env python3
"""
Generate MSA outputs from a multi-sequence Stockholm alignment.

This script takes a Stockholm alignment (from cmalign with multiple sequences)
and generates:
1. Stockholm file with per-sequence secondary structure annotations (#=GR SS)
2. Aligned FASTA file
3. TSV summary of structures

The secondary structures are projected from the consensus structure (SS_cons)
onto each sequence, handling gaps and insertions appropriately.
"""

import argparse
import json
import sys
from pathlib import Path


# Bracket pairs for structure notation
OPENING_BRACKETS = '(<[{'
CLOSING_BRACKETS = ')>]}'
BRACKET_PAIRS = dict(zip(OPENING_BRACKETS, CLOSING_BRACKETS))
REVERSE_BRACKET_PAIRS = dict(zip(CLOSING_BRACKETS, OPENING_BRACKETS))


def normalize_ss_cons(ss_cons: str) -> str:
    """Normalize SS_cons notation to standard dot-bracket."""
    result = []
    for char in ss_cons:
        if char in OPENING_BRACKETS or char in CLOSING_BRACKETS:
            result.append(char)
        elif char == '.':
            result.append('.')
        elif char in '_-:,~':
            result.append('.')
        elif char.isalpha():
            result.append(char)
        else:
            result.append('.')
    return ''.join(result)


def find_pairs(structure: str) -> dict[int, int]:
    """Find base pair positions in bracket notation."""
    pairs = {}
    stacks = {bracket: [] for bracket in OPENING_BRACKETS}
    letter_stacks = {chr(c): [] for c in range(ord('A'), ord('Z') + 1)}

    for i, char in enumerate(structure):
        if char in OPENING_BRACKETS:
            stacks[char].append(i)
        elif char in CLOSING_BRACKETS:
            opening = REVERSE_BRACKET_PAIRS[char]
            if stacks[opening]:
                j = stacks[opening].pop()
                pairs[i] = j
                pairs[j] = i
        elif char.isupper():
            letter_stacks[char].append(i)
        elif char.islower():
            upper = char.upper()
            if letter_stacks[upper]:
                j = letter_stacks[upper].pop()
                pairs[i] = j
                pairs[j] = i

    return pairs


def standardize_brackets(structure: str) -> str:
    """
    Convert all bracket types to standard parentheses ().

    Converts: <>, [], {} -> ()
    Preserves: Letter-based pseudoknots (A...a, B...b, etc.), dots, and dashes
    """
    result = []
    for char in structure:
        if char in OPENING_BRACKETS:
            result.append('(')
        elif char in CLOSING_BRACKETS:
            result.append(')')
        else:
            # Keep letters (pseudoknots), dots, and dashes as-is
            result.append(char)
    return ''.join(result)


def project_structure_aligned(aligned_seq: str, ss_cons: str) -> str:
    """
    Project consensus structure onto aligned sequence, keeping alignment gaps.

    Returns structure string with same length as aligned_seq, where:
    - Positions with residues get projected structure
    - Gap positions get '-' (to distinguish from '.' which means unpaired)
    """
    ss_cons = normalize_ss_cons(ss_cons)
    aln_pairs = find_pairs(ss_cons)

    projected = []
    for aln_col, (seq_char, ss_char) in enumerate(zip(aligned_seq, ss_cons)):
        if seq_char in '.-':
            # Gap in sequence - use '-' in structure (no nucleotide = no structure)
            projected.append('-')
        elif seq_char.islower():
            # Insertion state - unpaired
            projected.append('.')
        else:
            # Match state - check if pair is intact
            if aln_col in aln_pairs:
                partner_col = aln_pairs[aln_col]
                partner_char = aligned_seq[partner_col] if partner_col < len(aligned_seq) else '-'
                if partner_char in '.-':
                    # Partner is deleted - break pair
                    projected.append('.')
                else:
                    projected.append(ss_char)
            else:
                projected.append(ss_char)

    return ''.join(projected)


def parse_stockholm_msa(sto_path: Path) -> dict:
    """
    Parse a multi-sequence Stockholm alignment.

    Returns dict with:
    - sequences: {seq_id: aligned_sequence}
    - ss_cons: consensus structure
    - rf: reference annotation (if present)
    - gc_annotations: other GC lines
    """
    sequences = {}
    ss_cons_parts = []
    rf_parts = []
    gc_annotations = {}
    gf_annotations = []

    # For multi-block files
    seq_parts = {}

    with open(sto_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            if line.startswith('# STOCKHOLM') or line.startswith('//') or not line:
                continue

            if line.startswith('#=GF'):
                gf_annotations.append(line)
                continue

            if line.startswith('#=GC SS_cons'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    ss_cons_parts.append(parts[2])
                continue

            if line.startswith('#=GC RF'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    rf_parts.append(parts[2])
                continue

            if line.startswith('#=GC'):
                # Other GC annotations
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    key = parts[1]
                    if key not in gc_annotations:
                        gc_annotations[key] = []
                    gc_annotations[key].append(parts[2])
                continue

            if line.startswith('#'):
                continue

            # Sequence line
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id = parts[0]
                seq_data = parts[1].replace(' ', '')
                if seq_id not in seq_parts:
                    seq_parts[seq_id] = []
                seq_parts[seq_id].append(seq_data)

    # Concatenate parts
    for seq_id, parts in seq_parts.items():
        sequences[seq_id] = ''.join(parts)

    ss_cons = ''.join(ss_cons_parts) if ss_cons_parts else None
    rf = ''.join(rf_parts) if rf_parts else None

    # Concatenate other GC annotations
    for key in gc_annotations:
        gc_annotations[key] = ''.join(gc_annotations[key])

    return {
        'sequences': sequences,
        'ss_cons': ss_cons,
        'rf': rf,
        'gc_annotations': gc_annotations,
        'gf_annotations': gf_annotations
    }


def write_stockholm_with_ss(
    output_path: Path,
    sequences: dict[str, str],
    structures: dict[str, str],
    ss_cons: str,
    rf: str | None,
    rfam_id: str,
    gf_annotations: list[str] | None = None
):
    """Write Stockholm format with per-sequence SS annotations."""
    # Find max ID length for formatting
    max_id_len = max(len(seq_id) for seq_id in sequences.keys())
    max_id_len = max(max_id_len, len('SS_cons'), len('RF'))

    with open(output_path, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF AC   {rfam_id}\n")

        # Write any additional GF annotations
        if gf_annotations:
            for line in gf_annotations:
                if not line.startswith('#=GF AC'):
                    f.write(f"{line}\n")

        f.write("\n")

        # Write sequences with their SS annotations
        for seq_id, aligned_seq in sequences.items():
            f.write(f"{seq_id:<{max_id_len}}  {aligned_seq}\n")
            if seq_id in structures:
                f.write(f"#=GR {seq_id:<{max_id_len - 5}} SS  {structures[seq_id]}\n")

        f.write("\n")

        # Write consensus annotations
        f.write(f"#=GC {'SS_cons':<{max_id_len}}  {ss_cons}\n")
        if rf:
            f.write(f"#=GC {'RF':<{max_id_len}}  {rf}\n")

        f.write("//\n")


def write_aligned_fasta(output_path: Path, sequences: dict[str, str]):
    """Write aligned sequences in FASTA format."""
    with open(output_path, 'w') as f:
        for seq_id, aligned_seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Write in 80-char lines
            for i in range(0, len(aligned_seq), 80):
                f.write(f"{aligned_seq[i:i+80]}\n")


def write_structures_tsv(
    output_path: Path,
    sequences: dict[str, str],
    structures: dict[str, str]
):
    """Write structure summary TSV."""
    with open(output_path, 'w') as f:
        f.write("seq_id\tsequence_length\taligned_length\tstructure\tnum_pairs\tnum_unpaired\n")

        for seq_id, aligned_seq in sequences.items():
            ungapped_len = len(aligned_seq.replace('-', '').replace('.', ''))
            aligned_len = len(aligned_seq)

            if seq_id in structures:
                structure = structures[seq_id]
                # Count pairs (each pair counted once)
                pairs = find_pairs(structure)
                num_pairs = len(pairs) // 2
                num_unpaired = structure.count('.')
            else:
                structure = ''
                num_pairs = 0
                num_unpaired = 0

            # Get ungapped structure for output
            ungapped_structure = ''.join(
                s for s, c in zip(structures.get(seq_id, ''), aligned_seq)
                if c not in '.-'
            )

            f.write(f"{seq_id}\t{ungapped_len}\t{aligned_len}\t{ungapped_structure}\t{num_pairs}\t{num_unpaired}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MSA outputs from multi-sequence Stockholm alignment'
    )
    parser.add_argument('--input', type=Path, required=True, help='Input Stockholm file')
    parser.add_argument('--output-sto', type=Path, required=True, help='Output Stockholm file with SS')
    parser.add_argument('--output-fasta', type=Path, required=True, help='Output aligned FASTA file')
    parser.add_argument('--output-tsv', type=Path, required=True, help='Output structures TSV')
    parser.add_argument('--rfam-id', type=str, required=True, help='Rfam family ID')

    args = parser.parse_args()

    # Parse input Stockholm
    try:
        msa_data = parse_stockholm_msa(args.input)
    except Exception as e:
        sys.stderr.write(f"Error parsing Stockholm file: {e}\n")
        sys.exit(1)

    if not msa_data['sequences']:
        sys.stderr.write("Error: No sequences found in Stockholm file\n")
        sys.exit(1)

    if not msa_data['ss_cons']:
        sys.stderr.write("Error: No SS_cons line found in Stockholm file\n")
        sys.exit(1)

    # Project structure onto each sequence and standardize brackets
    structures = {}
    for seq_id, aligned_seq in msa_data['sequences'].items():
        projected = project_structure_aligned(aligned_seq, msa_data['ss_cons'])
        structures[seq_id] = standardize_brackets(projected)

    # Standardize SS_cons brackets as well
    ss_cons_standardized = standardize_brackets(normalize_ss_cons(msa_data['ss_cons']))

    # Write outputs
    write_stockholm_with_ss(
        args.output_sto,
        msa_data['sequences'],
        structures,
        ss_cons_standardized,
        msa_data['rf'],
        args.rfam_id,
        msa_data['gf_annotations']
    )

    write_aligned_fasta(args.output_fasta, msa_data['sequences'])
    write_structures_tsv(args.output_tsv, msa_data['sequences'], structures)

    sys.stderr.write(
        f"Generated MSA outputs for {args.rfam_id}: "
        f"{len(msa_data['sequences'])} sequences\n"
    )


if __name__ == '__main__':
    main()
