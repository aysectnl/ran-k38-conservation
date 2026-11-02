#!/usr/bin/env python3
"""Report the residue found in the same alignment column as a reference residue.

Usage:
    python residue_column.py data/all.aln.fasta RAN_HUMAN 38

Finds residue 38 (1-based, ungapped numbering) of the reference sequence,
maps it to its alignment column, and prints the residue, its own
sequence-specific residue number, and a +/-4 column window for every sequence.
"""
import sys


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main(path, ref_key, position):
    seqs = read_fasta(path)
    ref = next(k for k in seqs if ref_key in k)
    count, col = 0, None
    for i, aa in enumerate(seqs[ref]):
        if aa != "-":
            count += 1
            if count == position:
                col = i
                break
    if col is None:
        sys.exit(f"{ref} has fewer than {position} residues")
    print(f"{ref} residue {position} ({seqs[ref][col]}) is at alignment column {col + 1}\n")
    for name, seq in seqs.items():
        aa = seq[col]
        own_number = len(seq[: col + 1].replace("-", "")) if aa != "-" else "-"
        window = seq[max(col - 4, 0): col + 5]
        print(f"{name:28s} {aa}  residue #{own_number:<4} window {window}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
