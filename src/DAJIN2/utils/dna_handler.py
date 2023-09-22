from __future__ import annotations


def revcomp(sequence: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement[nt] for nt in sequence[::-1])
