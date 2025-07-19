from __future__ import annotations

import re
from typing import NamedTuple


class AlleleComponents(NamedTuple):
    """Components of an allele name."""
    allele_id: str  # e.g., "01", "02"
    allele_name: str  # e.g., "deletion01", "deletion_01_allele"
    suffix: str  # e.g., "intact", "indels", "SV"
    percent: str  # e.g., "5.0", "100"


def extract_allele_from_header(header: str) -> str:
    """
    Extract the allele name from a header string.

    Header format: allele{ID}_{allele_name}_{suffix}_{percent}%
    Examples:
        - allele01_25003_Tombola_TMF2635-2636_intact_100% -> 25003_Tombola_TMF2635-2636
        - allele01_control_intact_100% -> control
        - allele02_deletion01_SV_75% -> deletion01
        - allele01_deletion_01_allele_intact_5.0% -> deletion_01_allele

    Args:
        header (str): The header string to parse

    Returns:
        str: The extracted allele name
    """
    # Pattern to match: allele{digits}_{allele_name}_{suffix}_{percent}%
    # suffix can be: intact, indels, SV
    # percent can be integer or decimal
    pattern = r"^allele\d+_(.+)_(intact|indels|SV)_\d+(?:\.\d+)?%$"
    match = re.match(pattern, header)

    if match:
        return match.group(1)

    # Fallback to original method for unexpected formats
    return header.split("_")[1] if "_" in header else header


def parse_allele_name(allele_name: str) -> AlleleComponents | None:
    """
    Parse an allele name into its components.

    Args:
        allele_name (str): The full allele name to parse
            e.g., "allele01_deletion_01_allele_intact_5.0%"

    Returns:
        AlleleComponents or None: Parsed components if successful, None otherwise
    """
    # Pattern to match the full allele name structure
    pattern = r"^allele(\d+)_(.+)_(intact|indels|SV)_(\d+(?:\.\d+)?)%$"
    match = re.match(pattern, allele_name)

    if match:
        return AlleleComponents(
            allele_id=match.group(1),
            allele_name=match.group(2),
            suffix=match.group(3),
            percent=match.group(4)
        )

    return None


def build_allele_name(components: AlleleComponents) -> str:
    """
    Build an allele name from its components.

    Args:
        components (AlleleComponents): The components to build from

    Returns:
        str: The complete allele name
    """
    return f"allele{components.allele_id}_{components.allele_name}_{components.suffix}_{components.percent}%"


def get_allele_base_name(allele_name: str) -> str:
    """
    Get the base name of an allele for grouping purposes.

    Args:
        allele_name (str): The full allele name

    Returns:
        str: The base name (allele ID + allele name)
    """
    components = parse_allele_name(allele_name)
    if components:
        return f"allele{components.allele_id}_{components.allele_name}"

    # Fallback for unexpected formats
    parts = allele_name.split('_')
    if len(parts) >= 2:
        return '_'.join(parts[:2])
    return allele_name
