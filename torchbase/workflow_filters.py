"""Workflow filtering for suspect data handling (Issue #22).

This module provides functions for filtering alleles, loci, and profiles
based on quality data. Supports three levels of filtering:
1. Allele-level: remove specific suspect alleles
2. Locus-level: remove entire loci
3. Profile-level: remove loci used in suspect profiles

Design principles:
- Positive semantics: include suspect data by default
- Users opt-in to filtering with flags
- Filtering happens before computation (MinHash/alignment)
- Results document what was filtered for reproducibility
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any


def load_quality_data(torch_path: Path) -> Optional[Dict]:
    """Load quality.json from torch directory if present.

    Args:
        torch_path: Path to torch directory

    Returns:
        Parsed quality data dict, or None/empty dict if file not found

    Raises:
        json.JSONDecodeError: If quality.json is malformed
    """
    quality_path = Path(torch_path) / "quality.json"

    if not quality_path.exists():
        return None

    with open(quality_path) as f:
        return json.load(f)


def filter_alleles(
    alleles: List[str],
    suspect_alleles: List[str],
    include_suspect: bool = True
) -> List[str]:
    """Filter allele list based on suspect list.

    Args:
        alleles: List of allele names
        suspect_alleles: List of suspect allele names to exclude
        include_suspect: If True, keep all alleles. If False, remove suspects.

    Returns:
        Filtered allele list, preserving original order
    """
    if include_suspect or not suspect_alleles:
        return alleles

    suspect_set = set(suspect_alleles)
    return [a for a in alleles if a not in suspect_set]


def filter_alleles_by_locus(
    alleles_by_locus: Dict[str, List[str]],
    suspect_alleles: List[str],
    include_suspect: bool = True
) -> Dict[str, List[str]]:
    """Filter alleles by locus, removing suspect alleles.

    Args:
        alleles_by_locus: Dict mapping locus name to list of alleles
        suspect_alleles: List of suspect allele names to exclude
        include_suspect: If False, remove suspects

    Returns:
        Dict with same structure, filtered alleles
    """
    if include_suspect or not suspect_alleles:
        return alleles_by_locus

    result = {}
    suspect_set = set(suspect_alleles)

    for locus, alleles in alleles_by_locus.items():
        filtered = [a for a in alleles if a not in suspect_set]
        result[locus] = filtered

    return result


def filter_loci(
    alleles_by_locus: Dict[str, List[str]],
    suspect_loci: List[str],
    include_suspect: bool = True
) -> Dict[str, List[str]]:
    """Filter entire loci from allele database.

    Args:
        alleles_by_locus: Dict mapping locus name to list of alleles
        suspect_loci: List of suspect locus names to remove entirely
        include_suspect: If False, remove suspects

    Returns:
        Dict with suspect loci removed entirely
    """
    if include_suspect or not suspect_loci:
        return alleles_by_locus

    suspect_set = set(suspect_loci)
    return {
        locus: alleles
        for locus, alleles in alleles_by_locus.items()
        if locus not in suspect_set
    }


def filter_profiles(
    alleles_by_locus: Dict[str, List[str]],
    suspect_profiles: List[str],
    loci_in_profiles: Dict[str, List[str]],
    include_suspect: bool = True
) -> Dict[str, List[str]]:
    """Filter loci that appear in suspect profiles.

    Args:
        alleles_by_locus: Dict mapping locus name to list of alleles
        suspect_profiles: List of suspect profile identifiers
        loci_in_profiles: Dict mapping profile id to list of loci it uses
        include_suspect: If False, remove loci from suspect profiles

    Returns:
        Dict with loci from suspect profiles removed
    """
    if include_suspect or not suspect_profiles:
        return alleles_by_locus

    # Collect all loci used by suspect profiles
    suspect_profile_set = set(suspect_profiles)
    loci_to_remove = set()

    for profile_id, loci in loci_in_profiles.items():
        if profile_id in suspect_profile_set:
            loci_to_remove.update(loci)

    return {
        locus: alleles
        for locus, alleles in alleles_by_locus.items()
        if locus not in loci_to_remove
    }


def apply_filters(
    alleles_by_locus: Dict[str, List[str]],
    suspect_alleles: Optional[List[str]] = None,
    suspect_loci: Optional[List[str]] = None,
    suspect_profiles: Optional[List[str]] = None,
    profile_loci_map: Optional[Dict[str, List[str]]] = None,
    include_suspect_alleles: bool = True,
    include_suspect_loci: bool = True,
    include_suspect_profiles: bool = True
) -> Dict[str, List[str]]:
    """Apply all filtering levels in order (most restrictive last).

    Filtering order:
    1. Allele-level (finest grain)
    2. Locus-level (removes all alleles from locus)
    3. Profile-level (most restrictive)

    Args:
        alleles_by_locus: Dict mapping locus to list of alleles
        suspect_alleles: List of suspect allele names
        suspect_loci: List of suspect locus names
        suspect_profiles: List of suspect profile ids
        profile_loci_map: Dict mapping profile id to list of loci
        include_suspect_alleles: If False, remove suspect alleles
        include_suspect_loci: If False, remove suspect loci
        include_suspect_profiles: If False, remove suspect profiles

    Returns:
        Filtered alleles by locus
    """
    result = alleles_by_locus.copy()

    # Apply allele-level filtering first (finest grain)
    if suspect_alleles:
        result = filter_alleles_by_locus(
            result,
            suspect_alleles,
            include_suspect=include_suspect_alleles
        )

    # Apply locus-level filtering (overrides allele level)
    if suspect_loci:
        result = filter_loci(
            result,
            suspect_loci,
            include_suspect=include_suspect_loci
        )

    # Apply profile-level filtering (most restrictive)
    if suspect_profiles and profile_loci_map:
        result = filter_profiles(
            result,
            suspect_profiles,
            profile_loci_map,
            include_suspect=include_suspect_profiles
        )

    return result


def prepare_allele_database(
    alleles_by_locus: Dict[str, List[str]],
    suspect_alleles: Optional[List[str]] = None,
    suspect_loci: Optional[List[str]] = None,
    suspect_profiles: Optional[List[str]] = None,
    profile_loci_map: Optional[Dict[str, List[str]]] = None,
    include_suspect_alleles: bool = True,
    include_suspect_loci: bool = True,
    include_suspect_profiles: bool = True
) -> Dict[str, List[str]]:
    """Prepare allele database by applying filters before MinHash.

    This is the main entry point for workflow filtering.

    Args:
        alleles_by_locus: Dict mapping locus to list of alleles
        suspect_alleles: List of suspect allele names
        suspect_loci: List of suspect locus names
        suspect_profiles: List of suspect profile ids
        profile_loci_map: Dict mapping profile id to list of loci
        include_suspect_alleles: If False, remove suspect alleles
        include_suspect_loci: If False, remove suspect loci
        include_suspect_profiles: If False, remove suspect profiles

    Returns:
        Prepared allele database ready for MinHash
    """
    return apply_filters(
        alleles_by_locus,
        suspect_alleles=suspect_alleles,
        suspect_loci=suspect_loci,
        suspect_profiles=suspect_profiles,
        profile_loci_map=profile_loci_map,
        include_suspect_alleles=include_suspect_alleles,
        include_suspect_loci=include_suspect_loci,
        include_suspect_profiles=include_suspect_profiles
    )


def create_filtered_fasta(
    input_path: Path,
    output_path: Path,
    suspect_alleles: List[str],
    include_suspect: bool = True
) -> None:
    """Create filtered FASTA by removing suspect allele sequences.

    Args:
        input_path: Path to input FASTA file
        output_path: Path to output FASTA file
        suspect_alleles: List of suspect allele names to exclude
        include_suspect: If False, exclude suspect alleles
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if include_suspect or not suspect_alleles:
        # No filtering needed, copy input to output
        with open(input_path) as f_in:
            with open(output_path, 'w') as f_out:
                f_out.write(f_in.read())
        return

    suspect_set = set(suspect_alleles)

    with open(input_path) as f_in:
        with open(output_path, 'w') as f_out:
            in_suspect = False
            for line in f_in:
                if line.startswith('>'):
                    # Extract allele name from header
                    header = line[1:].strip().split()[0]  # Take first token
                    in_suspect = header in suspect_set

                if not in_suspect:
                    f_out.write(line)


def filter_with_stats(
    alleles: List[str],
    suspect_alleles: List[str],
    include_suspect: bool = True
) -> Tuple[List[str], Dict[str, int]]:
    """Filter alleles and return statistics.

    Args:
        alleles: List of allele names
        suspect_alleles: List of suspect allele names
        include_suspect: If False, remove suspects

    Returns:
        Tuple of (filtered alleles, stats dict)
    """
    filtered = filter_alleles(alleles, suspect_alleles, include_suspect)

    stats = {
        'total_alleles': len(alleles),
        'suspect_alleles': len(suspect_alleles),
        'kept_alleles': len(filtered),
        'removed_alleles': len(alleles) - len(filtered)
    }

    return filtered, stats


def create_results_with_filter_info(
    results: Dict[str, Any],
    filter_info: Dict[str, Any]
) -> Dict[str, Any]:
    """Wrap results with filter metadata for reproducibility.

    Args:
        results: Original typing results dict
        filter_info: Dict with filtering information

    Returns:
        Dict with typing_results and filter_metadata sections
    """
    return {
        'typing_results': results,
        'filter_metadata': filter_info
    }


class WorkflowFilterConfig:
    """Configuration for workflow-level filtering.

    Tracks include/exclude decisions for all three filtering levels.
    Default: include all suspect data (no filtering).
    """

    def __init__(
        self,
        include_suspect_alleles: bool = True,
        include_suspect_loci: bool = True,
        include_suspect_profiles: bool = True
    ):
        """Initialize filter configuration.

        Args:
            include_suspect_alleles: If True (default), keep suspect alleles
            include_suspect_loci: If True (default), keep suspect loci
            include_suspect_profiles: If True (default), keep suspect profiles
        """
        self.include_suspect_alleles = include_suspect_alleles
        self.include_suspect_loci = include_suspect_loci
        self.include_suspect_profiles = include_suspect_profiles

    def should_filter_alleles(self) -> bool:
        """Check if allele-level filtering is enabled."""
        return not self.include_suspect_alleles

    def should_filter_loci(self) -> bool:
        """Check if locus-level filtering is enabled."""
        return not self.include_suspect_loci

    def should_filter_profiles(self) -> bool:
        """Check if profile-level filtering is enabled."""
        return not self.include_suspect_profiles

    def get_filtering_level(self) -> str:
        """Get the most restrictive filtering level active.

        Returns:
            'allele', 'locus', 'profile', or 'none'
        """
        if self.should_filter_profiles():
            return 'profile'
        elif self.should_filter_loci():
            return 'locus'
        elif self.should_filter_alleles():
            return 'allele'
        return 'none'
