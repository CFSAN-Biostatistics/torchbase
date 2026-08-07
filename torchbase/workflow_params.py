"""Workflow parameter parsing and auto-provisioning for embedded WDL workflows."""

from pathlib import Path, PurePath
from typing import Dict, List, Optional, Tuple


# Parameter name patterns for torch data files
TORCH_PARAM_PATTERNS = {
    'allele_fasta': 'alleles',
    'allele_database': 'alleles',
    'profiles_table': 'profiles',
    'profiles': 'profiles',
}


def parse_torch_args(args: List[str]) -> Dict[str, str]:
    """Parse key=value pairs from CLI torch_args.

    Args:
        args: List of strings in format ['key=value', 'key2=value2']

    Returns:
        Dict mapping parameter names to values

    Raises:
        ValueError: If any argument is not in key=value format

    Examples:
        >>> parse_torch_args(['confidence_threshold=0.95', 'min_depth=5'])
        {'confidence_threshold': '0.95', 'min_depth': '5'}
    """
    params = {}

    for arg in args:
        if '=' not in arg:
            raise ValueError(
                f"Invalid parameter format: '{arg}'. "
                f"Parameters must be in key=value format."
            )

        key, value = arg.split('=', 1)  # Split on first = only
        params[key] = value

    return params


def auto_provision_torch_parameters(
    workflow_inputs: Dict[str, Tuple[str, Optional[str]]],
    allele_fasta_path: Path,
    profiles_table_path: Path
) -> Dict[str, str]:
    """Auto-provision torch data files for workflow parameters.

    Matches workflow parameter names to torch data files using heuristic patterns.
    Only provisions File-type parameters that match known patterns.

    Args:
        workflow_inputs: Dict from WDLParser.workflow_inputs mapping
            param_name -> (type, default_value)
        allele_fasta_path: Path to concatenated allele FASTA
        profiles_table_path: Path to profiles TSV

    Returns:
        Dict mapping parameter names to file paths as platform-neutral
        POSIX-style strings (WDL inputs JSON is consumed by engines that
        expect ``/`` separators regardless of host OS).

    Examples:
        >>> workflow_inputs = {
        ...     'allele_fasta': ('File', None),
        ...     'profiles_table': ('File', None),
        ...     'query_sequences': ('File', None),
        ...     'confidence_threshold': ('Float', '0.85')
        ... }
        >>> auto_provision_torch_parameters(workflow_inputs,
        ...     Path('/tmp/alleles.fasta'), Path('/tmp/profiles.tsv'))
        {'allele_fasta': '/tmp/alleles.fasta', 'profiles_table': '/tmp/profiles.tsv'}
    """
    provisioned = {}

    for param_name, (param_type, _) in workflow_inputs.items():
        # Skip non-File parameters
        if not param_type.startswith('File'):
            continue

        # Skip query_sequences (user provides via -c/-r flags)
        if param_name == 'query_sequences':
            continue

        # Check if parameter name matches a known pattern
        if param_name in TORCH_PARAM_PATTERNS:
            pattern = TORCH_PARAM_PATTERNS[param_name]

            if pattern == 'alleles':
                provisioned[param_name] = PurePath(allele_fasta_path).as_posix()
            elif pattern == 'profiles':
                provisioned[param_name] = PurePath(profiles_table_path).as_posix()

    return provisioned


def validate_required_parameters(
    workflow_inputs: Dict[str, Tuple[str, Optional[str]]],
    provided_params: Dict[str, str]
) -> List[str]:
    """Validate that all required workflow parameters are provided.

    A parameter is required if:
    - It has no default value (default is None)
    - It is not nullable (type doesn't contain '?')

    Args:
        workflow_inputs: Dict from WDLParser.workflow_inputs
        provided_params: Dict of all parameters (auto + user)

    Returns:
        List of error messages (empty if all required params present)

    Examples:
        >>> workflow_inputs = {
        ...     'required_file': ('File', None),
        ...     'optional_file': ('File?', None),
        ...     'default_string': ('String', '"default"')
        ... }
        >>> validate_required_parameters(workflow_inputs, {})
        ['Missing required parameter: required_file (type: File)']
        >>> validate_required_parameters(workflow_inputs, {'required_file': 'x'})
        []
    """
    errors = []

    for param_name, (param_type, default_value) in workflow_inputs.items():
        # Has default value -> not required
        if default_value is not None:
            continue

        # Is nullable (File? or String?) -> not required
        if '?' in param_type:
            continue

        # Check if provided
        if param_name not in provided_params:
            errors.append(
                f"Missing required parameter: {param_name} (type: {param_type})"
            )

    return errors
