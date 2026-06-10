"""Tests for workflow parameter parsing and auto-provisioning."""

import pytest
from pathlib import Path
from torchbase.workflow_params import (
    parse_torch_args,
    auto_provision_torch_parameters,
    validate_required_parameters,
)


class TestParseTorchArgs:
    """Tests for parse_torch_args function."""

    def test_parse_single_parameter(self):
        """Parse single key=value pair."""
        result = parse_torch_args(['key=value'])
        assert result == {'key': 'value'}

    def test_parse_multiple_parameters(self):
        """Parse multiple key=value pairs."""
        result = parse_torch_args(['a=1', 'b=2', 'c=3'])
        assert result == {'a': '1', 'b': '2', 'c': '3'}

    def test_parse_numeric_values(self):
        """Parse numeric values as strings."""
        result = parse_torch_args([
            'confidence_threshold=0.95',
            'min_depth=5',
            'max_depth=100'
        ])
        assert result == {
            'confidence_threshold': '0.95',
            'min_depth': '5',
            'max_depth': '100'
        }

    def test_parse_boolean_values(self):
        """Parse boolean values as strings."""
        result = parse_torch_args([
            'exclude_suspect_alleles=true',
            'exclude_suspect_loci=false'
        ])
        assert result == {
            'exclude_suspect_alleles': 'true',
            'exclude_suspect_loci': 'false'
        }

    def test_parse_file_paths(self):
        """Parse file paths with slashes."""
        result = parse_torch_args(['path=/home/user/file.txt'])
        assert result == {'path': '/home/user/file.txt'}

    def test_parse_value_with_equals(self):
        """Parse value containing equals sign (split on first = only)."""
        result = parse_torch_args(['query=field=value'])
        assert result == {'query': 'field=value'}

    def test_parse_empty_list(self):
        """Parse empty argument list."""
        result = parse_torch_args([])
        assert result == {}

    def test_parse_empty_value(self):
        """Parse parameter with empty value."""
        result = parse_torch_args(['key='])
        assert result == {'key': ''}

    def test_invalid_format_no_equals(self):
        """Raise error when argument has no equals sign."""
        with pytest.raises(ValueError, match="Invalid parameter format"):
            parse_torch_args(['no_equals_sign'])

    def test_invalid_format_mixed(self):
        """Raise error when some arguments are invalid."""
        with pytest.raises(ValueError, match="Invalid parameter format"):
            parse_torch_args(['valid=value', 'invalid'])


class TestAutoProvisionTorchParameters:
    """Tests for auto_provision_torch_parameters function."""

    def test_provision_allele_fasta(self):
        """Auto-provision allele_fasta parameter."""
        workflow_inputs = {
            'allele_fasta': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {'allele_fasta': '/tmp/alleles.fasta'}

    def test_provision_allele_database(self):
        """Auto-provision allele_database parameter (alternative name)."""
        workflow_inputs = {
            'allele_database': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {'allele_database': '/tmp/alleles.fasta'}

    def test_provision_profiles_table(self):
        """Auto-provision profiles_table parameter."""
        workflow_inputs = {
            'profiles_table': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {'profiles_table': '/tmp/profiles.tsv'}

    def test_provision_profiles(self):
        """Auto-provision profiles parameter (alternative name)."""
        workflow_inputs = {
            'profiles': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {'profiles': '/tmp/profiles.tsv'}

    def test_provision_multiple_parameters(self):
        """Auto-provision multiple torch parameters."""
        workflow_inputs = {
            'allele_fasta': ('File', None),
            'profiles_table': ('File', None),
            'confidence_threshold': ('Float', '0.85'),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {
            'allele_fasta': '/tmp/alleles.fasta',
            'profiles_table': '/tmp/profiles.tsv',
        }
        # confidence_threshold should NOT be provisioned (not a File)

    def test_skip_query_sequences(self):
        """Do not auto-provision query_sequences (user provides via -c/-r)."""
        workflow_inputs = {
            'query_sequences': ('File', None),
            'allele_fasta': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert 'query_sequences' not in result
        assert result == {'allele_fasta': '/tmp/alleles.fasta'}

    def test_skip_non_file_parameters(self):
        """Do not auto-provision non-File parameters."""
        workflow_inputs = {
            'confidence_threshold': ('Float', '0.85'),
            'min_depth': ('Int', '3'),
            'input_type': ('String', '"contigs"'),
            'exclude_suspect_alleles': ('Boolean', 'false'),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {}

    def test_skip_unknown_file_parameters(self):
        """Do not auto-provision File parameters with unknown names."""
        workflow_inputs = {
            'unknown_file': ('File', None),
            'custom_input': ('File', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {}

    def test_nullable_file_parameter(self):
        """Auto-provision nullable File? parameters."""
        workflow_inputs = {
            'allele_fasta': ('File?', None),
        }
        result = auto_provision_torch_parameters(
            workflow_inputs,
            Path('/tmp/alleles.fasta'),
            Path('/tmp/profiles.tsv')
        )
        assert result == {'allele_fasta': '/tmp/alleles.fasta'}


class TestValidateRequiredParameters:
    """Tests for validate_required_parameters function."""

    def test_all_required_provided(self):
        """No errors when all required parameters are provided."""
        workflow_inputs = {
            'required_param': ('File', None),
            'another_required': ('String', None),
        }
        provided = {
            'required_param': 'value1',
            'another_required': 'value2',
        }
        errors = validate_required_parameters(workflow_inputs, provided)
        assert errors == []

    def test_missing_required_parameter(self):
        """Error when required parameter is missing."""
        workflow_inputs = {
            'required_param': ('File', None),
        }
        provided = {}
        errors = validate_required_parameters(workflow_inputs, provided)
        assert len(errors) == 1
        assert 'required_param' in errors[0]
        assert 'File' in errors[0]

    def test_missing_multiple_required(self):
        """Errors for multiple missing required parameters."""
        workflow_inputs = {
            'param1': ('File', None),
            'param2': ('String', None),
            'param3': ('Int', None),
        }
        provided = {}
        errors = validate_required_parameters(workflow_inputs, provided)
        assert len(errors) == 3

    def test_optional_nullable_not_required(self):
        """Nullable parameters (File?, String?) are not required."""
        workflow_inputs = {
            'optional_file': ('File?', None),
            'optional_string': ('String?', None),
        }
        provided = {}
        errors = validate_required_parameters(workflow_inputs, provided)
        assert errors == []

    def test_parameters_with_defaults_not_required(self):
        """Parameters with defaults are not required."""
        workflow_inputs = {
            'confidence_threshold': ('Float', '0.85'),
            'min_depth': ('Int', '3'),
            'input_type': ('String', '"contigs"'),
        }
        provided = {}
        errors = validate_required_parameters(workflow_inputs, provided)
        assert errors == []

    def test_mixed_required_and_optional(self):
        """Only required parameters trigger errors."""
        workflow_inputs = {
            'required': ('File', None),
            'optional_nullable': ('File?', None),
            'optional_default': ('String', '"default"'),
        }
        provided = {}
        errors = validate_required_parameters(workflow_inputs, provided)
        assert len(errors) == 1
        assert 'required' in errors[0]

    def test_extra_parameters_allowed(self):
        """No errors when extra parameters are provided."""
        workflow_inputs = {
            'required': ('File', None),
        }
        provided = {
            'required': 'value',
            'extra1': 'value1',
            'extra2': 'value2',
        }
        errors = validate_required_parameters(workflow_inputs, provided)
        assert errors == []

    def test_partial_provision(self):
        """Error only for missing required parameters."""
        workflow_inputs = {
            'param1': ('File', None),
            'param2': ('String', None),
            'param3': ('Int', None),
        }
        provided = {
            'param1': 'value1',
            'param3': 'value3',
        }
        errors = validate_required_parameters(workflow_inputs, provided)
        assert len(errors) == 1
        assert 'param2' in errors[0]
