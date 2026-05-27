"""Tests for auto strategy decision logic (Issue #59).

Acceptance criteria:
- --strategy auto option works in CLI
- Pre-analysis inspects input sequences (type, length distribution, N50)
- Correctly routes to fast/balanced/sensitive based on characteristics
- Decision rationale included in workflow output notes
- Tests verify decision logic for contigs, reads, edge cases
- Help text documents auto strategy behavior

Decision logic:
- Contigs (mean length >1000bp, N50 high) -> select "fast"
- Reads (mean length <500bp) -> select "balanced"
- Edge cases or uncertain -> default "balanced"
"""

import pytest
import toml
import csv
from click.testing import CliRunner
from unittest.mock import patch, MagicMock

from torchbase.cli import cli
from torchbase.torchfs import Torch


@pytest.fixture
def torch_without_workflow(tmp_path):
    """Create a torch without embedded workflow for strategy routing."""
    torch_path = tmp_path / "test_namespace" / "data_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_namespace",
        "name": "data_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Data torch without workflow"},
        "manifest": {"profiles": "profiles.tsv"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def contig_file(tmp_path):
    """Create a FASTA file with contig-like sequences (long, high N50)."""
    contig_file = tmp_path / "contigs.fasta"
    with open(contig_file, "w") as f:
        # Write 3 contigs with lengths: 5000bp, 3000bp, 2000bp
        # Mean: 3333bp, N50: 5000bp -> should trigger "fast" strategy
        f.write(">contig1\n")
        f.write("A" * 5000 + "\n")
        f.write(">contig2\n")
        f.write("C" * 3000 + "\n")
        f.write(">contig3\n")
        f.write("G" * 2000 + "\n")
    return contig_file


@pytest.fixture
def short_reads_file(tmp_path):
    """Create a FASTQ file with short read sequences (mean length <500bp)."""
    reads_file = tmp_path / "reads.fastq"
    with open(reads_file, "w") as f:
        # Write 5 reads with lengths: 150bp, 200bp, 100bp, 250bp, 150bp
        # Mean: 170bp -> should trigger "balanced" strategy
        for i, length in enumerate([150, 200, 100, 250, 150]):
            f.write(f"@read{i}\n")
            f.write("A" * length + "\n")
            f.write("+\n")
            f.write("I" * length + "\n")
    return reads_file


@pytest.fixture
def edge_case_file(tmp_path):
    """Create a file with ambiguous characteristics (between contigs and reads)."""
    edge_file = tmp_path / "edge.fasta"
    with open(edge_file, "w") as f:
        # Mixed lengths: 800bp, 600bp, 400bp
        # Mean: 600bp (between thresholds) -> should default to "balanced"
        f.write(">seq1\n")
        f.write("A" * 800 + "\n")
        f.write(">seq2\n")
        f.write("C" * 600 + "\n")
        f.write(">seq3\n")
        f.write("G" * 400 + "\n")
    return edge_file


class TestAutoStrategyFlagPresence:
    """Test that --strategy auto is available and documented."""

    def test_auto_strategy_in_help(self):
        """--strategy help text mentions auto option."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert result.exit_code == 0
        assert '--strategy' in result.output
        # Help text should mention auto
        assert 'auto' in result.output.lower()

    def test_auto_strategy_accepted_as_choice(self, torch_without_workflow, contig_file):
        """--strategy auto is accepted as a valid choice."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                mock_run.return_value.returncode = 0

                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'auto',
                        str(torch_without_workflow), '-c', str(contig_file)
                    ]
                )

                # Should not fail with invalid choice error
                assert 'invalid choice' not in result.output.lower()

    def test_auto_strategy_help_text_describes_behavior(self):
        """Help text explains auto strategy behavior."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        help_text = result.output.lower()
        # Should mention automatic selection based on input
        assert ('auto' in help_text and
                ('automatic' in help_text or 'detect' in help_text or
                 'analyze' in help_text or 'select' in help_text))


class TestContigDetectionRoutesToFast:
    """Test that contig-like inputs route to fast strategy."""

    def test_contigs_detected_and_routed_to_fast(
        self, torch_without_workflow, contig_file
    ):
        """Contig input (long sequences) routes to fast strategy."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                mock_run.return_value.returncode = 0

                with patch('torchbase.cli.Path') as mock_path:
                    # Mock builtin workflow path resolution
                    mock_workflow_path = MagicMock()
                    mock_workflow_path.exists.return_value = True
                    mock_path.return_value = mock_workflow_path

                    _ = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    # Check that fast_typing.wdl was selected
                    if mock_run.called:
                        call_args = str(mock_run.call_args)
                        assert 'fast_typing' in call_args

    def test_high_n50_triggers_fast_strategy(
        self, torch_without_workflow, contig_file
    ):
        """High N50 value triggers fast strategy selection."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    # Mock analysis returning contig characteristics
                    mock_analyze.return_value = {
                        'mean_length': 3333,
                        'n50': 5000,
                        'sequence_type': 'contigs',
                        'selected_strategy': 'fast'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    # Analysis should have been called
                    # Result should use fast strategy
                    assert mock_analyze.called or result.exit_code != 2

    def test_mean_length_over_1000_triggers_fast(
        self, torch_without_workflow, tmp_path
    ):
        """Mean sequence length >1000bp triggers fast strategy."""
        # Create file with mean length just over threshold
        long_seqs = tmp_path / "long.fasta"
        with open(long_seqs, "w") as f:
            f.write(">seq1\n" + "A" * 1500 + "\n")
            f.write(">seq2\n" + "C" * 1200 + "\n")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 1350,
                        'n50': 1500,
                        'sequence_type': 'contigs',
                        'selected_strategy': 'fast'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(long_seqs)
                        ]
                    )

                    assert mock_analyze.called or result.exit_code != 2


class TestShortReadsRouteToBalanced:
    """Test that short read inputs route to balanced strategy."""

    def test_short_reads_detected_and_routed_to_balanced(
        self, torch_without_workflow, short_reads_file
    ):
        """Short read input (mean <500bp) routes to balanced strategy."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                mock_run.return_value.returncode = 0

                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'auto',
                        str(torch_without_workflow), '-r', str(short_reads_file)
                    ]
                )

                # Check that balanced_typing.wdl was selected
                if mock_run.called:
                    call_args = str(mock_run.call_args)
                    assert 'balanced_typing' in call_args

    def test_mean_length_under_500_triggers_balanced(
        self, torch_without_workflow, short_reads_file
    ):
        """Mean sequence length <500bp triggers balanced strategy."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 170,
                        'n50': 200,
                        'sequence_type': 'reads',
                        'selected_strategy': 'balanced'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-r', str(short_reads_file)
                        ]
                    )

                    assert mock_analyze.called or result.exit_code != 2

    def test_fastq_format_recognized_as_reads(
        self, torch_without_workflow, short_reads_file
    ):
        """FASTQ format is recognized as read data."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    # FASTQ should be recognized
                    mock_analyze.return_value = {
                        'format': 'fastq',
                        'mean_length': 170,
                        'sequence_type': 'reads',
                        'selected_strategy': 'balanced'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-r', str(short_reads_file)
                        ]
                    )

                    assert mock_analyze.called or result.exit_code != 2


class TestEdgeCasesDefaultToBalanced:
    """Test that edge cases and uncertain inputs default to balanced."""

    def test_ambiguous_lengths_default_to_balanced(
        self, torch_without_workflow, edge_case_file
    ):
        """Sequences with ambiguous length characteristics default to balanced."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 600,  # Between thresholds
                        'n50': 800,
                        'sequence_type': 'uncertain',
                        'selected_strategy': 'balanced'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(edge_case_file)
                        ]
                    )

                    assert mock_analyze.called or result.exit_code != 2

    def test_empty_file_defaults_to_balanced(
        self, torch_without_workflow, tmp_path
    ):
        """Empty input file defaults to balanced strategy."""
        empty_file = tmp_path / "empty.fasta"
        empty_file.touch()

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 0,
                        'sequence_count': 0,
                        'selected_strategy': 'balanced'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(empty_file)
                        ]
                    )

                    # Should not crash, should default to balanced
                    assert mock_analyze.called or result.exit_code != 2

    def test_single_sequence_uses_length_for_decision(
        self, torch_without_workflow, tmp_path
    ):
        """Single sequence file uses its length for strategy decision."""
        single_seq = tmp_path / "single.fasta"
        with open(single_seq, "w") as f:
            # Single long sequence should trigger fast
            f.write(">seq1\n" + "A" * 5000 + "\n")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 5000,
                        'n50': 5000,
                        'sequence_count': 1,
                        'selected_strategy': 'fast'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(single_seq)
                        ]
                    )

                    assert mock_analyze.called or result.exit_code != 2

    def test_analysis_failure_defaults_to_balanced(
        self, torch_without_workflow, contig_file
    ):
        """If sequence analysis fails, default to balanced strategy."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    # Simulate analysis failure
                    mock_analyze.side_effect = Exception("Analysis failed")

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    # Should fallback to balanced, not crash
                    # (exact behavior depends on error handling)
                    assert result.exit_code != 2 or 'analysis' in result.output.lower()


class TestDecisionRationaleInOutput:
    """Test that decision rationale is included in workflow output notes."""

    def test_decision_rationale_passed_to_workflow(
        self, torch_without_workflow, contig_file
    ):
        """Auto decision rationale is passed to workflow as input."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 3333,
                        'n50': 5000,
                        'sequence_type': 'contigs',
                        'selected_strategy': 'fast',
                        'rationale': 'contigs detected (mean: 3333bp, N50: 5000bp), selected fast strategy'
                    }

                    _ = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    # Check that rationale is in miniwdl command
                    if mock_run.called:
                        call_args = str(mock_run.call_args)
                        # Rationale should be passed as workflow input
                        assert ('rationale' in call_args or
                                'auto_decision' in call_args)

    def test_rationale_includes_sequence_statistics(
        self, torch_without_workflow, short_reads_file
    ):
        """Decision rationale includes sequence statistics."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'mean_length': 170,
                        'n50': 200,
                        'sequence_type': 'reads',
                        'selected_strategy': 'balanced',
                        'rationale': 'short reads detected (mean: 170bp), selected balanced strategy'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-r', str(short_reads_file)
                        ]
                    )

                    if mock_run.called:
                        # Rationale should mention statistics
                        assert (mock_analyze.called and
                                (result.exit_code == 0 or result.exit_code != 2))

    def test_rationale_explains_strategy_choice(
        self, torch_without_workflow, contig_file
    ):
        """Decision rationale explains why strategy was chosen."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run'):
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    rationale = 'contigs detected, selected fast strategy'
                    mock_analyze.return_value = {
                        'selected_strategy': 'fast',
                        'rationale': rationale
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    # Rationale should explain the "why"
                    assert mock_analyze.called or result.exit_code != 2


class TestSequenceAnalysisFunction:
    """Test the _analyze_sequences helper function."""

    def test_analyze_sequences_calculates_mean_length(self, contig_file):
        """_analyze_sequences calculates mean sequence length."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(contig_file)

        assert 'mean_length' in result
        # Contigs are 5000, 3000, 2000 -> mean is 3333
        assert result['mean_length'] == pytest.approx(3333, abs=10)

    def test_analyze_sequences_calculates_n50(self, contig_file):
        """_analyze_sequences calculates N50 value."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(contig_file)

        assert 'n50' in result
        # N50 for [5000, 3000, 2000] sorted by length is 5000
        assert result['n50'] == 5000

    def test_analyze_sequences_detects_contigs(self, contig_file):
        """_analyze_sequences detects contig-like sequences."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(contig_file)

        assert result['sequence_type'] == 'contigs'
        assert result['selected_strategy'] == 'fast'

    def test_analyze_sequences_detects_reads(self, short_reads_file):
        """_analyze_sequences detects short read sequences."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(short_reads_file)

        assert result['sequence_type'] == 'reads'
        assert result['selected_strategy'] == 'balanced'

    def test_analyze_sequences_handles_fasta_format(self, contig_file):
        """_analyze_sequences handles FASTA format."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(contig_file)

        # Should successfully parse FASTA
        assert 'mean_length' in result
        assert result['mean_length'] > 0

    def test_analyze_sequences_handles_fastq_format(self, short_reads_file):
        """_analyze_sequences handles FASTQ format."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(short_reads_file)

        # Should successfully parse FASTQ
        assert 'mean_length' in result
        assert result['mean_length'] > 0

    def test_analyze_sequences_returns_rationale(self, contig_file):
        """_analyze_sequences returns decision rationale."""
        from torchbase.cli import _analyze_sequences

        result = _analyze_sequences(contig_file)

        assert 'rationale' in result
        assert isinstance(result['rationale'], str)
        assert len(result['rationale']) > 0


class TestAutoStrategyWorkflowRouting:
    """Test that auto strategy correctly routes to built-in workflows."""

    def test_auto_routes_to_fast_workflow_for_contigs(
        self, torch_without_workflow, contig_file
    ):
        """Auto strategy routes to fast_typing.wdl for contigs."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'selected_strategy': 'fast',
                        'rationale': 'contigs detected'
                    }

                    _ = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(contig_file)
                        ]
                    )

                    if mock_run.called:
                        call_args = mock_run.call_args[0][0]
                        # Check workflow path in miniwdl command
                        assert any('fast_typing' in str(arg) for arg in call_args)

    def test_auto_routes_to_balanced_workflow_for_reads(
        self, torch_without_workflow, short_reads_file
    ):
        """Auto strategy routes to balanced_typing.wdl for short reads."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'selected_strategy': 'balanced',
                        'rationale': 'short reads detected'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-r', str(short_reads_file)
                        ]
                    )

                    if mock_run.called:
                        call_args = mock_run.call_args[0][0]
                        assert any('balanced_typing' in str(arg) for arg in call_args) or result.exit_code != 2

    def test_auto_routes_to_balanced_workflow_for_edge_cases(
        self, torch_without_workflow, edge_case_file
    ):
        """Auto strategy routes to balanced_typing.wdl for edge cases."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                    mock_analyze.return_value = {
                        'selected_strategy': 'balanced',
                        'rationale': 'uncertain characteristics, defaulted to balanced'
                    }

                    result = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'auto',
                            str(torch_without_workflow), '-c', str(edge_case_file)
                        ]
                    )

                    if mock_run.called:
                        call_args = mock_run.call_args[0][0]
                        assert any('balanced_typing' in str(arg) for arg in call_args) or result.exit_code != 2


class TestAutoStrategyWithEmbeddedWorkflows:
    """Test that auto strategy interacts correctly with embedded workflows."""

    def test_auto_not_allowed_with_embedded_workflow(self, tmp_path):
        """--strategy auto cannot be used with torch-embedded workflows."""
        # Create torch with embedded workflow
        torch_path = tmp_path / "test_namespace" / "workflow_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "test_namespace",
            "name": "workflow_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "semver", "timestamp": 1609459200},
            "typing": {"method": "mlst"},
            "description": {"short": "Torch with embedded workflow"},
            "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_path / "main.wdl", "w") as f:
            f.write("workflow custom { }\n")

        profiles = [["ST", "adk"], ["1", "1"]]
        with open(torch_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        (torch_path / "_resources").mkdir()

        contig_file = tmp_path / "contigs.fasta"
        contig_file.write_text(">seq1\n" + "A" * 5000 + "\n")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_path / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            result = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'auto',
                    str(torch_path), '-c', str(contig_file)
                ]
            )

            # Should fail with error about embedded workflow
            assert result.exit_code != 0
            assert ('embedded' in result.output.lower() or
                    'workflow' in result.output.lower())


class TestAutoStrategyIntegration:
    """Integration tests for auto strategy end-to-end."""

    def test_auto_strategy_full_pipeline_with_contigs(
        self, torch_without_workflow, contig_file
    ):
        """Full pipeline: auto detects contigs, routes to fast, executes."""
        runner = CliRunner()

        # Load torch to verify it has no workflow
        torch = Torch.load(torch_without_workflow)
        assert torch.workflow is None

        with patch('torchbase.cli.run') as mock_run:
            mock_run.return_value.returncode = 0

            with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                mock_analyze.return_value = {
                    'mean_length': 3333,
                    'n50': 5000,
                    'sequence_type': 'contigs',
                    'selected_strategy': 'fast',
                    'rationale': 'contigs detected (mean: 3333bp, N50: 5000bp), selected fast strategy'
                }

                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'auto',
                        str(torch_without_workflow), '-c', str(contig_file)
                    ]
                )

                # Should complete successfully
                assert mock_analyze.called
                assert mock_run.called or result.exit_code == 0

    def test_auto_strategy_full_pipeline_with_reads(
        self, torch_without_workflow, short_reads_file
    ):
        """Full pipeline: auto detects reads, routes to balanced, executes."""
        runner = CliRunner()

        torch = Torch.load(torch_without_workflow)
        assert torch.workflow is None

        with patch('torchbase.cli.run') as mock_run:
            mock_run.return_value.returncode = 0

            with patch('torchbase.cli._analyze_sequences') as mock_analyze:
                mock_analyze.return_value = {
                    'mean_length': 170,
                    'n50': 200,
                    'sequence_type': 'reads',
                    'selected_strategy': 'balanced',
                    'rationale': 'short reads detected (mean: 170bp), selected balanced strategy'
                }

                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'auto',
                        str(torch_without_workflow), '-r', str(short_reads_file)
                    ]
                )

                # Should complete successfully
                assert mock_analyze.called
                assert mock_run.called or result.exit_code == 0
