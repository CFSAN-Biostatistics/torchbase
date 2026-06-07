"""Real minimap2 execution tests for alignment-based allele calling.

Tests are skipped when minimap2 is not installed. When minimap2 is available,
these tests execute the real binary and assert on actual alignment output.
"""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

minimap2_available = pytest.mark.skipif(
    shutil.which("minimap2") is None,
    reason="minimap2 not installed",
)

REFERENCE_FASTA = """>adk_1
ATGAAAATCATTAATAAAAAAATTGAAGAAAAAATTCCCAAAAAAGAATTTCAGCATCTGCAAGCAAAAGCGATCATCAACCGCCATCTGCAGTCTGCAGC
>adk_2
ATGAAAATCATTAATAAAAAAATTGAAGAAAAAATTCCCAAAAAAGAATTTCAGCATCTGCAAGCAAAAGCGATCATCAACCGCCATCTGCAGTCTGCAGG
"""

# Query matches adk_1 exactly
QUERY_FASTA = """>query_contig
ATGAAAATCATTAATAAAAAAATTGAAGAAAAAATTCCCAAAAAAGAATTTCAGCATCTGCAAGCAAAAGCGATCATCAACCGCCATCTGCAGTCTGCAGC
"""


class TestMinimap2Alignment:
    @minimap2_available
    def test_minimap2_produces_sam_output(self):
        """minimap2 runs and produces SAM format output."""
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp) / "ref.fasta"
            qry = Path(tmp) / "query.fasta"
            ref.write_text(REFERENCE_FASTA)
            qry.write_text(QUERY_FASTA)

            result = subprocess.run(
                ["minimap2", "-a", "--eqx", "-x", "asm5", str(ref), str(qry)],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr
            # SAM output has header lines (@) and alignment lines
            lines = [l for l in result.stdout.splitlines() if not l.startswith("@")]
            assert len(lines) >= 1, "Expected at least one alignment"

    @minimap2_available
    def test_exact_match_has_high_identity(self):
        """Exact sequence match should produce near-100% identity in alignment."""
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp) / "ref.fasta"
            qry = Path(tmp) / "query.fasta"
            ref.write_text(REFERENCE_FASTA)
            qry.write_text(QUERY_FASTA)

            result = subprocess.run(
                ["minimap2", "-a", "--eqx", "-x", "asm5", str(ref), str(qry)],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr

            # Parse NM tag (edit distance) from SAM
            nm = None
            for line in result.stdout.splitlines():
                if line.startswith("@"):
                    continue
                fields = line.split("\t")
                for tag in fields[11:]:
                    if tag.startswith("NM:i:"):
                        nm = int(tag.split(":")[2])
                        break
                if nm is not None:
                    break

            assert nm is not None, "NM tag not found in SAM output"
            assert nm == 0, f"Expected 0 edit distance for exact match, got {nm}"

    @minimap2_available
    def test_snp_variant_has_nonzero_edit_distance(self):
        """Sequence with SNP vs reference should have NM >= 1."""
        ref_fasta = ">allele_1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n"
        # Single SNP at position 5
        qry_fasta = ">query\nATGCGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n"

        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp) / "ref.fasta"
            qry = Path(tmp) / "query.fasta"
            ref.write_text(ref_fasta)
            qry.write_text(qry_fasta)

            result = subprocess.run(
                ["minimap2", "-a", "--eqx", "-x", "asm5", str(ref), str(qry)],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr

            nm = None
            for line in result.stdout.splitlines():
                if line.startswith("@"):
                    continue
                fields = line.split("\t")
                for tag in fields[11:]:
                    if tag.startswith("NM:i:"):
                        nm = int(tag.split(":")[2])
                        break
                if nm is not None:
                    break

            assert nm is not None
            assert nm >= 1, f"Expected NM >= 1 for SNP variant, got {nm}"
