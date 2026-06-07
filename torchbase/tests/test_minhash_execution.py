"""Real sourmash execution tests for MinHash sketching and comparison.

Tests are skipped when sourmash is not installed. When sourmash is available,
these tests execute the real binary and assert on actual output values.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

sourmash_available = pytest.mark.skipif(
    shutil.which("sourmash") is None,
    reason="sourmash not installed",
)

SAMPLE_FASTA = """>seq1_1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq1_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""


class TestSourmashSketch:
    @sourmash_available
    def test_sketch_produces_signature_file(self):
        """sketch_sequences equivalent: sourmash sketch dna produces a .sig file."""
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "test.fasta"
            fasta.write_text(SAMPLE_FASTA)
            sig = Path(tmp) / "test.sig"

            result = subprocess.run(
                ["sourmash", "sketch", "dna", "-p", "k=21,scaled=100",
                 str(fasta), "-o", str(sig), "--singleton"],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr
            assert sig.exists()
            assert sig.stat().st_size > 0

    @sourmash_available
    def test_sketch_self_similarity_is_one(self):
        """A sequence compared to itself should have ANI similarity of 1.0."""
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "test.fasta"
            fasta.write_text(SAMPLE_FASTA)
            sig = Path(tmp) / "test.sig"
            csv_out = Path(tmp) / "similarity.csv"

            subprocess.run(
                ["sourmash", "sketch", "dna", "-p", "k=21,scaled=100",
                 str(fasta), "-o", str(sig), "--singleton"],
                capture_output=True, check=True
            )
            result = subprocess.run(
                ["sourmash", "compare", str(sig), str(sig),
                 "--csv", str(csv_out), "--ani", "-k", "21"],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr
            assert csv_out.exists()

            import csv
            with open(csv_out) as f:
                rows = list(csv.reader(f))
            # Diagonal (self-similarity) should be 1.0
            assert len(rows) >= 2
            data_row = rows[1]
            self_sim = float(data_row[0])
            assert self_sim == pytest.approx(1.0, abs=0.01)

    @sourmash_available
    def test_different_sequences_have_low_similarity(self):
        """Two unrelated sequences should have low ANI similarity."""
        fasta1_content = ">seq1\n" + "ATGC" * 100 + "\n"
        fasta2_content = ">seq2\n" + "CCCC" * 100 + "\n"

        with tempfile.TemporaryDirectory() as tmp:
            f1 = Path(tmp) / "seq1.fasta"
            f2 = Path(tmp) / "seq2.fasta"
            f1.write_text(fasta1_content)
            f2.write_text(fasta2_content)
            sig1 = Path(tmp) / "sig1.sig"
            sig2 = Path(tmp) / "sig2.sig"
            csv_out = Path(tmp) / "sim.csv"

            subprocess.run(
                ["sourmash", "sketch", "dna", "-p", "k=21,scaled=100",
                 str(f1), "-o", str(sig1)],
                capture_output=True, check=True
            )
            subprocess.run(
                ["sourmash", "sketch", "dna", "-p", "k=21,scaled=100",
                 str(f2), "-o", str(sig2)],
                capture_output=True, check=True
            )
            result = subprocess.run(
                ["sourmash", "compare", str(sig1), str(sig2),
                 "--csv", str(csv_out), "--ani", "-k", "21"],
                capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr
            import csv
            with open(csv_out) as f:
                rows = list(csv.reader(f))
            # Off-diagonal similarity should be < 0.9
            assert len(rows) >= 2
            sim = float(rows[1][1]) if len(rows[1]) > 1 else float(rows[2][1])
            assert sim < 0.9
