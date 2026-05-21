"""Acceptance tests for MinHash allele calling WDL task (Issue #15).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- WDL task signature: query sequences + allele FASTA → allele calls JSON
- Depth filtering for reads: histogram-based k-mer filtering, fallback ≥3x
- MinHash sketching via sourmash
- Best match per locus with similarity score
- Output JSON: {locus: {allele_id, similarity, confidence}}
- Containerized (Docker/Singularity with sourmash)
- miniwdl check validates syntax
- Test with synthetic data (known allele combinations)
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def allele_database_fasta():
    """Create temporary allele database (FASTA format).

    Creates MLST-style allele database with multiple loci and alleles:
    - adk: 3 alleles
    - fumC: 2 alleles
    - gyrB: 2 alleles
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "allele_db.fasta"

        # Create realistic MLST allele sequences (shorter for testing)
        fasta_content = """>adk_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>adk_2
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>adk_3
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGATCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>fumC_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCTTCTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAA
>fumC_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
>gyrB_2
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(db_path, "w") as f:
            f.write(fasta_content)

        yield db_path


@pytest.fixture
def query_contigs_fasta():
    """Create temporary query contigs file (FASTA format).

    Simulates contigs matching adk_1, fumC_2, gyrB_1.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "query_contigs.fasta"

        # Create contigs that match specific alleles
        fasta_content = """>contig1_contains_adk_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>contig2_contains_fumC_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>contig3_contains_gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


@pytest.fixture
def query_reads_fasta():
    """Create temporary query reads file (FASTA format).

    Simulates reads matching adk_1, fumC_2, gyrB_1 at varying depths.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        reads_path = tmpdir_path / "query_reads.fasta"

        # Create reads that match specific alleles (subset of sequences)
        fasta_content = """>read1_adk_1_depth1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGAC
>read2_adk_1_depth2
TGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTG
>read3_adk_1_depth3
GTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACAT
>read4_fumC_2_depth1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>read5_fumC_2_depth2
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>read6_gyrB_1_depth1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
"""

        with open(reads_path, "w") as f:
            f.write(fasta_content)

        yield reads_path


class TestMinHashWDLFileExists:
    """Test MinHash WDL file exists at expected location."""

    def test_wdl_file_exists(self):
        """MinHash allele calling WDL file exists"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        assert wdl_path.exists(), f"WDL file not found at {wdl_path}"

    def test_wdl_file_is_file(self):
        """MinHash WDL file is a regular file"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        assert wdl_path.is_file(), f"WDL path is not a file: {wdl_path}"


class TestMinHashWDLTaskSignature:
    """Test WDL task signature: query sequences + allele FASTA → allele calls JSON."""

    def test_wdl_has_task_definition(self):
        """WDL file contains task definition"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "task" in content, "WDL file does not contain task definition"

    def test_wdl_task_name_is_minhash_allele_calling(self):
        """WDL workflow/task is named minhash_allele_calling"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Accept either task or workflow with minhash_allele_calling name
        assert ("task minhash_allele_calling" in content or
                "workflow minhash_allele_calling" in content), "Task/workflow name is incorrect"

    def test_wdl_has_input_section(self):
        """WDL task has input section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "input {" in content, "Task does not have input section"

    def test_wdl_has_query_sequences_input(self):
        """WDL task accepts query_sequences input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "query_sequences" in content, "Task does not have query_sequences input"

    def test_wdl_has_allele_fasta_input(self):
        """WDL task accepts allele_fasta input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "allele_fasta" in content, "Task does not have allele_fasta input"

    def test_wdl_has_output_section(self):
        """WDL task has output section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "output {" in content, "Task does not have output section"

    def test_wdl_output_is_file(self):
        """WDL task output is a File type"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for File type in output section
        assert "File" in content, "Task output is not File type"


class TestMinHashWDLSyntaxValidation:
    """Test miniwdl check validates WDL syntax."""

    def test_miniwdl_check_passes(self):
        """miniwdl check validates WDL syntax without errors"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        # Run miniwdl check
        result = subprocess.run(
            ["miniwdl", "check", str(wdl_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed: {result.stderr}"


class TestMinHashWDLContainerization:
    """Test containerized (Docker/Singularity with sourmash)."""

    def test_wdl_has_runtime_section(self):
        """WDL task has runtime section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "runtime {" in content, "Task does not have runtime section"

    def test_wdl_has_docker_image(self):
        """WDL task specifies docker image"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "docker:" in content, "Task does not specify docker image"

    def test_wdl_uses_sourmash_container(self):
        """WDL task uses sourmash container image"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "sourmash" in content.lower(), "Task does not use sourmash container"


class TestMinHashWDLDepthFiltering:
    """Test depth filtering for reads: histogram-based k-mer filtering, fallback ≥3x."""

    def test_wdl_has_reads_file_input(self):
        """WDL task accepts optional reads_file input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for optional reads input
        assert "reads_file" in content or "reads" in content, "Task does not accept reads input"

    def test_wdl_has_min_coverage_parameter(self):
        """WDL task has min_coverage parameter for depth filtering"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "min_coverage" in content or "coverage" in content, "Task does not have coverage parameter"

    def test_wdl_command_includes_depth_filtering_logic(self):
        """WDL task command includes depth filtering logic"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for depth/coverage/histogram/filter keywords
        assert any(keyword in content.lower() for keyword in ["depth", "coverage", "histogram", "filter"]), \
            "Task command does not include depth filtering logic"


class TestMinHashWDLOutputFormat:
    """Test output JSON: {locus: {allele_id, similarity, confidence}}."""

    def test_wdl_produces_json_output(self):
        """WDL task produces JSON output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert ".json" in content.lower() or "json" in content.lower(), "Task does not produce JSON output"

    def test_wdl_command_writes_json_file(self):
        """WDL task command writes JSON output file"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for JSON writing logic
        assert "json.dump" in content or "json.dumps" in content or "> " in content, \
            "Task command does not write JSON file"


class TestMinHashWDLSourmashIntegration:
    """Test MinHash sketching via sourmash."""

    def test_wdl_command_uses_sourmash(self):
        """WDL task command uses sourmash"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "sourmash" in content, "Task command does not use sourmash"

    def test_wdl_command_creates_sketches(self):
        """WDL task command creates MinHash sketches"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for sketch creation
        assert "sketch" in content.lower(), "Task command does not create sketches"

    def test_wdl_command_compares_sketches(self):
        """WDL task command compares sketches for similarity"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for comparison logic
        assert "compare" in content.lower() or "similarity" in content.lower(), \
            "Task command does not compare sketches"


@pytest.mark.miniwdl
class TestMinHashWDLExecutionWithContigs:
    """Test with synthetic data (known allele combinations) - contigs."""

    def test_wdl_execution_with_contigs_produces_output(self, allele_database_fasta, query_contigs_fasta):
        """WDL task execution with contigs produces output JSON"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create input JSON for miniwdl
            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            # Run miniwdl
            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Check execution succeeded
            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            # Find output directory
            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

            # Read outputs
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            # Check output file exists
            assert "minhash_allele_calling.results" in outputs, "Output results not found"

            results_path = Path(outputs["minhash_allele_calling.results"])
            assert results_path.exists(), f"Results file does not exist: {results_path}"

            # Read and validate JSON structure
            with open(results_path) as f:
                results = json.load(f)

            # Check expected loci are present
            assert "adk" in results, "adk locus not in results"
            assert "fumC" in results, "fumC locus not in results"
            assert "gyrB" in results, "gyrB locus not in results"

    def test_wdl_output_has_allele_id_field(self, allele_database_fasta, query_contigs_fasta):
        """WDL output JSON contains allele_id field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # Check allele_id field
            for locus in results:
                assert "allele_id" in results[locus], f"allele_id not in {locus} results"

    def test_wdl_output_has_similarity_field(self, allele_database_fasta, query_contigs_fasta):
        """WDL output JSON contains similarity field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # Check similarity field
            for locus in results:
                assert "similarity" in results[locus], f"similarity not in {locus} results"
                assert 0.0 <= results[locus]["similarity"] <= 1.0, f"similarity out of range for {locus}"

    def test_wdl_output_has_confidence_field(self, allele_database_fasta, query_contigs_fasta):
        """WDL output JSON contains confidence field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # Check confidence field
            for locus in results:
                assert "confidence" in results[locus], f"confidence not in {locus} results"
                assert isinstance(results[locus]["confidence"], bool), f"confidence is not boolean for {locus}"


@pytest.mark.miniwdl
class TestMinHashWDLExecutionWithReads:
    """Test with synthetic data (known allele combinations) - reads with depth filtering."""

    def test_wdl_execution_with_reads_produces_output(self, allele_database_fasta, query_reads_fasta):
        """WDL task execution with reads produces output JSON"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create input JSON with reads
            input_json = {
                "minhash_allele_calling.query_sequences": str(query_reads_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta),
                "minhash_allele_calling.reads_file": str(query_reads_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            # Run miniwdl
            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Check execution succeeded
            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            # Find output directory
            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

            # Read outputs
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            # Check output file exists
            assert "minhash_allele_calling.results" in outputs, "Output results not found"

            results_path = Path(outputs["minhash_allele_calling.results"])
            assert results_path.exists(), f"Results file does not exist: {results_path}"

    def test_wdl_reads_mode_applies_depth_filtering(self, allele_database_fasta, query_reads_fasta):
        """WDL task with reads applies depth filtering"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Test with high min_coverage to trigger filtering
            input_json = {
                "minhash_allele_calling.query_sequences": str(query_reads_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta),
                "minhash_allele_calling.reads_file": str(query_reads_fasta),
                "minhash_allele_calling.min_coverage": 3
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Check execution succeeded
            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"


@pytest.mark.miniwdl
class TestMinHashWDLCorrectAlleleIdentification:
    """Test best match per locus with similarity score."""

    def test_wdl_identifies_correct_alleles_from_contigs(self, allele_database_fasta, query_contigs_fasta):
        """WDL task correctly identifies alleles from contigs"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # Check expected alleles (based on query_contigs_fasta fixture)
            # contig1 contains adk_1
            assert results["adk"]["allele_id"] == "1", f"Expected adk_1, got {results['adk']['allele_id']}"

            # contig2 contains fumC_2
            assert results["fumC"]["allele_id"] == "2", f"Expected fumC_2, got {results['fumC']['allele_id']}"

            # contig3 contains gyrB_1
            assert results["gyrB"]["allele_id"] == "1", f"Expected gyrB_1, got {results['gyrB']['allele_id']}"

    def test_wdl_high_similarity_for_exact_matches(self, allele_database_fasta, query_contigs_fasta):
        """WDL task reports high similarity for exact matches"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # For exact matches, similarity should be very high (> 0.9)
            for locus in results:
                assert results[locus]["similarity"] > 0.9, \
                    f"Similarity too low for exact match in {locus}: {results[locus]['similarity']}"

    def test_wdl_high_confidence_for_exact_matches(self, allele_database_fasta, query_contigs_fasta):
        """WDL task reports high confidence for exact matches"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "minhash_allele_calling.query_sequences": str(query_contigs_fasta),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # For exact matches, confidence should be True
            for locus in results:
                assert results[locus]["confidence"] is True, \
                    f"Confidence should be True for exact match in {locus}"


@pytest.mark.miniwdl
class TestMinHashWDLEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_wdl_handles_empty_query(self, allele_database_fasta):
        """WDL task handles empty query file gracefully"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create empty query file
            empty_query = tmpdir_path / "empty.fasta"
            empty_query.touch()

            input_json = {
                "minhash_allele_calling.query_sequences": str(empty_query),
                "minhash_allele_calling.allele_fasta": str(allele_database_fasta)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Should either succeed with empty results or fail gracefully
            # (implementation choice - document expected behavior)
            assert result.returncode in [0, 1], "Unexpected return code for empty query"

    def test_wdl_handles_single_locus(self):
        """WDL task handles allele database with single locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create single-locus database
            single_locus_db = tmpdir_path / "single_locus.fasta"
            with open(single_locus_db, "w") as f:
                f.write(">adk_1\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")
                f.write(">adk_2\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")

            # Create matching query
            query = tmpdir_path / "query.fasta"
            with open(query, "w") as f:
                f.write(">query\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")

            input_json = {
                "minhash_allele_calling.query_sequences": str(query),
                "minhash_allele_calling.allele_fasta": str(single_locus_db)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["minhash_allele_calling.results"])
            with open(results_path) as f:
                results = json.load(f)

            # Should have exactly one locus
            assert len(results) == 1, f"Expected 1 locus, got {len(results)}"
            assert "adk" in results, "Expected adk locus in results"
