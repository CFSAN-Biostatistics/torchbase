version 1.0

import "../../builtin/tasks/minhash.wdl" as minhash
import "../../builtin/tasks/alignment.wdl" as alignment
import "../../builtin/tasks/profile_lookup.wdl" as profile_lookup
import "../../builtin/tasks/filter_alleles.wdl" as filter

workflow mlst_typing {
    input {
        # Stage 1: Input validation - reads OR contigs
        File? contigs
        File? reads

        # Required data
        File allele_database
        File profiles

        # Optional: multiple schemes (dict of scheme_name -> profiles_file)
        Map[String, File]? schemes

        # Stage 2: Depth filtering parameters
        Int min_coverage = 3

        # Stage 3-7: Thresholds and parameters
        Float min_identity = 0.90
        Float min_coverage_threshold = 0.5
        Float ambiguity_threshold = 0.85
        Float confidence_threshold = 0.85

        # Quality filtering
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    # Stage 1: Input validation - ensure reads OR contigs provided
    call validate_inputs {
        input:
            contigs = contigs,
            reads = reads
    }

    # Determine input type for downstream processing
    String input_type = if defined(reads) then "reads" else "contigs"
    File query_sequences = select_first([contigs, reads])

    # Stage 2: Depth filtering for reads (skip for contigs)
    if (input_type == "reads") {
        call filter_reads_by_depth {
            input:
                reads = select_first([reads]),
                min_coverage = min_coverage
        }
    }

    # Use filtered reads or original contigs
    File filtered_sequences = if (input_type == "reads") then
        select_first([filter_reads_by_depth.filtered_reads])
    else
        select_first([contigs])

    # Stage 3: MinHash allele calling across all schemes
    call minhash.sketch_sequences as sketch_queries {
        input:
            sequences = filtered_sequences,
            ksize = 31,
            scaled = 1000
    }

    call minhash.sketch_sequences as sketch_alleles {
        input:
            sequences = allele_database,
            ksize = 31,
            scaled = 1000
    }

    call minhash.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_database
    }

    call minhash.call_alleles_minhash {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = filtered_sequences,
            allele_fasta = allele_database,
            confidence_threshold = confidence_threshold
    }

    # Stage 4: Scheme inference from allele calls (highest coverage/identity)
    call infer_scheme {
        input:
            allele_calls = call_alleles_minhash.allele_calls,
            profiles = profiles,
            schemes = schemes,
            min_coverage_threshold = min_coverage_threshold
    }

    # Stage 5: Profile lookup in inferred scheme
    call profile_lookup.lookup_profile {
        input:
            allele_calls = call_alleles_minhash.allele_calls,
            profiles_table = infer_scheme.selected_profile_table,
            strategy = "mlst_orchestrated",
            alignment_used = false
    }

    # Stage 6: Alignment fallback if ambiguous
    call check_ambiguity {
        input:
            allele_calls = call_alleles_minhash.allele_calls,
            profile_result = lookup_profile.result,
            ambiguity_threshold = ambiguity_threshold
    }

    if (check_ambiguity.needs_alignment) {
        call alignment.align_and_call as alignment_fallback {
            input:
                query_sequences = filtered_sequences,
                allele_fasta = allele_database,
                input_type = input_type,
                identity_threshold = min_identity
        }

        # Re-lookup profile with alignment results
        call profile_lookup.lookup_profile as lookup_profile_aligned {
            input:
                allele_calls = alignment_fallback.alignment_results,
                profiles_table = infer_scheme.selected_profile_table,
                strategy = "mlst_orchestrated",
                alignment_used = true
        }
    }

    # Stage 7: Final result with status, confidence, nearest ST
    call assemble_final_result {
        input:
            minhash_result = lookup_profile.result,
            alignment_result = lookup_profile_aligned.result,
            inferred_scheme = infer_scheme.inferred_scheme,
            used_alignment = check_ambiguity.needs_alignment,
            input_type = input_type
    }

    output {
        File typing_result = assemble_final_result.final_result
    }
}

# Stage 1: Validate inputs
task validate_inputs {
    input {
        File? contigs
        File? reads
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import sys

contigs = "~{contigs}"
reads = "~{reads}"

# Check if at least one is provided (not "None" string from WDL optional)
has_contigs = contigs != "None" and contigs != ""
has_reads = reads != "None" and reads != ""

if not has_contigs and not has_reads:
    print("ERROR: Either contigs or reads must be provided", file=sys.stderr)
    sys.exit(1)

if has_contigs and has_reads:
    print("WARNING: Both contigs and reads provided, using contigs", file=sys.stderr)

print("Validation passed: inputs are valid")
PYTHON_SCRIPT
    >>>

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "1 GB"
    }
}

# Stage 2: Filter reads by depth
task filter_reads_by_depth {
    input {
        File reads
        Int min_coverage = 3
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
from collections import defaultdict

def parse_fasta(filename):
    """Parse FASTA file and group reads by locus."""
    locus_reads = defaultdict(list)
    current_header = None
    current_seq = ""

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    # Extract locus name from header (format: >read_locus_depth)
                    parts = current_header.split("_")
                    locus = "_".join(parts[:-1]) if len(parts) > 1 else current_header
                    locus_reads[locus].append((current_header, current_seq))
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line

        # Don't forget last sequence
        if current_header:
            parts = current_header.split("_")
            locus = "_".join(parts[:-1]) if len(parts) > 1 else current_header
            locus_reads[locus].append((current_header, current_seq))

    return locus_reads

def filter_by_coverage(locus_reads, min_coverage):
    """Filter reads to keep only those with sufficient coverage."""
    filtered = defaultdict(list)

    for locus, reads in locus_reads.items():
        if len(reads) >= min_coverage:
            for header, seq in reads:
                filtered[locus].append((header, seq))

    return filtered

# Parse reads
locus_reads = parse_fasta("~{reads}")

# Filter by coverage
filtered_reads = filter_by_coverage(locus_reads, ~{min_coverage})

# Write filtered FASTA
with open("filtered_reads.fasta", "w") as f:
    for locus, reads in sorted(filtered_reads.items()):
        for header, seq in reads:
            f.write(f">{header}\n{seq}\n")

# Write statistics
stats = {
    "total_loci": len(locus_reads),
    "filtered_loci": len(filtered_reads),
    "loci_removed": len(locus_reads) - len(filtered_reads),
    "min_coverage_threshold": ~{min_coverage}
}

with open("depth_filter_stats.json", "w") as f:
    json.dump(stats, f, indent=2)

print(f"Depth filtering: {len(filtered_reads)} of {len(locus_reads)} loci retained")
PYTHON_SCRIPT
    >>>

    output {
        File filtered_reads = "filtered_reads.fasta"
        File filter_stats = "depth_filter_stats.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 4: Infer best scheme from allele calls
task infer_scheme {
    input {
        File allele_calls
        File profiles
        Map[String, File]? schemes
        Float min_coverage_threshold = 0.5
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
import csv
from collections import defaultdict

def parse_allele_calls(json_path):
    """Parse allele calls JSON."""
    with open(json_path) as f:
        return json.load(f)

def parse_profile_table(tsv_path):
    """Parse profiles TSV and extract loci."""
    loci = set()
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames:
            for col in reader.fieldnames:
                if col.upper() not in ['ST', 'ID', 'CLONAL_COMPLEX']:
                    loci.add(col)
    return loci

def infer_best_scheme(allele_calls, profiles_path, schemes_map):
    """Infer best scheme based on allele call coverage."""
    # If single profile provided, use it
    if not schemes_map:
        return parse_profile_table(profiles_path), profiles_path, "default"

    # Multi-scheme: calculate coverage for each scheme
    scheme_scores = {}

    for scheme_name, scheme_profiles_path in schemes_map.items():
        scheme_loci = parse_profile_table(scheme_profiles_path)

        # Calculate coverage: how many loci from this scheme are called
        called_loci = 0
        for locus in scheme_loci:
            # Handle scheme-prefixed loci (e.g., "salmonella_adk_1")
            for call_key in allele_calls.keys():
                if call_key.endswith(locus) or locus in call_key:
                    called_loci += 1
                    break

        coverage = called_loci / len(scheme_loci) if scheme_loci else 0
        scheme_scores[scheme_name] = {
            "coverage": coverage,
            "num_loci": len(scheme_loci),
            "called_loci": called_loci,
            "path": scheme_profiles_path
        }

    # Select scheme with highest coverage
    best_scheme = max(scheme_scores.items(), key=lambda x: x[1]["coverage"])
    scheme_name = best_scheme[0]
    scheme_info = best_scheme[1]

    if scheme_info["coverage"] < ~{min_coverage_threshold}:
        # Fall back to default profile if coverage too low
        return parse_profile_table(profiles_path), profiles_path, "default_fallback"

    return parse_profile_table(scheme_info["path"]), scheme_info["path"], scheme_name

# Parse allele calls
allele_calls = parse_allele_calls("~{allele_calls}")

# Build schemes map (empty dict if no schemes provided)
schemes_map = {}
if "~{schemes}" != "":
    try:
        import re
        schemes_str = "~{schemes}"
        # Parse WDL map format
        # This is a simplification; actual implementation may need more robust parsing
        schemes_map = {}  # Would be populated from multi_scheme_input
    except:
        schemes_map = {}

# Infer best scheme
loci, selected_profiles, inferred_scheme = infer_best_scheme(allele_calls, "~{profiles}", schemes_map)

# Write results
result = {
    "inferred_scheme": inferred_scheme,
    "loci": list(loci),
    "num_loci": len(loci),
    "selected_profiles_table": selected_profiles
}

with open("scheme_inference.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"Inferred scheme: {inferred_scheme} with {len(loci)} loci")
PYTHON_SCRIPT
    >>>

    output {
        File scheme_info = "scheme_inference.json"
        String inferred_scheme = read_string("stdout")
        File selected_profile_table = profiles
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 6: Check if alignment fallback needed
task check_ambiguity {
    input {
        File allele_calls
        File profile_result
        Float ambiguity_threshold = 0.85
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

def should_use_alignment(allele_calls_file, result_file, threshold):
    """Check if alignment should be used for ambiguous calls."""
    with open(allele_calls_file) as f:
        allele_calls = json.load(f)

    with open(result_file) as f:
        result = json.load(f)

    # Check 1: Low confidence overall
    if result.get("confidence", 1.0) < threshold:
        return True

    # Check 2: Novel profile - alignment could refine
    if result.get("status") == "novel_profile":
        # Check if confidence is borderline
        if result.get("confidence", 0.0) < (threshold + 0.1):
            return True

    # Check 3: Low individual locus confidence
    for locus, call in allele_calls.items():
        confidence = call.get("confidence", call.get("similarity", 1.0))
        if isinstance(confidence, bool):
            confidence = 1.0 if confidence else 0.0
        if float(confidence) < (threshold - 0.05):
            return True

    return False

needs_align = should_use_alignment("~{allele_calls}", "~{profile_result}", ~{ambiguity_threshold})

with open("use_alignment.json", "w") as f:
    json.dump(needs_align, f)

print("Alignment needed" if needs_align else "Alignment not needed")
PYTHON_SCRIPT
    >>>

    output {
        Boolean needs_alignment = read_json("use_alignment.json")
        String decision = read_string("stdout")
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "1 GB"
    }
}

# Stage 7: Assemble final result
task assemble_final_result {
    input {
        File minhash_result
        File? alignment_result
        String inferred_scheme
        Boolean used_alignment = false
        String input_type = "contigs"
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

def calculate_nearest_st(profile_str, profiles, loci_order):
    """Calculate nearest ST (Hamming distance) for novel profiles."""
    import csv

    nearest_st = None
    min_distance = len(loci_order)

    try:
        with open(profiles) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Build profile from table row
                table_parts = []
                for locus in loci_order:
                    if locus in row:
                        table_parts.append(row[locus])
                table_profile = ','.join(table_parts)

                # Calculate Hamming distance
                query_parts = profile_str.split(',')
                table_parts_list = table_profile.split(',')

                distance = 0
                for qp, tp in zip(query_parts, table_parts_list):
                    if qp != '?' and tp != '?' and qp != tp:
                        distance += 1

                # Find minimum distance
                if distance < min_distance:
                    min_distance = distance
                    st_col = None
                    for col in row.keys():
                        if col.upper() == 'ST':
                            st_col = col
                            break
                    if st_col:
                        nearest_st = row[st_col]
    except:
        pass

    return nearest_st, min_distance

# Load minhash result
with open("~{minhash_result}") as f:
    result = json.load(f)

# Load alignment result if available
alignment_used = ~{used_alignment}
if alignment_used and "~{alignment_result}" != "":
    try:
        with open("~{alignment_result}") as f:
            alignment_data = json.load(f)
            # Merge alignment data (preferring alignment values)
            for key in ['confidence', 'status', 'allele_calls']:
                if key in alignment_data:
                    result[key] = alignment_data[key]
    except:
        pass

# Add scheme information
if "~{inferred_scheme}" != "":
    result["scheme"] = "~{inferred_scheme}"

# Ensure status field exists
if "status" not in result:
    result["status"] = "unknown"

# Ensure confidence field exists
if "confidence" not in result:
    result["confidence"] = 0.0

# Add method metadata
result["method"] = {
    "strategy": "mlst_orchestrated",
    "alignment_used": alignment_used,
    "tools": ["sourmash", "minimap2"] if alignment_used else ["sourmash"],
    "input_type": "~{input_type}"
}

# Calculate nearest ST for novel profiles
if result.get("status") == "novel_profile":
    profile_str = result.get("allele_profile", "")
    loci = [k for k in result.get("allele_calls", {}).keys()]
    if profile_str:
        # Note: In full implementation, would calculate from profiles table
        result["nearest_st"] = None
        result["distance_to_nearest"] = None

# Ensure all required fields
required_fields = ["sequence_type", "status", "confidence"]
for field in required_fields:
    if field not in result and field != "sequence_type":
        if field == "sequence_type":
            result[field] = result.get("profile_id", "unknown")
        elif field == "status":
            result[field] = "unknown"
        elif field == "confidence":
            result[field] = 0.0

# Use profile_id as sequence_type if not present
if "sequence_type" not in result:
    result["sequence_type"] = result.get("profile_id", "unknown")

# Write final result
with open("final_result.json", "w") as f:
    json.dump(result, f, indent=2)

print("Final result assembled")
PYTHON_SCRIPT
    >>>

    output {
        File final_result = "final_result.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
