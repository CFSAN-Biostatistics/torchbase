version 1.0

# Profile lookup task for allelic typing
# Matches allele calls against known profiles in a profile table
# Outputs status indicating profile type (known, novel_profile, or novel_allele)

task lookup_profile {
    input {
        File allele_calls
        File profiles_table
    }

    command <<<
        set -e
        python3 <<'PYTHON_SCRIPT'
import json
import csv
from collections import defaultdict

def parse_allele_calls(json_path):
    """Parse allele calls from JSON output of allele calling task."""
    try:
        with open(json_path) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}

def parse_profiles_table(tsv_path):
    """Parse profiles table (TSV format) into searchable structure."""
    profiles_by_id = {}
    profiles_list = []

    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames is None:
                return profiles_by_id, profiles_list

            for row in reader:
                profiles_list.append(row)
                # Index by profile ID if present
                if 'ST' in row or 'profile_id' in row or 'id' in row:
                    profile_id = row.get('ST') or row.get('profile_id') or row.get('id')
                    profiles_by_id[profile_id] = row
    except FileNotFoundError:
        pass

    return profiles_by_id, profiles_list

def normalize_allele_calls(allele_calls):
    """Convert allele calls to comparable format."""
    if isinstance(allele_calls, dict):
        # Extract allele IDs from dict format
        normalized = {}
        for locus, data in allele_calls.items():
            if isinstance(data, dict):
                normalized[locus] = data.get('allele_id', str(data))
            else:
                normalized[locus] = str(data)
        return normalized
    return {}

def compare_profiles(call_dict, profile_dict):
    """Compare allele calls against a profile row."""
    # Compare locus by locus
    mismatches = 0
    total_loci = 0

    for locus, call_allele in call_dict.items():
        # Look for matching locus in profile
        profile_allele = profile_dict.get(locus)
        if profile_allele is None:
            # Try alternative naming: scheme_locus format
            for key in profile_dict.keys():
                if key.endswith('_' + locus) or key == locus:
                    profile_allele = profile_dict[key]
                    break

        if profile_allele is not None:
            total_loci += 1
            # Handle special values
            if profile_allele == '?':  # Wildcard in profile
                continue
            elif profile_allele == 'X':  # Exclusion marker
                mismatches += 1
            elif str(call_allele) != str(profile_allele):
                mismatches += 1

    if total_loci == 0:
        return False, 0.0

    similarity = (total_loci - mismatches) / total_loci
    return mismatches == 0, similarity

def detect_novel_allele(allele_calls, profiles_table_data):
    """Check if any allele in the call set is novel (not in any profile)."""
    called_alleles = set()
    for locus, data in allele_calls.items():
        if isinstance(data, dict):
            allele_id = data.get('allele_id')
        else:
            allele_id = data
        called_alleles.add((locus, str(allele_id)))

    known_alleles = set()
    for profile_row in profiles_table_data:
        for locus, allele_id in profile_row.items():
            if allele_id and allele_id not in ['?', 'X', 'ST', 'profile_id', 'id']:
                known_alleles.add((locus, str(allele_id)))

    # Check if any called allele is not in known set
    for locus, allele_id in called_alleles:
        if (locus, allele_id) not in known_alleles:
            return True

    return False

# Main execution
allele_calls = parse_allele_calls("~{allele_calls}")
profiles_by_id, profiles_list = parse_profiles_table("~{profiles_table}")

normalized_calls = normalize_allele_calls(allele_calls)

# Try to find exact matching profile
best_match_id = None
best_match_similarity = 0.0
best_match_exact = False

for profile_id, profile_row in profiles_by_id.items():
    is_exact, similarity = compare_profiles(normalized_calls, profile_row)
    if is_exact:
        best_match_id = profile_id
        best_match_exact = True
        best_match_similarity = 1.0
        break  # Found exact match
    elif similarity > best_match_similarity:
        best_match_similarity = similarity
        best_match_id = profile_id

# If no exact match found in indexed profiles, search full profile list
if not best_match_exact and profiles_list:
    for profile_row in profiles_list:
        is_exact, similarity = compare_profiles(normalized_calls, profile_row)
        if is_exact:
            profile_id = profile_row.get('ST') or profile_row.get('profile_id') or profile_row.get('id')
            best_match_id = profile_id
            best_match_exact = True
            best_match_similarity = 1.0
            break
        elif similarity > best_match_similarity:
            profile_id = profile_row.get('ST') or profile_row.get('profile_id') or profile_row.get('id')
            best_match_similarity = similarity
            best_match_id = profile_id

# Determine status
if best_match_exact:
    status = "known"
elif detect_novel_allele(allele_calls, profiles_list):
    status = "novel_allele"
else:
    status = "novel_profile"

result = {
    "status": status,
    "profile_id": best_match_id or "unknown",
    "similarity": best_match_similarity
}

# Write output JSON
with open('lookup_result.json', 'w') as f:
    json.dump(result, f, indent=2)

# Write status string for WDL output
with open('status.txt', 'w') as f:
    f.write(status)

PYTHON_SCRIPT
    >>>

    output {
        File lookup_result = "lookup_result.json"
        String status = read_string("status.txt")
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
