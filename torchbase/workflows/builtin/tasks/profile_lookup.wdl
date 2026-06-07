version 1.0

task lookup_profile {
    input {
        File allele_calls
        File profiles_table
        String strategy = "balanced"
        Boolean alignment_used = false
        String scheme = ""
    }

    command <<<
        set -e
        python3 <<'PYTHON_SCRIPT'
import json
import csv
from collections import defaultdict

def parse_allele_calls(json_path):
    """Parse JSON allele calls from minhash or alignment."""
    with open(json_path) as f:
        return json.load(f)

def parse_profiles_table(tsv_path):
    """Parse profiles TSV file."""
    profiles = []
    loci_order = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            profiles.append(row)
        # Get loci order from header (skip first column which is ID)
        if reader.fieldnames:
            loci_order = [col for col in reader.fieldnames if col.upper() != 'ST']
    return profiles, loci_order

def build_profile_string(allele_calls, loci_order):
    """Build profile string from allele calls."""
    profile_parts = []
    for locus in loci_order:
        if locus in allele_calls:
            allele_id = allele_calls[locus].get('allele_id', '?')
            profile_parts.append(str(allele_id))
        else:
            profile_parts.append('?')
    return ','.join(profile_parts)

def lookup_profile(profile_str, profiles, loci_order):
    """Look up profile in profile table."""
    st_col = None
    for profile in profiles:
        if st_col is None:
            for col in profile.keys():
                if col.upper() == 'ST':
                    st_col = col
                    break

        table_parts = []
        for locus in loci_order:
            if locus in profile:
                table_parts.append(profile[locus])
        table_profile = ','.join(table_parts)

        match = True
        query_parts = profile_str.split(',')
        table_parts_list = table_profile.split(',')

        if len(query_parts) != len(table_parts_list):
            continue

        for query_part, table_part in zip(query_parts, table_parts_list):
            if query_part == '?':
                continue
            if table_part == '?':
                continue
            if query_part != table_part:
                match = False
                break

        if match and st_col:
            return profile[st_col], "known_profile"

    return None, "novel_profile"


def find_nearest_st(profile_str, profiles, loci_order):
    """Find nearest ST by Hamming distance for novel profiles."""
    query_parts = profile_str.split(',')
    best_st = None
    best_distance = float('inf')
    st_col = None

    for profile in profiles:
        if st_col is None:
            for col in profile.keys():
                if col.upper() == 'ST':
                    st_col = col
                    break
        if st_col is None:
            continue

        table_parts = [profile.get(locus, '?') for locus in loci_order]
        if len(table_parts) != len(query_parts):
            continue

        distance = sum(
            0 if q == '?' or t == '?' or q == t else 1
            for q, t in zip(query_parts, table_parts)
        )
        if distance < best_distance:
            best_distance = distance
            best_st = profile[st_col]

    return best_st, best_distance


# Main
allele_calls = parse_allele_calls("~{allele_calls}")
profiles, loci_order = parse_profiles_table("~{profiles_table}")

# Infer scheme from profiles_table path if not provided
import os
provided_scheme = "~{scheme}"
if provided_scheme:
    scheme_name = provided_scheme
else:
    scheme_name = os.path.basename(os.path.dirname("~{profiles_table}"))

# Build profile string
profile_str = build_profile_string(allele_calls, loci_order)

# Lookup profile
profile_id, status = lookup_profile(profile_str, profiles, loci_order)

# For novel profiles, find nearest ST
nearest_st = None
nearest_st_distance = None
if status == "novel_profile" and profiles:
    nearest_st, nearest_st_distance = find_nearest_st(profile_str, profiles, loci_order)

# Calculate confidence from allele calls
confidences = []
for locus, call in allele_calls.items():
    if 'confidence' in call:
        confidences.append(1.0 if call['confidence'] else 0.0)
    elif 'similarity' in call or 'identity' in call:
        score = call.get('similarity', call.get('identity', 0.0))
        confidences.append(float(score))

overall_confidence = sum(confidences) / len(confidences) if confidences else 0.0

# Determine if alignment was used
alignment_used = ~{alignment_used}

# Build result
result = {
    "profile_id": profile_id if profile_id else "unknown",
    "profile_type": "sequence_type",
    "scheme": scheme_name,
    "status": status,
    "confidence": max(0.0, min(1.0, overall_confidence)),
    "allele_profile": profile_str,
    "allele_calls": allele_calls,
    "method": {
        "strategy": "~{strategy}",
        "alignment_used": alignment_used,
        "tools": ["sourmash", "minimap2"] if alignment_used else ["sourmash"]
    },
    "notes": {
        "num_loci": len(loci_order),
        "num_called": len(allele_calls),
        "mean_confidence": overall_confidence
    }
}

if nearest_st is not None:
    result["nearest_st"] = nearest_st
    result["nearest_st_distance"] = nearest_st_distance

# Write result
with open('profile_result.json', 'w') as f:
    json.dump(result, f, indent=2)

# Write individual outputs for easier access
with open('status.txt', 'w') as f:
    f.write(status)
with open('profile_id.txt', 'w') as f:
    f.write(profile_id if profile_id else "unknown")

PYTHON_SCRIPT
    >>>

    output {
        File result = "profile_result.json"
        String status = read_string("status.txt")
        String profile_id = read_string("profile_id.txt")
    }

    runtime {
        cpu: 1
        memory: "2 GB"
    }
}
