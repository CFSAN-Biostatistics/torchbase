version 1.0

task lookup_profile {
    input {
        File allele_calls
        File profiles
        String strategy = "sensitive"
        Boolean alignment_used = true
    }

    command <<<
        set -e
        python3 <<'PYTHON_SCRIPT'
import json

def parse_tsv_profiles(profiles_path):
    """Parse profile table TSV."""
    profiles_dict = {}
    with open(profiles_path) as f:
        lines = f.readlines()

    if len(lines) < 2:
        return profiles_dict

    # First line is header
    header = lines[0].strip().split('\t')
    st_idx = 0
    locus_indices = {}

    for i, col in enumerate(header):
        if col.lower() in ['st', 'id', 'profile_id']:
            st_idx = i
        else:
            locus_indices[col] = i

    # Parse data rows
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) <= st_idx:
            continue

        st = parts[st_idx]
        profile = {}
        for locus, idx in locus_indices.items():
            if idx < len(parts):
                profile[locus] = parts[idx]

        profiles_dict[st] = profile

    return profiles_dict, locus_indices

def find_matching_profile(allele_calls, profiles_dict):
    """Find matching profile from allele calls."""
    for st, profile in profiles_dict.items():
        all_match = True
        for locus, expected_allele in profile.items():
            if locus not in allele_calls:
                all_match = False
                break

            called_allele = allele_calls[locus]
            if isinstance(called_allele, dict):
                called_allele = called_allele.get('allele_id', str(called_allele))

            if str(called_allele) != str(expected_allele):
                all_match = False
                break

        if all_match:
            return st

    return None

# Load inputs
with open('~{allele_calls}') as f:
    allele_calls_data = json.load(f)

profiles_dict, locus_indices = parse_tsv_profiles('~{profiles}')

# Extract allele IDs from calls
allele_ids = {}
for locus, call_info in allele_calls_data.items():
    if isinstance(call_info, dict):
        allele_ids[locus] = call_info.get('allele_id', '0')
    else:
        allele_ids[locus] = str(call_info)

# Find matching profile
matched_st = find_matching_profile(allele_ids, profiles_dict)

if matched_st is None:
    status = "novel_profile"
    confidence = 0.0
else:
    status = "known"
    # Calculate confidence from average identity
    identities = []
    for call_info in allele_calls_data.values():
        if isinstance(call_info, dict) and 'identity' in call_info:
            identities.append(call_info['identity'])

    confidence = sum(identities) / len(identities) if identities else 0.0

# Assemble allele profile string
profile_parts = []
for locus in sorted(locus_indices.keys()):
    if locus in allele_ids:
        profile_parts.append(f"{locus}_{allele_ids[locus]}")

allele_profile = ','.join(profile_parts)

# Construct output JSON
result = {
    "profile_id": matched_st or "novel",
    "status": status,
    "confidence": float(confidence),
    "allele_profile": allele_profile,
    "allele_calls": allele_calls_data,
    "method": {
        "strategy": "~{strategy}",
        "alignment_used": ~{alignment_used}
    },
    "notes": {
        "alignment_metrics": "Enhanced metrics from strict alignment parameters"
    }
}

# Write output
with open('typing_result.json', 'w') as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File typing_result = "typing_result.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
