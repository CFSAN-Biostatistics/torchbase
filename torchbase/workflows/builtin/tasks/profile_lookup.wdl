version 1.0

task lookup_profile {
    input {
        File allele_calls
        File profiles_table
        String strategy = "fast"
    }

    command <<<
        set -e
        python3 <<CODE
import json
import csv

def parse_profiles_table(tsv_path):
    """Parse TSV profiles table into dict of profiles."""
    profiles = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # First column is typically ST (sequence type) or profile_id
            profile_id = row[list(row.keys())[0]]
            # Rest of the columns are locus: allele mappings
            allele_profile = {}
            for locus, allele in row.items():
                if locus != list(row.keys())[0]:  # Skip the profile_id column
                    allele_profile[locus] = allele
            profiles[profile_id] = allele_profile
    return profiles

def parse_allele_calls(json_path):
    """Parse allele calls JSON."""
    with open(json_path) as f:
        return json.load(f)

def normalize_locus_name(locus):
    """Normalize locus name by removing scheme prefix if present."""
    # Handle scheme-prefixed locus names like "salmonella_adk"
    parts = locus.rsplit('_', 1)
    if len(parts) == 2:
        scheme_part, base_part = parts
        # Check if this looks like a scheme prefix (contains organism name)
        if any(org in scheme_part.lower() for org in ['salmonella', 'ecoli', 'listeria', 'campylobacter']):
            return locus  # Keep as-is if it has a recognizable scheme
    return locus

def find_matching_profile(query_alleles, profiles):
    """Find matching profile from allele calls.

    Returns: (profile_id, status, confidence)
    status can be "known", "novel_profile", or "novel_allele"
    """
    query_loci = set(query_alleles.keys())

    for profile_id, profile_alleles in profiles.items():
        profile_loci = set(profile_alleles.keys())

        # Check if this profile matches
        if query_loci == profile_loci:
            # Check if all alleles match
            all_match = True
            for locus, allele_id in query_alleles.items():
                expected_allele = str(profile_alleles.get(locus, ''))
                actual_allele = str(allele_id)
                if expected_allele != actual_allele:
                    all_match = False
                    break

            if all_match:
                return profile_id, "known", 0.95

    # No exact match found - novel profile
    return None, "novel_profile", 0.80

# Load data
profiles = parse_profiles_table("~{profiles_table}")
allele_calls = parse_allele_calls("~{allele_calls}")

# Find matching profile
profile_id, status, base_confidence = find_matching_profile(allele_calls, profiles)

# Build result JSON
result = {
    "profile_id": profile_id or "unknown",
    "st": profile_id or "unknown",
    "sequence_type": profile_id or "unknown",
    "status": status,
    "confidence": base_confidence,
    "allele_profile": ",".join([f"{locus}_{allele}" for locus, allele in sorted(allele_calls.items())]),
    "allele_calls": allele_calls,
    "method": {
        "strategy": "~{strategy}",
        "alignment_used": False,
        "tools": ["sourmash"]
    }
}

# Write result JSON
with open('typing_result.json', 'w') as f:
    json.dump(result, f, indent=2)

CODE
    >>>

    output {
        File typing_result = "typing_result.json"
    }
}
