version 1.0

task filter_alleles {
    input {
        File allele_fasta
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    command <<<
        set -e
        python3 <<'PYTHON_SCRIPT'
import json
from pathlib import Path

def parse_fasta(fasta_path):
    """Parse FASTA file into list of (header, sequence) tuples."""
    entries = []
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            entries.append((current_header, ''.join(current_seq)))
    return entries

def extract_locus_and_allele(header):
    """Extract locus name and allele ID from FASTA header.

    Expects format like: locus_name_allele_id
    """
    parts = header.split('_')
    if len(parts) >= 2:
        allele_id = parts[-1]
        locus = '_'.join(parts[:-1])
        return locus, allele_id
    return header, "unknown"

def load_quality_json(quality_path):
    """Load quality.json and extract suspect data."""
    if not quality_path or not Path(quality_path).exists():
        return {
            'suspect_alleles': set(),
            'suspect_loci': set(),
            'suspect_profiles': set()
        }

    with open(quality_path) as f:
        quality_data = f.read().strip()
        if not quality_data:
            return {
                'suspect_alleles': set(),
                'suspect_loci': set(),
                'suspect_profiles': set()
            }
        data = json.loads(quality_data)

    suspect_alleles = set()
    suspect_loci = set()

    # Extract suspect alleles and loci from quality data
    if 'loci' in data:
        for locus, locus_data in data['loci'].items():
            # Check if locus is flagged as suspect
            if locus_data.get('suspect', False):
                suspect_loci.add(locus)

            # Check for suspect alleles within this locus
            if 'alleles' in locus_data:
                for allele, allele_data in locus_data['alleles'].items():
                    if allele_data.get('suspect', False):
                        suspect_alleles.add(f"{locus}_{allele}")

            # Also check similarities for low-quality alleles
            if 'similarities' in locus_data:
                threshold = locus_data.get('threshold', 90.0)
                for pair, similarity in locus_data['similarities'].items():
                    if similarity < threshold:
                        # Mark alleles in low-similarity pairs as suspect
                        allele1, allele2 = pair.split('-')
                        suspect_alleles.add(f"{locus}_{allele1}")
                        suspect_alleles.add(f"{locus}_{allele2}")

    # Suspect profiles would mark entire loci
    suspect_profiles = set()
    if 'profiles' in data:
        for profile_id, profile_data in data['profiles'].items():
            if profile_data.get('suspect', False):
                # Get loci used in this profile
                if 'loci' in profile_data:
                    suspect_profiles.update(profile_data['loci'])

    return {
        'suspect_alleles': suspect_alleles,
        'suspect_loci': suspect_loci,
        'suspect_profiles': suspect_profiles
    }

# Parse inputs
allele_fasta = "~{allele_fasta}"
quality_json = "~{quality_json}" if "~{quality_json}" and "~{quality_json}" != "" else None
exclude_alleles = ~{true='True' false='False' exclude_suspect_alleles}
exclude_loci = ~{true='True' false='False' exclude_suspect_loci}
exclude_profiles = ~{true='True' false='False' exclude_suspect_profiles}

# Load alleles
alleles = parse_fasta(allele_fasta)

# Load quality data
suspect_data = load_quality_json(quality_json)

# Determine what to exclude
excluded_alleles = []
excluded_loci = set()
filtered_alleles = []

# Build exclusion sets based on flags
to_exclude_alleles = set()
to_exclude_loci = set()

if exclude_profiles:
    # Most aggressive: exclude all loci from suspect profiles
    to_exclude_loci.update(suspect_data['suspect_profiles'])
    to_exclude_loci.update(suspect_data['suspect_loci'])
    to_exclude_alleles.update(suspect_data['suspect_alleles'])
elif exclude_loci:
    # Medium: exclude suspect loci and suspect alleles
    to_exclude_loci.update(suspect_data['suspect_loci'])
    to_exclude_alleles.update(suspect_data['suspect_alleles'])
elif exclude_alleles:
    # Least aggressive: only exclude specific suspect alleles
    to_exclude_alleles.update(suspect_data['suspect_alleles'])

# Filter alleles
for header, sequence in alleles:
    locus, allele_id = extract_locus_and_allele(header)
    full_allele_name = f"{locus}_{allele_id}"

    # Check if this allele should be excluded
    exclude = False

    if locus in to_exclude_loci:
        exclude = True
        excluded_loci.add(locus)
    elif full_allele_name in to_exclude_alleles:
        exclude = True
        excluded_alleles.append(full_allele_name)

    if not exclude:
        filtered_alleles.append((header, sequence))

# Write filtered FASTA
with open('filtered_alleles.fasta', 'w') as f:
    for header, sequence in filtered_alleles:
        f.write(f'>{header}\n')
        f.write(f'{sequence}\n')

# Write exclusion metadata
exclusion_data = {
    'excluded_alleles': list(excluded_alleles),
    'excluded_loci': list(excluded_loci),
    'num_excluded_alleles': len(excluded_alleles),
    'num_excluded_loci': len(excluded_loci),
    'total_input_alleles': len(alleles),
    'total_output_alleles': len(filtered_alleles)
}

with open('exclusions.json', 'w') as f:
    json.dump(exclusion_data, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File filtered_fasta = "filtered_alleles.fasta"
        File exclusions = "exclusions.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
