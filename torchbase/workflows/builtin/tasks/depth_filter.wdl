version 1.0

# Filter low-coverage sequences before allele calling.
# For contig inputs (input_type="contigs") this is a pass-through.
# For read inputs it removes sequences with mean k-mer depth below min_depth.

task depth_filter {
    input {
        File sequences
        String input_type = "contigs"
        Int min_depth = 3
        Int ksize = 21
    }

    command <<<
        set -e
        python3 <<'PYTHON_SCRIPT'
import sys

input_type = "~{input_type}"
min_depth = int("~{min_depth}")
ksize = int("~{ksize}")

def parse_fasta(path):
    records = []
    header, seq_lines = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header, seq_lines = line, []
            else:
                seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records

def kmer_depth(seq, k):
    """Estimate mean k-mer depth via simple frequency table."""
    if len(seq) < k:
        return 0.0
    counts = {}
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        counts[km] = counts.get(km, 0) + 1
    return sum(counts.values()) / len(counts)

records = parse_fasta("~{sequences}")
kept = []
removed = 0

# Contigs: no depth filtering — pass through unchanged
if input_type == "contigs":
    kept = records
else:
    for header, seq in records:
        depth = kmer_depth(seq, ksize)
        if depth >= min_depth:
            kept.append((header, seq))
        else:
            removed += 1

with open("filtered_sequences.fasta", "w") as out:
    for header, seq in kept:
        out.write(f"{header}\n{seq}\n")

import json
with open("depth_filter_stats.json", "w") as f:
    json.dump({
        "input_sequences": len(records),
        "kept_sequences": len(kept),
        "removed_sequences": removed,
        "input_type": input_type,
        "min_depth": min_depth,
        "ksize": ksize,
    }, f, indent=2)

sys.stderr.write(
    f"depth_filter: kept {len(kept)}/{len(records)} sequences "
    f"(removed {removed} below depth {min_depth})\n"
)
PYTHON_SCRIPT
    >>>

    output {
        File filtered_sequences = "filtered_sequences.fasta"
        File stats = "depth_filter_stats.json"
    }

    runtime {
        cpu: 1
        memory: "2 GB"
    }
}
