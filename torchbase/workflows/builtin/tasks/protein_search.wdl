version 1.0

# Protein-space search: tblastn of the operon's subunit reference set against
# translated nucleotide contigs (docs/operon-strategy-plan.md §4, mirrors
# stxtyper.cpp:887 including -gapextend 2). Normalizes BLAST tabular output
# to an internal HSP JSON list at this task boundary — nothing downstream
# ever sees BLAST format, so swapping the search engine later (e.g. for a
# Phraya protein-mode aligner once it gains affine gaps) only touches this
# task.
task search_subunits {
    input {
        File contigs
        File subunit_reference
        File operon_config_json
        Int gapextend = 2
        Float evalue = 1e-10
    }

    command <<<
        set -e
        makeblastdb -in ~{contigs} -dbtype nucl -out contigs_db

        tblastn \
            -query ~{subunit_reference} \
            -db contigs_db \
            -gapextend ~{gapextend} \
            -evalue ~{evalue} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen sstrand qseq sseq" \
            -out hits.tsv

        python3 <<'PYTHON_SCRIPT'
import json

HEADER_FORMAT_DEFAULT = "accession|subunit_role|reference_subtype|class"
with open("~{operon_config_json}") as f:
    _operon_cfg = json.load(f)
HEADER_FORMAT = _operon_cfg.get("reference", {}).get(
    "header_format", HEADER_FORMAT_DEFAULT
).split("|")


def parse_ref_header(qseqid):
    """Split a reference header (default:
    accession|subunit_role|reference_subtype|class) into a dict keyed by
    HEADER_FORMAT field names. `subunit_role` is the structural position
    ("A"/"B", matched against operon.subunit_order); `class` is the
    already-collapsed operon-level class token used for threshold lookup
    and cross-subunit class agreement (e.g. both stx2a and stx2c subunit-A
    references carry class "2")."""
    parts = qseqid.split("|")
    return dict(zip(HEADER_FORMAT, parts))


hsps = []
try:
    with open("hits.tsv") as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]
except FileNotFoundError:
    lines = []

for line in lines:
    fields = line.split("\t")
    (
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, nident,
        qlen, slen, sstrand, qseq, sseq,
    ) = fields

    ref = parse_ref_header(qseqid)
    subunit_role = ref.get("subunit_role", qseqid)
    reference_subtype = ref.get("reference_subtype", qseqid)
    class_label = ref.get("class", reference_subtype)

    sstart_i, send_i = int(sstart), int(send)
    start, stop = min(sstart_i, send_i), max(sstart_i, send_i)
    slen_i = int(slen)

    hsps.append({
        "subunit": subunit_role,
        "contig": sseqid,
        "strand": "+" if sstrand == "plus" else "-",
        "start": start,
        "stop": stop,
        "reference_accession": ref.get("accession", qseqid),
        "ref_subtype": reference_subtype,
        "ref_class": class_label,
        "nident": int(nident),
        "aln_len": int(length),
        "ref_len": int(qlen),
        # 0-based reference (query) offset where the aligned block starts,
        # plus the gapped query/subject alignment strings — together these
        # let operon_call.wdl project a fixed reference coordinate (e.g.
        # A312) back to the observed subject residue (§3.1, §7 risk 1-2).
        "qstart": int(qstart) - 1,
        "qseq": qseq,
        "sseq": sseq,
        # Internal-stop detection: tblastn reports a translated subject
        # alignment (`sseq`); an in-frame stop shows up as a literal "*".
        "internal_stop": "*" in sseq,
        # Contig-end proximity, small buffer for BLAST edge trimming.
        "contig_end": start <= 2 or stop >= slen_i - 2,
        # Frameshift detection requires stitching co-linear HSPs split by a
        # frame change (§7 risk 1-2) — left for operon_assembly.wdl, which
        # has visibility across all HSPs for a subunit/contig pair.
        "frameshift": False,
    })

with open("hsps.json", "w") as f:
    json.dump(hsps, f, indent=2)
PYTHON_SCRIPT
    >>>

    output {
        File hsps = "hsps.json"
    }

    runtime {
        docker: "ncbi/blast:2.16.0"
        cpu: 2
        memory: "4 GB"
    }
}
