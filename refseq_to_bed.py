#!/usr/bin/env python3
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        prog="refseq_to_bed",
        description="Convert RefSeq GTF or BED into region BED for bam_statistics_v4.py"
    )
    p.add_argument("--gtf", help="RefSeq GTF file", default=None)
    p.add_argument("--bed", help="RefSeq BED file (already exon/CDS)", default=None)
    p.add_argument("--genes", help="Comma-separated gene symbols to keep", default=None)
    p.add_argument("--transcripts", help="Comma-separated transcript IDs (e.g. NM_...) to keep", default=None)
    p.add_argument("--coding-only", action="store_true", help="Keep only CDS features")
    p.add_argument("--output", "-o", required=True, help="Output BED")
    return p.parse_args()

def load_filters(args):
    genes = set()
    transcripts = set()
    if args.genes:
        genes = {g.strip() for g in args.genes.split(",") if g.strip()}
    if args.transcripts:
        transcripts = {t.strip() for t in args.transcripts.split(",") if t.strip()}
    return genes, transcripts

def from_gtf(args, genes_filter, tx_filter):
    out = open(args.output, "w")
    with open(args.gtf) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr = cols
            if feature not in ("exon", "CDS"):
                continue
            if args.coding_only and feature != "CDS":
                continue
            # parse attributes
            attrs = {}
            for part in attr.split(";"):
                part = part.strip()
                if not part:
                    continue
                if " " in part:
                    k, v = part.split(" ", 1)
                    attrs[k] = v.strip().strip('"')
            gene_name = attrs.get("gene_name") or attrs.get("gene") or ""
            transcript_id = attrs.get("transcript_id") or ""
            # filters
            if genes_filter and gene_name not in genes_filter:
                continue
            if tx_filter and transcript_id not in tx_filter:
                continue
            info_fields = [
                "type=REFSEQ",
                f"gene={gene_name}" if gene_name else "gene=",
                f"transcript={transcript_id}" if transcript_id else "transcript=",
                f"feature={feature}"
            ]
            out.write("\t".join([chrom, start, end, ";".join(info_fields)]) + "\n")
    out.close()

def from_bed(args, genes_filter, tx_filter):
    out = open(args.output, "w")
    with open(args.bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[:3]
            info = cols[3:] if len(cols) > 3 else []
            info_str = "\t".join(info)
            # simple pass-through, optional filtering via gene= / transcript= tags in extra cols if present
            out.write("\t".join([chrom, start, end, info_str]) + "\n")
    out.close()

def main():
    args = parse_args()
    if not args.gtf and not args.bed:
        sys.stderr.write("Error: provide --gtf or --bed\n")
        sys.exit(1)
    genes_filter, tx_filter = load_filters(args)
    if args.gtf:
        from_gtf(args, genes_filter, tx_filter)
    else:
        from_bed(args, genes_filter, tx_filter)

if __name__ == "__main__":
    main()
