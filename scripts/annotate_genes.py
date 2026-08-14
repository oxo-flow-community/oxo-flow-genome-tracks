#!/usr/bin/env python3
"""Port of the upstream genome_tracks Snakefile load-time gene annotation logic.

Upstream (workflow/Snakefile + workflow/rules/common.smk):
  * parse_gene():   scan the gzipped 12-column genome BED for every line whose
                    name (col 4) equals the gene, take the min start and max
                    end across isoforms, extend by base_buffer, count isoforms.
  * parse_region(): for "chr:start-end" entries take the coordinates verbatim,
                    count = 1 (no BED lookup, no base_buffer).
  * genes not found in the BED are dropped and written to
                    {result_path}/genes_not_found.csv.
  * gene names with ":" are replaced by "-" (output file naming).

Output: tab-separated table with one row per gene:
        name<TAB>chr<TAB>start<TAB>end<TAB>count<TAB>ymax
The ymax column comes from the gene_list CSV (config key 'gene_list').

Usage:
  annotate_genes.py --gene-list FILE.csv --genome-bed FILE.bed.gz \
      --base-buffer N --out OUT.tsv --not-found-out OUT.csv
"""
import argparse
import csv
import gzip
import re
import sys

REGION_RE = re.compile(r"^chr[0-9XY]+:[0-9]+-[0-9]+$")


def parse_gene(gene, bed_path, base_buffer):
    """Return (chr, start, end, count) of the gene or None if not found."""
    count = 0
    start = end = chrom = None
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            parsed_line = line.split()
            if parsed_line[3] == gene:
                count += 1
                tmp_chrom, tmp_start, tmp_end = parsed_line[:3]
                if count == 1:
                    chrom, start, end = tmp_chrom, tmp_start, tmp_end
                else:
                    if int(tmp_start) < int(start):
                        start = tmp_start
                    if int(tmp_end) > int(end):
                        end = tmp_end
        if count == 0:
            return None
    return chrom, int(start) - base_buffer, int(end) + base_buffer, count


def parse_region(region):
    chrom, start, end = region.replace("-", ":").split(":")
    return chrom, int(start), int(end), 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-list", required=True)
    ap.add_argument("--genome-bed", required=True)
    ap.add_argument("--base-buffer", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--not-found-out", required=True)
    args = ap.parse_args()

    gene_dict = {}
    with open(args.gene_list, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_dict[row["gene_region"]] = row["ymax"]

    genes = list(gene_dict.keys())

    # find gtrack parameters for each gene and handle genes not found
    gene_annot_list = []
    remove_genes = []
    for gene in genes:
        if REGION_RE.match(gene):
            gene_annot_list.append(parse_region(gene) + (gene_dict[gene],))
        else:
            tmp_val = parse_gene(gene, args.genome_bed, args.base_buffer)
            if tmp_val is None:
                remove_genes.append(gene)
            else:
                gene_annot_list.append(tmp_val + (gene_dict[gene],))

    # save text file containing genes that were not found in the genome BED
    if len(remove_genes) > 0:
        with open(args.not_found_out, "w", newline="") as f:
            for g in remove_genes:
                f.write(g + "\n")
        genes = [g for g in genes if g not in remove_genes]

    # ':' -> '-' in gene names (output file naming, as upstream)
    genes = [g.replace(":", "-") for g in genes]

    with open(args.out, "w") as f:
        for name, (chrom, start, end, count, ymax) in zip(
            genes, gene_annot_list
        ):
            f.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    name, chrom, start, end, count, ymax
                )
            )


if __name__ == "__main__":
    sys.exit(main())
