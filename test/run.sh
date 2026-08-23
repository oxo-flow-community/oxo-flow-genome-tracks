#!/usr/bin/env bash
# Acceptance test for oxo-flow-genome-tracks port (epigen/genome_tracks v2.0.5).
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> dry-run: single-cell path active (fixture default sc_enabled = true)"
for inst in split_sc_bam_sc_sample_sc1_sc_group_g1 split_sc_bam_sc_sample_sc1_sc_group_g2 \
            split_sc_bam_sc_sample_sc2_sc_group_g1 split_sc_bam_sc_sample_sc2_sc_group_g2 \
            merge_sc_bams_sc_group_g1 merge_sc_bams_sc_group_g2 \
            coverage_sc_sc_group_g1 coverage_sc_sc_group_g2; do
  grep -q "$inst" /tmp/oxo-dryrun-$$.txt || { echo "missing sc instance: $inst"; exit 1; }
done

echo "==> dry-run: IGV rules deactivated by default (upstream parity)"
grep -q "make_bed  \[skip: when condition false\]" /tmp/oxo-dryrun-$$.txt || { echo "make_bed not skipped by default"; exit 1; }
grep -q "igv_report  \[skip: when condition false\]" /tmp/oxo-dryrun-$$.txt || { echo "igv_report not skipped by default"; exit 1; }

echo "==> dry-run: IGV path opt-in (igv_report_enabled = true, -t igv_report)"
"$OXO" dry-run main.oxoflow -t igv_report igv_report_enabled=true > /tmp/oxo-igv-dryrun-$$.txt 2>&1
grep -q "make_bed  \[run:" /tmp/oxo-igv-dryrun-$$.txt || { echo "make_bed not in igv run"; exit 1; }
grep -q "igv_report  \[run:" /tmp/oxo-igv-dryrun-$$.txt || { echo "igv_report not in igv run"; exit 1; }

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1
grep -qE '\{(sample|pair_id|sc_sample|sc_group|config\.)' /tmp/oxo-debug-$$.txt && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "==> debug: every ported rule instance present"
for inst in merge_bams_cohort_untreated merge_bams_cohort_treated \
            coverage_cohort_untreated coverage_cohort_treated \
            plot_tracks_chr1-1000-2000 ucsc_hub \
            split_sc_bam_sc_sample_sc1_sc_group_g1 split_sc_bam_sc_sample_sc2_sc_group_g2 \
            merge_sc_bams_sc_group_g1 merge_sc_bams_sc_group_g2 \
            coverage_sc_sc_group_g1 coverage_sc_sc_group_g2 \
            annotate_genes annot_export gene_list_export config_export; do
  grep -q "$inst" /tmp/oxo-debug-$$.txt || { echo "missing instance: $inst"; exit 1; }
done

echo "==> debug: IGV shell renders the sc + bulk track list (igv_report_enabled = true)"
"$OXO" debug main.oxoflow -r igv_report > /tmp/oxo-igv-debug-$$.txt 2>&1
grep -q "create_report results/genome_tracks/genes.bed --genome mm10 --tracks" /tmp/oxo-igv-debug-$$.txt || { echo "igv_report shell not rendered"; exit 1; }

echo "PASS"
