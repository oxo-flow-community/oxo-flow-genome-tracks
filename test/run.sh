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

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1
grep -qE '\{(sample|pair_id|config\.)' /tmp/oxo-debug-$$.txt && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "==> debug: every ported rule instance present"
for inst in merge_bams_cohort_untreated merge_bams_cohort_treated \
            coverage_cohort_untreated coverage_cohort_treated \
            plot_tracks_Tmem26 plot_tracks_chr1-1000-2000 ucsc_hub \
            annotate_genes annot_export gene_list_export config_export; do
  grep -q "$inst" /tmp/oxo-debug-$$.txt || { echo "missing instance: $inst"; exit 1; }
done

echo "PASS"
