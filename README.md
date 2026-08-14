# Genome Browser Track Visualization Workflow (oxo-flow port)

[![CI](https://github.com/oxo-flow-community/oxo-flow-genome-tracks/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-genome-tracks/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

A port of the [epigen/genome_tracks](https://github.com/epigen/genome_tracks)
Snakemake workflow to oxo-flow: generation of genome browser tracks of
aligned/mapped BAM files (e.g., RNA-seq, ATAC-seq). BAM files are merged per
group with `samtools merge`, coverage is computed with deepTools
`bamCoverage` into bigWig files, and genome tracks per gene / genomic region
are plotted with the `gtracks` wrapper for pyGenomeTracks. A UCSC genome
browser track hub is created for online sharing.

## Source

Ported from **[epigen/genome_tracks](https://github.com/epigen/genome_tracks)**,
version `v2.0.5` (commit `e7016b746d98be1e824ed6a3b767574458ce7cf7`, MIT
license). This port is maintained independently and **may lag the upstream** —
check the `v2.0.5` tag above and the fidelity table below for the exact
ported state. Port date: 2026-08-15.

## Fidelity

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `merge_bams` | `merge_bams` | samtools 1.19.2 | identical command (`samtools merge -@ N` + `samtools index -@ N -b`); the per-group BAM list comes from `{config.bam_dir}/{group}/*.bam` glob instead of the annotation CSV's `bam` column (which is still copied verbatim by `annot_export`); `threads: 4 × config.threads` baked in as `threads = 4` |
| `coverage` | `coverage` | deepTools 3.5.5 | identical command incl. `-p max --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2407883318` default and `> {bw}.log 2>&1` redirect |
| Snakefile load-time gene annotation (`parse_gene`/`parse_region`, `gene_annot_df`) | `annotate_genes` | python3 (stdlib) | new single-instance rule; same algorithm (BED scan, min start / max end across isoforms, `base_buffer` extension for genes, no buffer for `chr:start-end` regions, `genes_not_found.csv`, `:`→`-` name replacement); upstream computes it in the Snakemake base env (numpy/pandas) — the port script uses only stdlib, so the upstream `global.yaml` env is not needed |
| `plot_tracks` | `plot_tracks` | gtracks 1.12.6, pyGenomeTracks 3.8 | identical `gtracks` invocation (coordinates, `--genes`, optional `--max ymax`, `--gene-rows`/`--genes-height` = isoform count, `--x-axis`, `--width`, `--color-palette` with `#000000` default); per-gene fan-out uses `[[pairs]]` `pair_id` (oxo-flow has no gene wildcard source); `depends_on = ["coverage"]` added because `expand_inputs` input lists do not form DAG edges in oxo-flow 0.11.0 |
| `ucsc_hub` | `ucsc_hub` | python3 (stdlib) | identical hub content (hub.txt, genomes.txt, trackDb.txt with hex→RGB colors, `../{group}.bw` relative symlinks) ported from the Python run block to `scripts/ucsc_hub.py`; the per-group symlinks are side effects (outputs declared only for the three text files) |
| `env_export` | not ported | — | upstream requests `conda env export` for the pygenometracks/igv_reports/sinto envs; needs a conda runtime and documents envs of branches not ported — the checked-in `envs/pygenometracks.yaml` serves the same reproducibility role |
| `config_export` | `config_export` | python3 (stdlib) | `json.dump(config)` equivalent: `scripts/export_config.py` dumps the workflow's `[config]` table |
| `annot_export` | `annot_export` | cp | identical (`cp` of the annotation CSV) |
| `gene_list_export` | `gene_list_export` | cp | identical (`cp` of the gene list CSV) |
| `make_bed` | not ported | — | only feeds the deactivated `igv_report` rule (not in the default target) |
| `split_sc_bam` | not ported | — | single-cell branch (sinto 0.10.0); not on the default path (no `.tsv` `group` entries in the default annotation) |
| `igv_report` | not ported | — | **temporarily deactivated upstream** (commented out of `rule all` at v2.0.5) |
| Snakemake `report()` wrappers | — | — | no equivalent in oxo-flow; the report artifacts are written as plain files |

Configuration mapping: upstream `config/config.yaml` keys became `[config]`
keys with upstream defaults, except `result_path` (placeholder path →
`results`), `mem`/`threads` (→ per-rule `[rules.resources]`; upstream's
`4 × threads` for merge/coverage baked in as `threads = 4`), and
`track_colors` (YAML dict → comma-joined `group=#hex` string with the same
`#000000` default). Group fan-out uses `[[sample_groups]]` (one `{sample}`
per annotation group), gene fan-out uses `[[pairs]]`. Sample annotation,
gene list, genome BED and BAM files must be kept in sync with `[[sample_groups]]`/
`[[pairs]]` and `config.bam_dir`; the annotation CSV itself remains the
documentation record (`annot_export`).

## Quickstart

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data: BAMs under <bam_dir>/<group>/*.bam + annotation.csv +
#    genes.csv + genome BED (see test/fixtures/ and config/README.md upstream)
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 8
# 5. run a subset (one group / one gene)
oxo-flow run main.oxoflow --samples first:1
oxo-flow run main.oxoflow -t plot_tracks
```

Results land under `results/genome_tracks/`: `merged_bams/`, `bigWigs/`
(+ UCSC hub), `tracks/{gene|region}.pdf`, `configs/`, `genes_annotated.tsv`.

## Requirements

- **oxo-flow ≥ 0.11.0** — install the prebuilt binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

- Conda users may alternatively `conda install -c bioconda oxo-flow-cli`
  (note: the bioconda package currently lags the release binary at 0.10.2 —
  some 0.11.0 format features may not validate).
- Docker/Singularity/conda at runtime, per the environments declared in
  `main.oxoflow` (`envs/pygenometracks.yaml` provides samtools, deepTools,
  pyGenomeTracks and gtracks with exact pins; the helper rules need only a
  system `python3`).

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md).

## Community

https://oxo-flow-community.github.io/
