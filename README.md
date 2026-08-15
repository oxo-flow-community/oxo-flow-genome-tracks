# oxo-flow-genome-tracks — Genome browser tracks: coverage, gene plots and UCSC hub

> ★ Verified · ⇄ Official port of [`epigen/genome_tracks`](https://github.com/epigen/genome_tracks) @ `v2.0.5` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-genome-tracks/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-genome-tracks/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Turn aligned BAM files (RNA-seq, ATAC-seq or any aligned/mapped data) into
publication-ready genome browser tracks. BAM files are merged per experimental
group with `samtools merge` and indexed; per-group coverage is computed with
deepTools `bamCoverage` into bigWig files (RPGC-normalized by default); genome
tracks per gene or genomic region — with isoform-aware layout, configurable
colors, axis and width — are plotted with the `gtracks` wrapper for
pyGenomeTracks; and a UCSC genome browser track hub is generated so all tracks
can be shared online. You get merged, indexed BAMs, bigWig coverage files,
per-gene/per-region track plots, and a UCSC track hub.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.11.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note: the
conda package may lag behind releases; other platform binaries are available
on the [releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-genome-tracks.git
cd oxo-flow-genome-tracks
```

### 3. Requirements

Derived from `main.oxoflow`:

- **Reference data (user-provided)** — point the `[config]` block of
  `main.oxoflow` at your files:
  - BAM files per group at `<bam_dir>/<group>/*.bam` (e.g.
    `test/fixtures/bams/untreated/*.bam`). Input BAMs need no index —
    `merge_bams` produces merged, indexed BAMs (`samtools merge` + `samtools index -b`).
  - A sample annotation CSV with a `group` column (`sample_annotation`, e.g.
    `test/fixtures/annotation.csv`); the group values are the fan-out unit for
    merging, coverage and the hub.
  - A gene list CSV with `gene_region,ymax` columns (`gene_list`, e.g.
    `test/fixtures/genes.csv`) — one gene symbol or `chr:start-end` region per row.
  - A 12-column genome BED for gene annotation (`genome_bed`, e.g.
    `test/fixtures/genome_bed/ref.bed.gz`). No genome FASTA or annotation GTF
    is required: genes are annotated from the BED, and bamCoverage normalizes
    via the built-in `--effectiveGenomeSize` (default config is mm10,
    `2407883318`).
- **Compute** — up to 4 CPUs / 4 GB per rule: `merge_bams` (samtools) and
  `coverage` (bamCoverage) run with `threads = 4` / `4000M`; all other rules
  need 1 CPU / 1 GB.
- **Tools** — conda environments with pinned versions, declared per rule in
  `main.oxoflow` (`[rules.environment]`): `envs/pygenometracks.yaml` pins
  samtools 1.19.2, deepTools 3.5.5, pyGenomeTracks 3.8, python 3.10.13,
  pip 24.0 and gtracks 1.12.6 (pip). Conda/mamba is required at runtime to
  build the environment; the helper rules (`annotate_genes`, `ucsc_hub`,
  `config_export`) need only a system `python3` (stdlib).
- **Disk** — results land under `results/genome_tracks/`; merged BAMs plus
  per-group bigWigs are the main outputs, modest in size for typical groups.

## Usage

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

Configuration: all settings live in the `[config]` table of `main.oxoflow` —
`project_name`, `result_path`, `sample_annotation`, `bam_dir`, `gene_list`,
`genome_bed`, `bamCoverage_parameters` (default `-p max --binSize 10
--normalizeUsing RPGC --effectiveGenomeSize 2407883318` for mm10),
`track_colors` (comma-joined `group=#hex`, `#000000` default), `x_axis`,
`width`, `base_buffer`, `file_type`. Group fan-out is declared in
`[[sample_groups]]`, per-gene fan-out in `[[pairs]]`; keep the annotation CSV,
gene list, genome BED and BAMs in sync with those tables and `config.bam_dir`.

## Source

Upstream: **[epigen/genome_tracks](https://github.com/epigen/genome_tracks)**
@ `v2.0.5` (commit `e7016b746d98be1e824ed6a3b767574458ce7cf7`), MIT license.
Created 2026-08-15; this workflow may lag behind upstream releases. See
[NOTICE.md](NOTICE.md) for attribution.

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

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` (warnings acceptable, errors not) + `dry-run` +
a debug-instance check against `main.oxoflow`. Expects `oxo-flow` on `PATH`
(or set `OXO=/path/to/oxo-flow`); CI runs the same script on every push.

## License

Apache-2.0 for the workflow (see [LICENSE](LICENSE)). Copyright (c) 2026
oxo-flow-community. The upstream workflow is MIT licensed; attribution and
the upstream license text are in [NOTICE.md](NOTICE.md) and
[LICENSE.upstream](LICENSE.upstream).
