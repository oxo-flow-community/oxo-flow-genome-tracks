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

Additional modes are ported (upstream v2.0.5 parity):
- **Single-cell mode** (`sc_enabled`, on by default): sc samples are declared
  with a BAM + 2-column barcode TSV (`sc_bam_dir` / `sc_metadata`); sinto
  `filterbarcodes` splits each sc BAM per cell-barcode group, the per-group
  BAMs are merged and covered exactly like bulk groups, and the sc groups join
  the same plots and hub.
- **IGV report** (`igv_report_enabled`, off by default — deactivated upstream
  at v2.0.5): `make_bed` + `igv_report` build a self-contained `igv-report.html`
  of all merged BAMs (bulk + sc) over the annotated gene regions with
  igv-reports `create_report`. Never in the default graph; opt in with
  `igv_report_enabled = true` and `-t igv_report`.
- **Conda env export** (`env_export_enabled`, off by default — upstream runs
  `env_export` in `rule all`, the port keeps it opt-in so the default graph is
  unchanged): `env_export_pygenometracks` / `env_export_sinto` /
  `env_export_igv_reports` run the upstream `conda env export` shell verbatim
  inside the env each rule declares, writing the resolved environment
  (versions + builds) to `results/genome_tracks/envs/*.yaml`. The checked-in
  `envs/*.yaml` serve the same reproducibility role, so these rules are
  documentation-only; they need the conda runtime (already a port
  requirement) and the envs built (guaranteed by the engine for the env each
  rule declares). DRAFT — mechanics live-verified with conda 26.1.1, the
  three real env builds not yet run.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
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
- **Single-cell data (optional)** — one aligned BAM per sc sample at
  `<sc_bam_dir>/<sc_id>.bam`, each read carrying a `CB` cell-barcode tag, plus
  one 2-column barcode TSV (`barcode<TAB>group`, no header) per sample at
  `<sc_metadata>/<sc_id>.tsv`. The unique group values of TSV column 2 become
  the sc groups and must be declared in `config.sc_groups` (and in
  `[[values]] sc_group`). `sc_enabled = false` removes the sc rules from the
  graph entirely.
- **Compute** — up to 4 CPUs / 4 GB per rule: `merge_bams`/`coverage`
  (samtools / bamCoverage) and their sc variants (`split_sc_bam` sinto,
  `merge_sc_bams`, `coverage_sc`) run with `threads = 4` / `4000M`; all other
  rules need 1 CPU / 1 GB — except `igv_report`, which is pinned to the
  upstream 8000 MB minimum (`igv_report_memory`, static — oxo-flow resources
  cannot express upstream's dynamic `max(2 * size_mb, 8000)`).
- **Tools** — conda environments with pinned versions, declared per rule in
  `main.oxoflow` (`[rules.environment]`): `envs/pygenometracks.yaml` pins
  samtools 1.19.2, deepTools 3.5.5, pyGenomeTracks 3.8, python 3.10.13,
  pip 24.0 and gtracks 1.12.6 (pip); `envs/sinto.yaml` pins sinto 0.10.0
  (sc split); `envs/igv_reports.yaml` pins igv-reports 1.14.1 / python 3.8 /
  pysam 0.22.0 (IGV report, opt-in only). Conda/mamba is required at runtime
  to build the environments; the helper rules (`annotate_genes`, `ucsc_hub`,
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
# 6. single-cell mode (sc_enabled, on by default): per-group sinto split +
#    merge + coverage of the sc groups
oxo-flow run main.oxoflow -t coverage_sc
# 7. IGV report (opt-in, deactivated by default like upstream)
oxo-flow run main.oxoflow -t igv_report igv_report_enabled=true
# 8. conda env export docs (opt-in, off by default — upstream runs env_export
#    in rule all; the checked-in envs/*.yaml serve the same role)
oxo-flow run main.oxoflow -t env_export_pygenometracks env_export_enabled=true
```

Results land under `results/genome_tracks/`: `merged_bams/` (bulk + sc),
`bigWigs/` (+ UCSC hub), `sc_bams/<sc_sample>/` (sinto splits), `tracks/{gene|region}.pdf`,
`igv-report.html` (opt-in), `configs/`, `genes_annotated.tsv`.

Configuration: all settings live in the `[config]` table of `main.oxoflow` —
`project_name`, `result_path`, `sample_annotation`, `bam_dir`, `gene_list`,
`genome_bed`, `bamCoverage_parameters` (default `-p max --binSize 10
--normalizeUsing RPGC --effectiveGenomeSize 2407883318` for mm10),
`track_colors` (comma-joined `group=#hex`, `#000000` default), `x_axis`,
`width`, `base_buffer`, `file_type`; single-cell keys `sc_enabled`,
`sc_bam_dir`, `sc_metadata`, `sc_groups`; IGV keys `igv_report_enabled`,
`igv_report_memory` (informational mirror of the rule's fixed 8000M);
env-export key `env_export_enabled` (opt-in, default off — the three
`env_export_*` rules write `results/genome_tracks/envs/*.yaml`).
Group fan-out is declared in `[[sample_groups]]` (bulk) and `[[values]]`
`sc_sample` × `sc_group` (single-cell), per-gene fan-out in `[[pairs]]`;
keep the annotation CSV, gene list, genome BED, metadata TSVs and BAMs in
sync with those tables and `config.bam_dir`. `config.samples_list` is the
single track-group namespace (upstream `plot_groups = sorted(sc_groups +
bulk_groups)`): declare the union of bulk and sc groups in it — the engine
re-consolidates it from the declared value plus `[[sample_groups]]` and
`[[pairs]]` names, so `[[pairs]]` ids also appear in it (pre-existing port
behavior; those groups have no BAM/bigWig and must be absent from
`track_colors`).

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
| `plot_tracks` | `plot_tracks` | gtracks 1.12.6, pyGenomeTracks 3.8 | identical `gtracks` invocation (coordinates, `--genes`, optional `--max ymax`, `--gene-rows`/`--genes-height` = isoform count, `--x-axis`, `--width`, `--color-palette` with `#000000` default); per-gene fan-out uses `[[pairs]]` `pair_id` (oxo-flow has no gene wildcard source); `depends_on = ["coverage"]` added because `expand_inputs` input lists do not form DAG edges in oxo-flow 0.12.0 |
| `ucsc_hub` | `ucsc_hub` | python3 (stdlib) | identical hub content (hub.txt, genomes.txt, trackDb.txt with hex→RGB colors, `../{group}.bw` relative symlinks) ported from the Python run block to `scripts/ucsc_hub.py`; the per-group symlinks are side effects (outputs declared only for the three text files) |
| `env_export` | `env_export_pygenometracks` / `env_export_sinto` / `env_export_igv_reports` | conda | upstream fans `{env}` over the three envs; oxo-flow cannot wildcard `[rules.environment]`, so one rule per env — identical `conda env export` shell, each rule exports its own activated env (the engine's `conda run -n <env>` wrapper + `conda env export` = the upstream semantics, verified live with conda 26.1.1). Upstream runs it in `rule all`; the port gates it on `env_export_enabled` (default off) so the default graph is unchanged — the checked-in `envs/*.yaml` serve the same reproducibility role. DRAFT (mechanics live-verified; the three real env builds not yet run) |
| `config_export` | `config_export` | python3 (stdlib) | `json.dump(config)` equivalent: `scripts/export_config.py` dumps the workflow's `[config]` table |
| `annot_export` | `annot_export` | cp | identical (`cp` of the annotation CSV) |
| `gene_list_export` | `gene_list_export` | cp | identical (`cp` of the gene list CSV) |
| `split_sc_bam` | `split_sc_bam` | sinto 0.10.0 | DRAFT (validated, not yet live-run): same `sinto filterbarcodes -b -c --outdir -p` command + upstream's touch-empty-bam fallback for groups absent in a sample; fan-out via `[[values]]` `sc_sample` × `sc_group` (upstream derives them from the metadata TSVs at load time; oxo-flow declares them — keep `[[values]] sc_sample`/`sc_group` in sync with `sc_bam_dir`/`sc_metadata`/`sc_groups`); upstream's `{sample}` = BAM-path md5 is replaced by readable sc ids |
| `merge_bams` (sc variant) | `merge_sc_bams` | samtools 1.19.2 | DRAFT: upstream switches `merge_bams` inputs per wildcard (sc groups read `sc_bams/`, bulk groups the annotation BAM column); oxo-flow cannot switch inputs per wildcard, so the sc variant is a separate rule writing the same `merged_bams/` namespace, gated on `sc_enabled` |
| `coverage` (sc variant) | `coverage_sc` | deepTools 3.5.5 | DRAFT: same `bamCoverage` command as bulk `coverage`; sc groups' bigWigs join `plot_tracks`/`ucsc_hub` via `config.samples_list` |
| `make_bed` | `make_bed` | awk | DRAFT: upstream projects `gene_annot_df` to `chr,start,end,name` in Python; the port uses an awk projection of `genes_annotated.tsv` (name,chr,start,end → BED4); gated on `igv_report_enabled` like the rule it feeds |
| `igv_report` | `igv_report` | igv-reports 1.14.1 | DRAFT: **temporarily deactivated upstream** (commented out of `rule all` at v2.0.5), ported as opt-in (`igv_report_enabled = true` + `-t igv_report`); same `create_report --genome --tracks --output` + the upstream `Variants`→`Genes and genomic regions` sed; track list = `config.samples_list` BAMs; memory fixed at the upstream 8000 MB minimum (oxo-flow resources are static) |
| Snakemake `report()` wrappers | — | — | no equivalent in oxo-flow; the report artifacts are written as plain files |

Configuration mapping: upstream `config/config.yaml` keys became `[config]`
keys with upstream defaults, except `result_path` (placeholder path →
`results`), `mem`/`threads` (→ per-rule `[rules.resources]`; upstream's
`4 × threads` for merge/coverage baked in as `threads = 4`), and
`track_colors` (YAML dict → comma-joined `group=#hex` string with the same
`#000000` default). New keys for the ported sc/IGV branches: `sc_enabled`,
`sc_bam_dir`, `sc_metadata` (directory keys — oxo-flow rule inputs cannot
index comma-joined config lists), `sc_groups` (merged + sorted into
`samples_list`), `igv_report_enabled`, `igv_report_memory` (documented mirror;
the rule's `memory` is the fixed upstream minimum), `env_export_enabled`
(opt-in, default off — upstream runs `env_export` in `rule all`; the checked-in
`envs/*.yaml` serve the same reproducibility role, so the port keeps the
default graph unchanged). Group fan-out uses
`[[sample_groups]]` (one `{sample}` per annotation group) + `[[values]]`
`sc_sample` × `sc_group`, gene fan-out uses `[[pairs]]`. Sample annotation,
gene list, genome BED, metadata TSVs and BAM files must be kept in sync with
those tables and `config.bam_dir`; the annotation CSV itself remains the
documentation record (`annot_export`).

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` (warnings acceptable, errors not) + `dry-run` +
a debug-instance check against `main.oxoflow`, including the single-cell
path (all `split_sc_bam`/`merge_sc_bams`/`coverage_sc` instances in the
default dry-run, `make_bed` + `igv_report` asserted to skip until
`igv_report_enabled = true` + `-t igv_report` brings them in), the IGV
shell rendering, and the `env_export_*` rules (skipped by default,
running with `env_export_enabled = true`). Expects `oxo-flow` on `PATH`
(or set `OXO=/path/to/oxo-flow`); CI runs the same script on every push.

The sc + IGV + env_export rules are DRAFT: `validate`/`lint`/`dry-run` and
the generated fixture BAMs (test/fixtures/sc_bams/*.bam, generated by
`test/fixtures/make_sc_fixtures.py` — stdlib-only, verified with pysam +
htslib index build) pass, but the sinto / igv-reports conda envs have not
been built and run live yet. Minimal live scope: on a Linux box with conda,
`conda env create -f envs/sinto.yaml` + `conda env create -f
envs/pygenometracks.yaml`, then `oxo-flow run main.oxoflow -t coverage_sc`
(fixtures already on disk; sinto needs a `CB` tag per read — the generated
BAMs have it); for the IGV report, additionally `conda env create -f
envs/igv_reports.yaml` and `oxo-flow run main.oxoflow -t igv_report
igv_report_enabled=true` (verify the `Variants`→`Genes and genomic regions`
label swap in the HTML); for env_export, `oxo-flow run main.oxoflow -t
env_export_pygenometracks env_export_enabled=true` (verify the resolved
YAML lands in `results/genome_tracks/envs/`).

## License

Apache-2.0 for the workflow (see [LICENSE](LICENSE)). Copyright (c) 2026
oxo-flow-community. The upstream workflow is MIT licensed; attribution and
the upstream license text are in [NOTICE.md](NOTICE.md) and
[LICENSE.upstream](LICENSE.upstream).
