# NuRepress: inferring transcriptional repressors from phased nucleosome architecture

NuRepress is an end-to-end workflow for identifying well-phased nucleosome arrays in repressive chromatin, stratifying array subtypes by accessibility patterns, performing motif enrichment on subtype-defined regions, and prioritizing candidate transcriptional repressors by integrating chromatin structure, sequence evidence, and transcriptional effects.

## Workflow

![NuRepress workflow](docs/figures/workflow_overview.png)

## Software dependencies

NuRepress is intended primarily for Linux or Linux-like HPC environments.

### External tools
Required for the full workflow:
- `bash`
- `python`
- `Rscript`
- `samtools`
- `danpos.py` from DANPOS3
- `sicer` from SICER2, when `--peak_mode sicer2` is used
- `findMotifsGenome.pl` from HOMER, when motif enrichment is run

Conditionally required:
- `bedtools`, especially when SICER2 processes BAM input
- HOMER genome packages

### Python packages
The bundled Python scripts use:
- `pandas`
- `numpy`
- `tqdm`
- `pyBigWig`

### R packages
CRAN packages:
- `readr`
- `dplyr`
- `tidyr`
- `stringr`
- `ggplot2`
- `data.table`
- `purrr`
- `cluster`
- `mclust`

Bioconductor packages:
- `GenomicRanges`
- `GenomicFeatures`
- `IRanges`
- `S4Vectors`
- `rtracklayer`
- `Rsamtools`
- `Biostrings`

---

## Recommended environment setup

### Option 1. Conda / mamba

```bash
mamba env create -f environment.yml
mamba activate nurepress
cd /path/to/NuRepress
```

### Option 2. Python packages only

```bash
pip install -r requirements.txt
```

### Option 3. Install R packages inside an existing R environment

```bash
Rscript install_r_packages.R
```

### PATH configuration

```bash
sh setup_permissions.sh
export PATH="/path/to/NuRepress/scripts:$PATH"
```


---

## Overview

NuRepress is driven directly from the integrated CLI script and **does not require any sample sheet**. All sample-specific inputs are provided through command-line arguments such as `--samples`, `--treatment_map`, `--atac_bw_map`, and related map-style options.

At a high level, NuRepress takes:
- repressive chromatin BAM files for nucleosome calling and array identification
- accessibility signal tracks for subtype clustering
- genome reference and annotation files for downstream annotation and motif analysis
- an expression matrix for final candidate scoring

It produces:
- nucleosome positioning outputs
- well-phased array calls
- accessibility-based subtype assignments
- TSS/promoter annotation results
- motif enrichment results
- ranked candidate transcriptional repressors

---

## Quick start

### Full workflow with control BAMs

```bash
run_nurepress_cli.py run \
  --out_dir /path/to/output \
  --steps nucleosome,array,cluster,describe,motif,score \
  --samples sample1,sample2,sample3 \
  --atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw;sample3=/path/to/sample3.atac.bw" \
  --treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam" \
  --control_map "sample1=/path/to/sample1.input.bam;sample2=/path/to/sample2.input.bam;sample3=/path/to/sample3.input.bam" \
  --genome hg38 \
  --genome_fasta /path/to/genome.fa \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --annotation_gtf /path/to/annotation.gtf \
  --expr_tsv /path/to/expression.tsv
```

### Full workflow without control BAMs

```bash
run_nurepress_cli.py run \
  --out_dir /path/to/output \
  --steps nucleosome,array,cluster,describe,motif,score \
  --samples sample1,sample2,sample3 \
  --atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw;sample3=/path/to/sample3.atac.bw" \
  --treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam" \
  --genome hg38 \
  --genome_fasta /path/to/genome.fa \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --annotation_gtf /path/to/annotation.gtf \
  --expr_tsv /path/to/expression.tsv
```

### Rerun only motif and score using existing outputs

```bash
run_nurepress_cli.py run \
  --out_dir /path/to/output \
  --steps motif,score \
  --samples sample1,sample2,sample3 \
  --annotation_gtf /path/to/annotation.gtf \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --genome hg38 \
  --expr_tsv /path/to/expression.tsv
```

### Rerun only score using existing describe and motif outputs

```bash
run_nurepress_cli.py score \
  --out_dir /path/to/output \
  --expr_tsv /path/to/expression.tsv
```

---


## Repository contents

### Top-level integrated CLI
- `run_nurepress_cli.py`  
  Main entry point for the integrated workflow controller and stepwise subcommands.

### Bundled step scripts
- `00_nucleosome/call_nucleosomes.sh`
- `01_array/call_peak_and_identify_array.py`
- `02_cluster/cluster_array_subtype.R`
- `02_cluster/describe_array_subtype.R`
- `03_motif/run_array_motif_enrichment.py`
- `04_score/score_repressor_candidates.R`

### Environment / installation files
- `environment.yml`
- `requirements.txt`
- `install_r_packages.R`
- `setup_permissions.sh`

---


## Command structure

NuRepress uses the following subcommands:

```bash
run_nurepress_cli.py run
run_nurepress_cli.py nucleosome
run_nurepress_cli.py array
run_nurepress_cli.py cluster
run_nurepress_cli.py describe
run_nurepress_cli.py motif
run_nurepress_cli.py score
```

### `run`
Top-level controller for running multiple steps in sequence.

### `nucleosome`
Run only the DANPOS wrapper step.

### `array`
Run only the broad-peak / well-phased array calling step.

### `cluster`
Prepare clustering inputs and run accessibility-based subtype clustering.

### `describe`
Run TSS/promoter annotation and subtype description.

### `motif`
Run motif target/background construction and motif enrichment.

### `score`
Run candidate repressor scoring.

---

## Input format overview

### 1. Sample list

All workflows start from a comma-separated sample list:

```bash
--samples sample1,sample2,sample3
```

The sample names defined here are the keys used everywhere else in the mapping arguments.

---

### 2. Map-style sample arguments

Several arguments use the same format:

```bash
sample1=/path/to/file1;sample2=/path/to/file2;sample3=/path/to/file3
```

For example:

```bash
--treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam"
```

NuRepress parses the following mapping arguments:

- `--atac_bw_map`
- `--positions_map`
- `--treatment_map`
- `--control_map`
- `--peak_bed_map`
- `--dpos_map`
- `--insertion_bw_map`

### Sample naming rule
The sample names in every map must exactly match the names in `--samples`.

For example, if `--samples sample1,sample2` is used, then a valid map looks like:

```bash
--atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw"
```

and an invalid map would be any map using names not listed in `--samples`.

---

## Detailed input requirements

### A. H3K27me3 treatment BAM: `--treatment_map`

Used for:
- `nucleosome`
- `array` when `--peak_mode sicer2`
- `run` whenever the selected steps require nucleosome calling or SICER2 peak calling

Format requirements:
- one BAM per sample
- coordinate-sorted
- BAM index recommended
- map format: `sample=/path/to/file.bam;...`

Example:

```bash
--treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam"
```

---

### B. Control BAM: `--control_map` (optional)

Used only when matched control/background BAMs are available for nucleosome calling or SICER2.

Format requirements:
- one BAM per sample
- sample names must match those in `--samples`

Example:

```bash
--control_map "sample1=/path/to/sample1.input.bam;sample2=/path/to/sample2.input.bam;sample3=/path/to/sample3.input.bam"
```

Important behavior:
- in the `nucleosome` step, control must be provided for **either all samples or none**
- mixed control availability across samples is not allowed

---

### C. ATAC bigWig: `--atac_bw_map`

Used for:
- `cluster`
- `run` when clustering is included
- fallback insertion signal for some downstream motif settings

Format requirements:
- one bigWig per sample
- map format: `sample=/path/to/file.bw;...`
- should represent Tn5 cut-site signal or a consistently generated accessibility-centered signal track

Example:

```bash
--atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw;sample3=/path/to/sample3.atac.bw"
```

---

### D. Precomputed nucleosome position table: `--positions_map` (optional)

Used when nucleosome calling has already been completed and you want to start from array or cluster stages without rerunning DANPOS.

Format requirements:
- one file per sample
- should point to a DANPOS position table such as `*.positions.ref_adjust.xls` or `*.positions.xls`

Example:

```bash
--positions_map "sample1=/path/to/sample1.positions.ref_adjust.xls;sample2=/path/to/sample2.positions.ref_adjust.xls"
```

Automatic behavior:
- if `--positions_map` is omitted, NuRepress will try to infer position files from:
  - `OUT_DIR/00_nucleosome_call/pooled/*.positions.ref_adjust.xls`
  - or `OUT_DIR/00_nucleosome_call/pooled/*.positions.xls`

This means that after running the nucleosome step under the same `--out_dir`, later steps usually do not need `--positions_map`.

---

### E. Precomputed broad peaks: `--peak_bed_map` (optional)

Used only when:

```bash
--peak_mode existing
```

Format requirements:
- one BED file per sample
- map format: `sample=/path/to/file.bed;...`

Example:

```bash
--peak_bed_map "sample1=/path/to/sample1.broadPeak.bed;sample2=/path/to/sample2.broadPeak.bed"
```

When `--peak_mode existing` is used, `--peak_bed_map` is required.

When `--peak_mode sicer2` is used, NuRepress ignores `--peak_bed_map` and instead calls peaks from the treatment/control BAM inputs.

---

### F. DPOS file override: `--dpos_map` (optional)

Used mainly by motif analysis in `internal_linker` mode.

Format requirements:
- one file per sample
- map format: `sample=/path/to/file;...`

Example:

```bash
--dpos_map "sample1=/path/to/sample1.positions.ref_adjust.xls;sample2=/path/to/sample2.positions.ref_adjust.xls"
```

Automatic behavior:
- if `--dpos_map` is not provided, NuRepress falls back to `positions_xls`

---

### G. Insertion bigWig override: `--insertion_bw_map` (optional)

Used mainly by motif analysis in `edge_outside` mode.

Format requirements:
- one bigWig per sample
- map format: `sample=/path/to/file.bw;...`

Example:

```bash
--insertion_bw_map "sample1=/path/to/sample1.insertion.bw;sample2=/path/to/sample2.insertion.bw"
```

Automatic behavior:
- if `--insertion_bw_map` is not provided, NuRepress falls back to `--atac_bw_map`

---

### H. Shared reference arguments

These arguments are shared across samples and therefore take a single value rather than a map.

#### `--genome`
Used by:
- `array` when `--peak_mode sicer2`
- `motif`
- `run` when these stages are included

Example:

```bash
--genome hg38
```

#### `--genome_fasta`
Used mainly by:
- `describe`

Optional but recommended when sequence-based summaries are needed.

Example:

```bash
--genome_fasta /path/to/genome.fa
```

#### `--chrom_sizes`
Used by:
- `motif`

Example:

```bash
--chrom_sizes /path/to/genome.chrom.sizes
```

#### `--annotation_gtf`
Used by:
- `describe`
- `motif`

Format requirements:
- GTF or GTF.gz
- should be compatible with downstream gene/TSS annotation needs

Example:

```bash
--annotation_gtf /path/to/annotation.gtf
```

---

### I. Expression matrix: `--expr_tsv`

Used by:
- `score`
- `run` when scoring is included

Format requirements:
- tab-separated file
- one row per gene
- expression columns should correspond to sample-level expression values
- by default, sample columns are expected to be identifiable from the sample prefix

Recommended format:

```text
gene_id	sample1_1	sample1_2	sample2_1	sample2_2
ENSG00000111111	5.20	4.98	2.14	2.31
ENSG00000122222	8.31	8.05	6.77	6.59
ENSG00000133333	1.42	1.18	4.91	5.07
```

Example:

```bash
--expr_tsv /path/to/expression.tsv
```

Optional scoring-related controls:
- `--expr_id_col`
- `--expr_sample_regex_map`

Use these when the expression matrix does not follow the default identifier or sample-column naming conventions.

---

## What each subcommand needs

### 1. `run`

The `run` subcommand executes multiple stages in one controller.

Example:

```bash
run_nurepress_cli.py run \
  --out_dir /path/to/output \
  --steps nucleosome,array,cluster,describe,motif,score \
  --samples sample1,sample2,sample3 \
  --atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw;sample3=/path/to/sample3.atac.bw" \
  --treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam" \
  --control_map "sample1=/path/to/sample1.input.bam;sample2=/path/to/sample2.input.bam;sample3=/path/to/sample3.input.bam" \
  --genome hg38 \
  --genome_fasta /path/to/genome.fa \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --annotation_gtf /path/to/annotation.gtf \
  --expr_tsv /path/to/expression.tsv
```

#### `run` required inputs depend on `--steps`

| Step included | Required inputs |
|---|---|
| `nucleosome` | `--treatment_map` |
| `array` with `--peak_mode sicer2` | `--treatment_map`, `--genome` |
| `array` with `--peak_mode existing` | `--peak_bed_map` |
| `cluster` | `--atac_bw_map` |
| `describe` | `--annotation_gtf` |
| `motif` | `--annotation_gtf`, `--chrom_sizes`, `--genome` |
| `score` | `--expr_tsv`, plus TSS annotation from `describe` or `--tss_anno_tsv` |

Important note:
- in the current integrated CLI, downstream steps can reuse outputs generated earlier under the same `--out_dir`

---

### 2. `nucleosome`

Runs only the DANPOS wrapper.

Required:
- `--out_dir`
- `--samples`
- `--treatment_map`

Optional:
- `--control_map`
- `--nucleosome_command`
- `--nucleosome_threads`
- `--danpos_py`
- `--skip_nucleosome_index`
- `--nucleosome_extra_args`

Example:

```bash
run_nurepress_cli.py nucleosome \
  --out_dir /path/to/output \
  --samples sample1,sample2,sample3 \
  --treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam"
```

Output:
- `00_nucleosome_call/`
- `logs/nucleosome.log`

Key result:
- pooled DANPOS position tables under `00_nucleosome_call/pooled/`

---

### 3. `array`

Runs broad-peak calling and well-phased array identification.

Required:
- `--out_dir`
- `--samples`
- either:
  - `--peak_mode sicer2` with `--treatment_map` and `--genome`
  - or `--peak_mode existing` with `--peak_bed_map`

Optional:
- `--positions_map`
- `--control_map`
- `--dpos_map`
- `--insertion_bw_map`
- `--array_method`
- `--array_k`
- `--min_overlap_frac`
- `--chrom_allow_regex`
- `--array_extra_args`

Example using SICER2:

```bash
run_nurepress_cli.py array \
  --out_dir /path/to/output \
  --samples sample1,sample2,sample3 \
  --treatment_map "sample1=/path/to/sample1.H3K27me3.bam;sample2=/path/to/sample2.H3K27me3.bam;sample3=/path/to/sample3.H3K27me3.bam" \
  --genome hg38 \
  --peak_mode sicer2
```

Example using existing broad peaks:

```bash
run_nurepress_cli.py array \
  --out_dir /path/to/output \
  --samples sample1,sample2 \
  --positions_map "sample1=/path/to/sample1.positions.ref_adjust.xls;sample2=/path/to/sample2.positions.ref_adjust.xls" \
  --peak_bed_map "sample1=/path/to/sample1.broadPeak.bed;sample2=/path/to/sample2.broadPeak.bed" \
  --peak_mode existing
```

Output:
- `01_array_call/<sample>/`
- `logs/array.<sample>.log`

Key results per sample:
- `01_sicer/` if SICER2 is used
- `02_array_raw/`
- `03_array_peak_filtered/`
- `run_metadata.txt`

---

### 4. `cluster`

Runs accessibility-based subtype clustering.

Required:
- `--out_dir`
- `--samples`
- `--atac_bw_map`

Optional:
- `--positions_map`
- `--dpos_map`
- `--insertion_bw_map`
- `--ks`
- `--best_k_source`
- `--sample_order`
- `--cluster_extra_args`

Example:

```bash
run_nurepress_cli.py cluster \
  --out_dir /path/to/output \
  --samples sample1,sample2,sample3 \
  --atac_bw_map "sample1=/path/to/sample1.atac.bw;sample2=/path/to/sample2.atac.bw;sample3=/path/to/sample3.atac.bw"
```

Automatic behavior:
- if `--positions_map` is omitted, NuRepress tries to recover DANPOS position files from the nucleosome output under the same `--out_dir`
- internal intermediate tables needed for clustering are generated automatically inside `intermediate/`

Output:
- `02_cluster/`
- `intermediate/positions_sample_sheet.tsv`
- `intermediate/cluster_sample_sheet.tsv`
- `intermediate/region_stats_sheet.tsv`
- `intermediate/peak_sheet.tsv`
- `logs/cluster.log`

Key results:
- `02_cluster/best_k_selection.tsv`
- `02_cluster/cluster_tables/k{best_k}_cluster_assignment.z1z2.*.tsv`

---

### 5. `describe`

Runs subtype annotation and TSS/promoter association analysis.

Required:
- `--out_dir`
- `--annotation_gtf`

Optional:
- `--genome_fasta`
- `--sample_order`
- `--describe_use_transcript_tss` (default behavior)
- `--describe_gene_tss`
- `--describe_extra_args`

Example:

```bash
run_nurepress_cli.py describe \
  --out_dir /path/to/output \
  --annotation_gtf /path/to/annotation.gtf \
  --genome_fasta /path/to/genome.fa
```

Input dependency:
- this step expects clustering outputs and the automatically generated intermediate files to already exist under the same `--out_dir`

Output:
- `03_describe/`
- `logs/describe.log`

Key result:
- `03_describe/03_TSS_annotation/region_TSS_association.detail.transcriptTSS.tsv`

If `--describe_gene_tss` is used, the output switches to the gene-TSS version.

---

### 6. `motif`

Runs motif region construction and motif enrichment.

Required:
- `--out_dir`
- `--annotation_gtf`
- `--chrom_sizes`
- `--genome`

Optional:
- `--motif_mode`
- `--target_clusters`
- `--tss_near_bp`
- `--cluster_dir_default`
- `--cluster_dir_map`
- `--rel_range`
- `--homer_auto_bg`
- `--build_only`
- `--motif_extra_args`

Example:

```bash
run_nurepress_cli.py motif \
  --out_dir /path/to/output \
  --annotation_gtf /path/to/annotation.gtf \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --genome hg38
```

Mode-specific behavior:

#### `--motif_mode internal_linker`
- uses `dpos_xls`
- falls back to DANPOS position files if no explicit `--dpos_map` is given

#### `--motif_mode edge_outside`
- uses insertion/accessibility signal
- falls back to `--atac_bw_map` if no explicit `--insertion_bw_map` is given

Input dependency:
- this step expects clustering outputs to already exist under the same `--out_dir`

Output:
- `04_motif/`
- `logs/motif.log`

Key subdirectories:
- `04_motif/bed_target/`
- `04_motif/bed_bg/`
- `04_motif/meta/`
- `04_motif/homer_out/`
- `04_motif/homer_logs/`

---

### 7. `score`

Runs candidate repressor scoring.

Required:
- `--out_dir`
- `--expr_tsv`

Optional:
- `--tss_anno_tsv`
- `--sample_order`
- `--expr_id_col`
- `--expr_sample_regex_map`
- `--promoter_flag_col`
- `--restrict_clusters`
- `--gene_id_candidates`
- `--motif_q_cutoff`
- `--min_group_n`
- `--top_n_plot`
- `--score_extra_args`

Example:

```bash
run_nurepress_cli.py score \
  --out_dir /path/to/output \
  --expr_tsv /path/to/expression.tsv
```

Automatic behavior:
- if `--tss_anno_tsv` is omitted, NuRepress tries to use:
  - `03_describe/03_TSS_annotation/region_TSS_association.detail.transcriptTSS.tsv`

Output:
- `05_score/`
- `logs/score.log`

This directory contains the final TF ranking tables and related plotting outputs.

---

## Standard output structure

NuRepress writes outputs into a standardized directory tree under `--out_dir`:

```text
OUT_DIR/
├── 00_nucleosome_call/
├── 01_array_call/
├── 02_cluster/
├── 03_describe/
├── 04_motif/
├── 05_score/
├── intermediate/
├── logs/
└── pipeline_run_metadata.txt
```

### Main output directories

#### `00_nucleosome_call/`
DANPOS wrapper outputs and pooled nucleosome position tables.

#### `01_array_call/`
Per-sample well-phased array calling results.

#### `02_cluster/`
Subtype clustering results, best-k selection, and cluster assignment tables.

#### `03_describe/`
Annotation and TSS association outputs.

#### `04_motif/`
Motif target/background BEDs, metadata, HOMER outputs, and motif logs.

#### `05_score/`
Ranked candidate TF tables and score-related plots.

#### `intermediate/`
Automatically generated internal tables used to bridge steps. These are created by the integrated CLI and do **not** need to be prepared manually.

#### `logs/`
Step-level execution logs.

#### `pipeline_run_metadata.txt`
A simple key-value metadata summary written by the `run` controller.

---

## Most important outputs at a glance

If you only want the main biological deliverables, focus on these files and directories:

- **nucleosome position outputs**  
  `00_nucleosome_call/pooled/*.positions.ref_adjust.xls`

- **well-phased array outputs**  
  `01_array_call/<sample>/03_array_peak_filtered/`

- **best clustering solution**  
  `02_cluster/best_k_selection.tsv`

- **array subtype assignments**  
  `02_cluster/cluster_tables/k{best_k}_cluster_assignment.z1z2.*.tsv`

- **TSS/promoter association table**  
  `03_describe/03_TSS_annotation/region_TSS_association.detail.transcriptTSS.tsv`

- **motif enrichment results**  
  `04_motif/homer_out/`

- **final TF ranking results**  
  `05_score/`

---

## Automatic step-to-step reuse

A major feature of the integrated CLI is that later steps can automatically reuse outputs produced earlier under the same `--out_dir`.

Examples:

- `array` can infer DANPOS position files from `00_nucleosome_call/`
- `cluster` can infer position files and will automatically generate the internal intermediate tables it needs
- `score` can automatically use the TSS annotation produced by `describe`

Because of this design, a typical workflow is:
1. run upstream steps once
2. adjust parameters
3. rerun only the downstream steps that need to change

---



## Troubleshooting

### 1. DANPOS position files cannot be found downstream
Check:
- whether `nucleosome` has already been run under the same `--out_dir`
- whether `00_nucleosome_call/pooled/` contains `*.positions.ref_adjust.xls` or `*.positions.xls`
- whether `--positions_map` needs to be provided explicitly

### 2. `array` fails under `--peak_mode existing`
Check:
- whether `--peak_bed_map` was provided
- whether every sample in `--samples` has a corresponding BED file

### 3. `array` fails under `--peak_mode sicer2`
Check:
- whether `--treatment_map` was provided
- whether `--genome` was provided
- whether `sicer` is installed and on `PATH`

### 4. `cluster` cannot proceed
Check:
- whether `--atac_bw_map` was provided
- whether nucleosome and array outputs already exist under the same `--out_dir`
- whether the previous step produced the expected per-sample array outputs

### 5. `describe` cannot find upstream inputs
Check:
- whether clustering has already been completed under the same `--out_dir`
- whether the `intermediate/region_stats_sheet.tsv` file was generated

### 6. `motif` produces no useful motif results
Check:
- whether `findMotifsGenome.pl` is on `PATH`
- whether the genome package has been configured in HOMER
- whether target BED files contain enough valid regions

### 7. `score` reports no usable genes
Check:
- whether the TSS annotation exists
- whether expression gene IDs are compatible with the annotation-derived gene identifiers
- whether `--expr_id_col` or `--expr_sample_regex_map` should be specified
- whether the promoter flag column is correctly set

---

## Contact

If you have any questions, feel free to contact:

**xiang.qianming@hsc.pku.edu.cn**
