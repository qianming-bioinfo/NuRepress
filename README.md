# NuRepress: inferring transcriptional repressors from phased nucleosome architecture

NuRepress is an end-to-end workflow for identifying well-phased nucleosome arrays in repressive chromatin, stratifying array subtypes by accessibility patterns, performing motif enrichment on subtype-defined regions, and prioritizing candidate transcriptional repressors by integrating chromatin structure, sequence evidence, and transcriptional effects.

## Workflow

![NuRepress workflow](docs/figures/workflow_overview.png)


## Abstract

This repository provides a connected implementation of the NuRepress analysis workflow. The pipeline can start from ChIP/input BAM files, call nucleosome positions through DANPOS3, identify well-phased arrays within broad repressive domains, cluster arrays into accessibility-defined subtypes, describe subtype organization, run motif enrichment on internal linker or edge-outside regions, and score candidate repressors using subtype-associated expression shifts and motif support.

The repository is organized so that each stage can run independently, but the recommended usage is through the top-level orchestrator `run_full_nurepress_pipeline.py`, which keeps all intermediate paths consistent and allows downstream steps to inherit the optimal clustering result (`best_k`) automatically.

---

## Software dependencies

This repository is primarily intended for **Linux** or Linux-like HPC environments, because the workflow depends on Bash wrappers, common command-line genomics tools, and R/Python executables available on `$PATH`.

### External tools

Required for the full workflow:

- `bash`
- `python` (recommended: Python 3.10+)
- `Rscript` (recommended: R 4.3+)
- `samtools`
- `danpos.py` from DANPOS3
- `sicer` from SICER2, if `--peak_mode sicer2` is used
- `findMotifsGenome.pl` from HOMER, if motif enrichment is run

Conditionally required:

- `bedtools`  
  Needed by SICER2 when treatment/control input files are BAM rather than BED.
- HOMER genome packages  
  HOMER installation and genome configuration are managed through `configureHomer.pl`.

### Python packages

The bundled Python scripts use:

- `pandas`
- `numpy`
- `tqdm`
- `pyBigWig`

For SICER2 itself, the official repository notes that Python 3, `numpy`, `scipy`, and a C compiler are required, and BAM input additionally requires `bedtools`.

### R packages

The bundled R scripts use:

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
- `IRanges`
- `S4Vectors`
- `rtracklayer`
- `Rsamtools`
- `Biostrings`

---

## Recommended environment setup

### Option 1. Conda / mamba environment

```bash
mamba env create -f environment.yml
mamba activate nurepress
```

This environment installs the Python stack, common command-line dependencies, and the R/Bioconductor packages required by the bundled scripts.

### Option 2. Install Python packages only

```bash
pip install -r requirements.txt
```

This installs the Python packages used directly by the bundled Python scripts, but it does **not** install DANPOS3, SICER2, HOMER, `samtools`, `bedtools`, or the R/Bioconductor stack.

### Option 3. Install R packages inside an existing R environment

```bash
Rscript install_r_packages.R
```

This installs the CRAN and Bioconductor packages required by the bundled R scripts.

### Tool-specific installation notes

#### DANPOS3

The DANPOS3 repository documents installation through `git clone` plus `pip install -r requirements.txt`, and lists Python 3.7.6, R 4.0.1, `samtools 1.7`, `numpy`, `pysam`, and `rpy2` in the tested environment.

#### SICER2

The SICER2 repository states that the easiest installation is:

```bash
pip install SICER2
```

and notes requirements including Python 3, `numpy`, `scipy`, a C compiler, and optionally `bedtools` for BAM input.

#### HOMER

HOMER is installed/configured through `configureHomer.pl`, and additional genome packages are managed through the same configuration interface.

---


## Workflow summary

### Step 0. Nucleosome calling

`call_nucleosomes.sh` wraps DANPOS3 in command-line form. DANPOS3 is a Python 3 update of the original DANPOS codebase and its repository documents Python, R, `samtools`, and Python-library requirements.

### Step 1. Broad-peak calling and well-phased array identification

`call_peak_and_identify_array.py` supports two peak modes:

- `sicer2`: call broad peaks with SICER2
- `existing`: use precomputed broad-peak BED files

For SICER2, the official command is `sicer`; required inputs are treatment and species, and BAM input additionally requires `bedtools`. Official docs also describe default parameters such as `-w 200`, `-g 600`, and `-fdr 0.01`.

Array identification supports two methods:

- `seedextend` (default)
- `mergewindow`

Both methods use `k=3` by default in this repository.

### Step 2. Array subtype clustering

`cluster_array_subtype.R` uses accessibility-derived features around each array to evaluate multiple `k` values and then selects the best clustering solution by silhouette score. Downstream steps read `best_k_selection.tsv` and automatically inherit the selected `best_k`.

### Step 3. Subtype description and annotation

`describe_array_subtype.R` can summarize subtype composition, perform TSS annotation, calculate optional GC-content summaries, and describe optional union-region remodeling behavior. Bioconductor packages such as `GenomicRanges`, `rtracklayer`, `Rsamtools`, and `Biostrings` are used for genomic interval and sequence operations; these packages are documented through Bioconductor.

### Step 4. Motif enrichment

`run_array_motif_enrichment.py` directly connects to the clustering output directory. It reads `best_k_selection.tsv`, locates the matching `k{best_k}_cluster_assignment...tsv`, constructs target/background regions, and optionally runs HOMER.

HOMER's `findMotifsGenome.pl` is the official wrapper for genome-region motif analysis, and when no explicit background is provided, HOMER automatically selects background sequences.

### Step 5. Candidate repressor scoring

`score_repressor_candidates.R` parses HOMER known-motif results, combines them with subtype-specific expression effects, and produces ranked TF tables and plots.

---

## Input table design: sample sheet

Two template files are provided:

- `sample_sheet.minimal.tsv`
- `sample_sheet.full.tsv`

A field-by-field explanation is available in `sample_sheet_templates_README.md`.

### Minimum required columns for the top-level pipeline

The top-level pipeline always requires:

- `sample`
- `atac_bw`

### Common optional columns

- `positions_xls`
- `treatment`
- `control`
- `peak_bed`
- `species`
- `genome`
- `dpos_xls`
- `insertion_bw`
- `genome_fasta`
- `chrom_sizes`
- `annotation_gtf`
- `motif_annotation_file`

### Shared-reference inheritance rules

If you do not pass some shared references on the command line, `run_full_nurepress_pipeline.py` will inherit them from the sample sheet **only when exactly one non-empty shared value exists across all samples**.

This applies to:

- `--gtf` ← `annotation_gtf`
- `--ref_fa` ← `genome_fasta`
- `--motif_annotation_file` ← `motif_annotation_file`, otherwise `annotation_gtf`
- `--chrom_sizes` ← `chrom_sizes`
- `--genome` ← `genome`, otherwise `species`

---

## Quick start

### A. Full run from treatment/control BAMs

```bash
python run_full_nurepress_pipeline.py   --sample_sheet /path/to/sample_sheet.full.tsv   --out_dir /path/to/pipeline_out   --steps nucleosome,array,cluster,describe,motif,score   --call_nucleosomes_script /path/to/call_nucleosomes.sh   --danpos_py /path/to/danpos.py   --nucleosome_threads 32   --nucleosome_extra_args "-m 1 --mifrsz 80 --mafrsz 800"   --peak_mode sicer2   --array_method seedextend   --array_k 3   --ks 2,3,4,5,6   --best_k_source within_sample   --gtf /path/to/annotation.gtf.gz   --motif_annotation_file /path/to/annotation.gtf.gz   --chrom_sizes /path/to/genome.chrom.sizes   --genome hg38   --motif_mode internal_linker   --target_clusters ALL   --expr_tsv /path/to/expression_matrix.tsv
```

### B. Skip DANPOS and start from existing `positions_xls`

```bash
python run_full_nurepress_pipeline.py   --sample_sheet /path/to/sample_sheet.full.tsv   --out_dir /path/to/pipeline_out   --steps array,cluster,describe,motif,score   --peak_mode existing   --array_method seedextend   --array_k 3   --ks 2,3,4,5,6   --expr_tsv /path/to/expression_matrix.tsv
```

### C. Run motif-region construction only

```bash
python run_full_nurepress_pipeline.py   --sample_sheet /path/to/sample_sheet.full.tsv   --out_dir /path/to/pipeline_out   --steps motif   --motif_mode edge_outside   --target_clusters ALL   --build_only
```

---

## Step-by-step input requirements

### `nucleosome`

Required per sample:
- `treatment`
- `control`

Outputs passed downstream:
- per-sample `*.positions.ref_adjust.xls`, rediscovered under `00_nucleosome_call/pooled/`

### `array`

Required per sample:
- `positions_xls` (either pre-existing or generated by `nucleosome`)
- and either:
  - `treatment` + `species`, when `--peak_mode sicer2`, or
  - `peak_bed`, when `--peak_mode existing`

### `cluster`

Required per sample:
- filtered array BED generated by the array step
- `atac_bw`

### `describe`

Required:
- cluster output directory
- `--gtf` or a single shared `annotation_gtf`

Optional but useful:
- `--ref_fa` or shared `genome_fasta`

### `motif`

Required:
- cluster output directory
- motif annotation file (`--motif_annotation_file` or shared value)
- `--chrom_sizes`
- `--genome`

Additional mode-specific requirements:
- `internal_linker`: `positions_xls` or `dpos_xls`
- `edge_outside`: `atac_bw` or `insertion_bw`

### `score`

Required:
- motif output directory
- expression matrix (`--expr_tsv`)
- TSS annotation TSV (either passed explicitly by `--tss_anno_tsv` or inherited from the describe step)

---

## Output directory structure

```text
<out_dir>/
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

### Key files by stage

#### `00_nucleosome_call/`
- DANPOS wrapper output
- DANPOS logs in `00_nucleosome_call/logs/`
- `pooled/*.positions.ref_adjust.xls`

#### `01_array_call/<sample>/`
- `01_sicer/` if `--peak_mode sicer2`
- `02_array_raw/`
- `03_array_peak_filtered/`
- `run_metadata.txt`

#### `02_cluster/`
- `best_k_selection.tsv`
- `cluster_tables/k{best_k}_cluster_assignment.z1z2.*.tsv`
- silhouette summaries and clustering figures

#### `03_describe/`
- `03_TSS_annotation/region_TSS_association.detail.transcriptTSS.tsv`
- optional GC-content and remodeling summaries

#### `04_motif/`
- `bed_target/`
- `bed_bg/`
- `meta/`
- `homer_out/`
- `homer_logs/`

#### `05_score/`
- TF ranking tables by sample and cluster
- priority summaries
- plotting outputs for top candidates

---

## Resume and partial rerun strategy

The repository is designed so that steps can be rerun independently.

### Recommended practice

- Use the full pipeline script for the first run.
- For revisions, rerun only the necessary downstream steps.
- Keep each run in a dedicated `out_dir` to preserve consistent provenance.

### Typical examples

Re-run only clustering and all downstream steps:

```bash
python run_full_nurepress_pipeline.py   --sample_sheet /path/to/sample_sheet.full.tsv   --out_dir /path/to/pipeline_out   --steps cluster,describe,motif,score   --expr_tsv /path/to/expression_matrix.tsv
```

Re-run only motif and scoring after adjusting motif settings:

```bash
python run_full_nurepress_pipeline.py   --sample_sheet /path/to/sample_sheet.full.tsv   --out_dir /path/to/pipeline_out   --steps motif,score   --motif_mode edge_outside   --target_clusters C2   --expr_tsv /path/to/expression_matrix.tsv
```

---


## Repository contents

### Top-level scripts

- `run_full_nurepress_pipeline.py`  
  Top-level orchestrator for the complete workflow.
- `call_nucleosomes.sh`  
  Wrapper around `danpos.py dpos` with BAM integrity checks, indexing, metadata logging, and pass-through DANPOS arguments.
- `call_peak_and_identify_array.py`  
  Calls SICER2 or uses existing peaks, then identifies well-phased arrays.
- `cluster_array_subtype.R`  
  Computes subtype-defining accessibility features and selects the optimal number of clusters using silhouette scores.
- `describe_array_subtype.R`  
  Generates subtype summaries, TSS annotation, optional GC-content summaries, and optional union-region remodeling summaries.
- `run_array_motif_enrichment.py`  
  Builds motif target/background BED files and optionally runs HOMER.
- `score_repressor_candidates.R`  
  Integrates motif support, subtype-specific regulatory effect, and specificity into TF ranking tables.

### Template and environment files

- `sample_sheet.minimal.tsv`
- `sample_sheet.full.tsv`
- `sample_sheet_templates_README.md`
- `environment.yml`
- `requirements.txt`
- `install_r_packages.R`

---

## Troubleshooting

### 1. `positions_xls` cannot be found downstream

Check whether:
- `nucleosome` was included in `--steps`, and
- `00_nucleosome_call/pooled/` contains the expected `*.positions.ref_adjust.xls` files.

### 2. SICER2 fails on BAM input

Make sure:
- `sicer` is installed,
- `bedtools` is installed,
- `species` is correctly provided.

### 3. HOMER runs but no motif result is produced

Check whether:
- `findMotifsGenome.pl` is on `$PATH`,
- the requested genome package has been configured in HOMER,
- target BED files contain enough regions.

### 4. The scoring step reports no usable genes

Check:
- promoter annotation flags in the TSS table,
- gene identifier compatibility between TSS annotation and expression matrix,
- whether `--expr_id_col` or `--expr_sample_regex_map` needs to be specified.

---

## Citation

If you use this repository in a manuscript, please cite:

1. The NuRepress manuscript or software paper, once available.
2. The underlying external tools used in your workflow, especially DANPOS, SICER/SICER2, HOMER, and the Bioconductor packages used for genomic interval and annotation handling.

A suggested repository citation block can be added here once the final manuscript title, author list, DOI, and release tag are available.

---

## Notes

- Public examples in this repository intentionally use generic placeholders such as `sample1`, `sample2`, and `/path/to/...`.
- For reproducibility, keep your final repository release synchronized with the exact script bundle used in the manuscript.
