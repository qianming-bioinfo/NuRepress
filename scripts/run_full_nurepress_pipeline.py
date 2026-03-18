#!/usr/bin/env python3
import argparse
import glob
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def split_csv(x: Optional[str]) -> List[str]:
    if x is None:
        return []
    return [i.strip() for i in str(x).split(",") if i.strip()]


def shlex_list(x: Optional[str]) -> List[str]:
    if x is None or str(x).strip() == "":
        return []
    return shlex.split(str(x))


class PipelineError(RuntimeError):
    pass


def detect_delim(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    with open(path, "r", newline="") as f:
        head = f.readline()
    if head.count("\t") >= head.count(","):
        return "\t"
    return ","


def read_sample_sheet(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise PipelineError(f"Sample sheet not found: {path}")
    delim = detect_delim(path)
    df = pd.read_csv(path, sep=delim, dtype=str).fillna("")
    df.columns = [c.strip() for c in df.columns]

    sample_col = None
    for cand in ["sample", "cell_type", "sample_id"]:
        if cand in df.columns:
            sample_col = cand
            break
    if sample_col is None:
        raise PipelineError("Sample sheet must contain one of: sample, cell_type, sample_id")
    df = df.rename(columns={sample_col: "sample"})

    required = ["sample", "atac_bw"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise PipelineError(
            "Sample sheet missing required columns: " + ", ".join(missing)
            + ". Required columns are: sample, atac_bw"
        )

    optional_cols = [
        "positions_xls",
        "treatment",
        "control",
        "peak_bed",
        "species",
        "dpos_xls",
        "insertion_bw",
        "genome",
        "genome_fasta",
        "chrom_sizes",
        "annotation_gtf",
        "motif_annotation_file",
    ]
    for c in optional_cols:
        if c not in df.columns:
            df[c] = ""

    df["sample"] = df["sample"].astype(str).str.strip()
    if (df["sample"] == "").any():
        raise PipelineError("Sample sheet contains empty sample names")
    if df["sample"].duplicated().any():
        dup = df.loc[df["sample"].duplicated(), "sample"].tolist()
        raise PipelineError(f"Duplicated sample names in sample sheet: {dup}")

    for col in ["sample", *optional_cols, "atac_bw"]:
        df[col] = df[col].astype(str).str.strip()

    return df


def repo_script_path(base_dir: Path, name: str) -> Path:
    p = base_dir / name
    if not p.exists():
        raise PipelineError(f"Required script not found: {p}")
    return p


def run_cmd(cmd: List[str], log_path: Path, cwd: Optional[Path] = None):
    ensure_dir(log_path.parent)
    eprint("[RUN]", " ".join(shlex.quote(x) for x in cmd))
    with open(log_path, "w") as log:
        log.write("COMMAND=\n")
        log.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=str(cwd) if cwd else None)
    if proc.returncode != 0:
        raise PipelineError(f"Command failed with exit code {proc.returncode}. See log: {log_path}")


def parse_kv_file(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.exists():
        return out
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def discover_single_file(pattern: str) -> str:
    hits = sorted(glob.glob(pattern))
    if len(hits) == 0:
        raise PipelineError(f"No file matched pattern: {pattern}")
    return hits[0]


def write_pipeline_metadata(path: Path, rows: List[Tuple[str, str]]):
    ensure_dir(path.parent)
    with open(path, "w") as f:
        for k, v in rows:
            f.write(f"{k}={v}\n")


def first_unique_nonempty(df: pd.DataFrame, column: str) -> str:
    if column not in df.columns:
        return ""
    vals = sorted({str(x).strip() for x in df[column].tolist() if str(x).strip()})
    if len(vals) == 0:
        return ""
    if len(vals) == 1:
        return vals[0]
    raise PipelineError(
        f"Column '{column}' contains multiple non-empty values. "
        f"Please pass the corresponding command-line argument explicitly."
    )


def resolve_shared_value(cli_value: Optional[str], samples: pd.DataFrame, column: str, fallback_columns: Optional[List[str]] = None) -> str:
    if cli_value and str(cli_value).strip():
        return str(cli_value).strip()
    tried = [column] + (fallback_columns or [])
    for col in tried:
        val = first_unique_nonempty(samples, col)
        if val:
            return val
    return ""


def discover_danpos_position_file(nucleosome_out_dir: Path, sample: str) -> str:
    pooled = nucleosome_out_dir / "pooled"
    if not pooled.exists():
        raise PipelineError(f"DANPOS pooled directory not found: {pooled}")

    patterns = [
        str(pooled / f"{sample}*.positions.ref_adjust.xls"),
        str(pooled / f"{sample}*.positions.xls"),
        str(pooled / f"{sample}*positions.ref_adjust.xls"),
        str(pooled / f"{sample}*positions.xls"),
    ]
    hits: List[str] = []
    for pat in patterns:
        hits.extend(sorted(glob.glob(pat)))
    hits = sorted(dict.fromkeys(hits))
    if not hits:
        raise PipelineError(
            f"Cannot find DANPOS position output for sample '{sample}' under: {pooled}. "
            f"Expected files like {sample}*.positions.ref_adjust.xls"
        )
    return hits[0]


def enrich_positions_from_nucleosome(samples: pd.DataFrame, nucleosome_out_dir: Path) -> pd.DataFrame:
    df = samples.copy()
    positions: List[str] = []
    for row in df.itertuples(index=False):
        if str(row.positions_xls).strip():
            positions.append(str(row.positions_xls).strip())
        else:
            positions.append(discover_danpos_position_file(nucleosome_out_dir, row.sample))
    df["positions_xls"] = positions
    needs_dpos = df["dpos_xls"].astype(str).str.strip() == ""
    df.loc[needs_dpos, "dpos_xls"] = df.loc[needs_dpos, "positions_xls"]
    needs_ins = df["insertion_bw"].astype(str).str.strip() == ""
    df.loc[needs_ins, "insertion_bw"] = df.loc[needs_ins, "atac_bw"]
    return df


def discover_array_outputs(sample_out_dir: Path) -> Tuple[str, str, str]:
    meta = parse_kv_file(sample_out_dir / "run_metadata.txt")
    bed = meta.get("filtered_bed", "")
    stats = meta.get("filtered_stats", "")
    peak = meta.get("peak_bed", "")

    if not bed:
        bed = discover_single_file(str(sample_out_dir / "03_array_peak_filtered" / "*.bed"))
    if not stats:
        stats = discover_single_file(str(sample_out_dir / "03_array_peak_filtered" / "*.stats.tsv"))
    if not peak:
        peak = discover_single_file(str(sample_out_dir / "01_sicer" / "**" / "*.bed"))

    return bed, stats, peak


def build_cluster_sample_sheet(samples: pd.DataFrame, array_root: Path, out_path: Path) -> pd.DataFrame:
    rows = []
    for row in samples.itertuples(index=False):
        sample = row.sample
        sample_dir = array_root / sample
        bed, stats_tsv, peak_bed = discover_array_outputs(sample_dir)
        rows.append({
            "cell_type": sample,
            "regions_bed": bed,
            "atac_bw": row.atac_bw,
            "stats_tsv": stats_tsv,
            "peak_bed": peak_bed,
            "positions_xls": row.positions_xls,
            "dpos_xls": row.dpos_xls if row.dpos_xls else row.positions_xls,
            "insertion_bw": row.insertion_bw if row.insertion_bw else row.atac_bw,
        })
    out_df = pd.DataFrame(rows)
    ensure_dir(out_path.parent)
    out_df.to_csv(out_path, sep="\t", index=False)
    return out_df


def locate_best_k_cluster_tsv(cluster_run_dir: Path, fit_space_preference: str = "auto") -> Tuple[int, str]:
    bk_path = cluster_run_dir / "best_k_selection.tsv"
    if not bk_path.exists():
        raise PipelineError(f"Missing best_k_selection.tsv under: {cluster_run_dir}")
    bk = pd.read_csv(bk_path, sep="\t")
    if "best_k" not in bk.columns or bk.empty:
        raise PipelineError(f"Invalid best_k_selection.tsv: {bk_path}")
    best_k = int(bk.loc[0, "best_k"])
    fit_space = None
    if "fit_space" in bk.columns and pd.notna(bk.loc[0, "fit_space"]):
        fit_space = str(bk.loc[0, "fit_space"])

    cluster_dir = cluster_run_dir / "cluster_tables"
    if not cluster_dir.exists():
        raise PipelineError(f"Missing cluster_tables under: {cluster_run_dir}")

    candidates: List[Path] = []
    if fit_space_preference != "auto":
        candidates.append(cluster_dir / f"k{best_k}_cluster_assignment.z1z2.{fit_space_preference}.tsv")
    if fit_space:
        candidates.append(cluster_dir / f"k{best_k}_cluster_assignment.z1z2.{fit_space}.tsv")
    candidates.extend(sorted(cluster_dir.glob(f"k{best_k}_cluster_assignment.z1z2.*.tsv")))

    for c in candidates:
        if Path(c).exists():
            return best_k, str(c)
    raise PipelineError(f"Cannot locate best-k cluster assignment TSV under: {cluster_dir}")


def infer_tss_output(describe_out_dir: Path, use_transcript_tss: bool = True) -> str:
    suffix = "transcriptTSS" if use_transcript_tss else "geneTSS"
    path = describe_out_dir / "03_TSS_annotation" / f"region_TSS_association.detail.{suffix}.tsv"
    if not path.exists():
        raise PipelineError(f"Expected TSS annotation output not found: {path}")
    return str(path)


def bool_str(x: bool) -> str:
    return "true" if x else "false"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "End-to-end orchestration for the NuRepress workflow: optional DANPOS nucleosome calling, "
            "array calling/peak filtering, subtype clustering, optional subtype description/TSS annotation, "
            "motif enrichment, and repressor scoring."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io = ap.add_argument_group("Input / output")
    io.add_argument(
        "--sample_sheet",
        required=True,
        help=(
            "TSV/CSV sample sheet. Required columns: sample, atac_bw. "
            "Supported optional columns: positions_xls, treatment, control, peak_bed, species, dpos_xls, insertion_bw, "
            "genome, genome_fasta, chrom_sizes, annotation_gtf, motif_annotation_file."
        ),
    )
    io.add_argument("--out_dir", required=True, help="Pipeline output directory")
    io.add_argument(
        "--steps",
        default="nucleosome,array,cluster,describe,motif,score",
        help="Comma-separated steps to run. Choices: nucleosome,array,cluster,describe,motif,score",
    )

    tools = ap.add_argument_group("Executables and script paths")
    tools.add_argument("--python", default=sys.executable, help="Python executable")
    tools.add_argument("--rscript", default="Rscript", help="Rscript executable")
    tools.add_argument("--bash", default="bash", help="Bash executable")
    tools.add_argument("--call_nucleosomes_script", default=None, help="Path to call_nucleosomes.sh")
    tools.add_argument("--call_array_script", default=None, help="Path to call_peak_and_identify_array.py")
    tools.add_argument("--cluster_script", default=None, help="Path to cluster_array_subtype.R")
    tools.add_argument("--describe_script", default=None, help="Path to describe_array_subtype.R")
    tools.add_argument("--motif_script", default=None, help="Path to run_array_motif_enrichment.py")
    tools.add_argument("--score_script", default=None, help="Path to score_repressor_candidates.R")

    nucg = ap.add_argument_group("Nucleosome-calling step (DANPOS wrapper)")
    nucg.add_argument("--nucleosome_out_subdir", default="00_nucleosome_call", help="Subdirectory for DANPOS wrapper output")
    nucg.add_argument("--nucleosome_command", default="dpos", help="DANPOS subcommand passed to call_nucleosomes.sh")
    nucg.add_argument("--nucleosome_threads", type=int, default=1, help="Threads passed to call_nucleosomes.sh")
    nucg.add_argument("--danpos_py", default="danpos.py", help="Path to danpos.py passed to call_nucleosomes.sh")
    nucg.add_argument("--skip_nucleosome_index", action="store_true", default=False, help="Pass --skip-index to call_nucleosomes.sh")
    nucg.add_argument("--nucleosome_extra_args", default="", help="Extra arguments passed to call_nucleosomes.sh after --")

    arrayg = ap.add_argument_group("Array-calling step")
    arrayg.add_argument("--array_method", choices=["seedextend", "mergewindow"], default="seedextend")
    arrayg.add_argument("--array_k", type=int, default=3)
    arrayg.add_argument("--peak_mode", choices=["sicer2", "existing"], default="sicer2")
    arrayg.add_argument("--min_overlap_frac", type=float, default=0.30)
    arrayg.add_argument("--chrom_allow_regex", default=r"^chr([0-9]+|X|Y)$")
    arrayg.add_argument("--array_extra_args", default="", help="Extra arguments passed to call_peak_and_identify_array.py")

    clg = ap.add_argument_group("Cluster step")
    clg.add_argument("--sample_order", default=None, help="Comma-separated sample order")
    clg.add_argument("--ks", default="2,3,4,5,6")
    clg.add_argument("--best_k_source", default="within_sample")
    clg.add_argument("--cluster_extra_args", default="", help="Extra arguments passed to cluster_array_subtype.R")

    desg = ap.add_argument_group("Describe / TSS annotation step")
    desg.add_argument("--gtf", default=None, help="GTF/GFF used for TSS annotation in describe step")
    desg.add_argument("--ref_fa", default=None, help="Indexed reference FASTA for optional GC module in describe step")
    desg.add_argument("--describe_use_transcript_tss", action="store_true", default=True)
    desg.add_argument("--describe_gene_tss", dest="describe_use_transcript_tss", action="store_false")
    desg.add_argument("--describe_extra_args", default="", help="Extra arguments passed to describe_array_subtype.R")

    motifg = ap.add_argument_group("Motif step")
    motifg.add_argument("--motif_annotation_file", default=None, help="Annotation file for motif step. Supports BED/BED.gz/GTF/GTF.gz")
    motifg.add_argument("--chrom_sizes", default=None)
    motifg.add_argument("--genome", default=None)
    motifg.add_argument("--motif_mode", choices=["internal_linker", "edge_outside"], default="internal_linker")
    motifg.add_argument("--target_clusters", default="ALL")
    motifg.add_argument("--tss_near_bp", type=int, default=3000)
    motifg.add_argument("--cluster_dir_default", choices=["neutral", "inside_low", "inside_high"], default="neutral")
    motifg.add_argument("--cluster_dir_map", default=None, help="Example: C1=inside_low;C2=neutral")
    motifg.add_argument("--rel_range", default="+275:+475")
    motifg.add_argument("--homer_auto_bg", action="store_true", default=False)
    motifg.add_argument("--build_only", action="store_true", default=False)
    motifg.add_argument("--motif_extra_args", default="", help="Extra arguments passed to run_array_motif_enrichment.py")

    scoreg = ap.add_argument_group("Scoring step")
    scoreg.add_argument("--expr_tsv", default=None, help="Expression matrix TSV for scoring step")
    scoreg.add_argument("--tss_anno_tsv", default=None, help="Optional existing TSS annotation TSV. Overrides describe-step output when provided")
    scoreg.add_argument("--promoter_flag_col", default="overlap_TSSwin_2000")
    scoreg.add_argument("--restrict_clusters", default="ALL")
    scoreg.add_argument("--expr_id_col", default=None)
    scoreg.add_argument("--expr_sample_regex_map", default=None)
    scoreg.add_argument("--gene_id_candidates", default="nearest_TSS_gene_name,gene_symbol_mapped,nearest_TSS_id")
    scoreg.add_argument("--motif_q_cutoff", type=float, default=0.05)
    scoreg.add_argument("--min_group_n", type=int, default=5)
    scoreg.add_argument("--top_n_plot", type=int, default=20)
    scoreg.add_argument("--score_extra_args", default="", help="Extra arguments passed to score_repressor_candidates.R")

    return ap


def main():
    args = build_parser().parse_args()

    steps = [s.lower() for s in split_csv(args.steps)]
    allowed_steps = {"nucleosome", "array", "cluster", "describe", "motif", "score"}
    unknown_steps = [s for s in steps if s not in allowed_steps]
    if unknown_steps:
        raise PipelineError(f"Unknown steps: {unknown_steps}")

    sample_sheet = Path(args.sample_sheet).resolve()
    samples = read_sample_sheet(sample_sheet)

    out_dir = Path(args.out_dir).resolve()
    ensure_dir(out_dir)
    log_dir = out_dir / "logs"
    ensure_dir(log_dir)

    script_base = Path(__file__).resolve().parent
    call_nucleosomes_script = Path(args.call_nucleosomes_script).resolve() if args.call_nucleosomes_script else repo_script_path(script_base, "call_nucleosomes.sh")
    call_array_script = Path(args.call_array_script).resolve() if args.call_array_script else repo_script_path(script_base, "call_peak_and_identify_array.py")
    cluster_script = Path(args.cluster_script).resolve() if args.cluster_script else repo_script_path(script_base, "cluster_array_subtype.R")
    describe_script = Path(args.describe_script).resolve() if args.describe_script else repo_script_path(script_base, "describe_array_subtype.R")
    motif_script = Path(args.motif_script).resolve() if args.motif_script else repo_script_path(script_base, "run_array_motif_enrichment.py")
    score_script = Path(args.score_script).resolve() if args.score_script else repo_script_path(script_base, "score_repressor_candidates.R")

    resolved_gtf = resolve_shared_value(args.gtf, samples, "annotation_gtf")
    resolved_ref_fa = resolve_shared_value(args.ref_fa, samples, "genome_fasta")
    resolved_motif_annotation = resolve_shared_value(args.motif_annotation_file, samples, "motif_annotation_file", fallback_columns=["annotation_gtf"])
    resolved_chrom_sizes = resolve_shared_value(args.chrom_sizes, samples, "chrom_sizes")
    resolved_genome = resolve_shared_value(args.genome, samples, "genome", fallback_columns=["species"])

    nucleosome_root = out_dir / args.nucleosome_out_subdir
    array_root = out_dir / "01_array_call"
    cluster_root = out_dir / "02_cluster"
    describe_root = out_dir / "03_describe"
    motif_root = out_dir / "04_motif"
    score_root = out_dir / "05_score"

    intermediate_dir = out_dir / "intermediate"
    ensure_dir(intermediate_dir)
    positions_sample_sheet = intermediate_dir / "sample_sheet.with_positions.tsv"
    region_stats_sheet = intermediate_dir / "region_stats_sheet.tsv"
    peak_sheet = intermediate_dir / "peak_sheet.tsv"
    cluster_sample_sheet = intermediate_dir / "cluster_sample_sheet.tsv"

    if "nucleosome" in steps:
        missing_treat = samples["treatment"].astype(str).str.strip() == ""
        if missing_treat.any():
            bad = samples.loc[missing_treat, "sample"].tolist()
            raise PipelineError(f"The nucleosome step requires treatment for all samples. Missing in: {bad}")
        missing_bg = samples["control"].astype(str).str.strip() == ""
        if missing_bg.any():
            bad = samples.loc[missing_bg, "sample"].tolist()
            raise PipelineError(f"The nucleosome step requires control for all samples. Missing in: {bad}")

        input_pairs = ",".join(f"{r.sample}:{r.treatment}" for r in samples.itertuples(index=False))
        bg_pairs = ",".join(f"{r.sample}:{r.control}" for r in samples.itertuples(index=False))
        ensure_dir(nucleosome_root)
        cmd = [
            args.bash, str(call_nucleosomes_script),
            "-i", input_pairs,
            "-b", bg_pairs,
            "-o", str(nucleosome_root),
            "-x", args.nucleosome_command,
            "-p", args.danpos_py,
            "-t", str(args.nucleosome_threads),
        ]
        if args.skip_nucleosome_index:
            cmd.append("--skip-index")
        extra = shlex_list(args.nucleosome_extra_args)
        if extra:
            cmd.append("--")
            cmd.extend(extra)
        run_cmd(cmd, log_dir / "nucleosome.log")

    if "nucleosome" in steps or (samples["positions_xls"].astype(str).str.strip() == "").any():
        samples = enrich_positions_from_nucleosome(samples, nucleosome_root)
    else:
        needs_dpos = samples["dpos_xls"].astype(str).str.strip() == ""
        samples.loc[needs_dpos, "dpos_xls"] = samples.loc[needs_dpos, "positions_xls"]
        needs_ins = samples["insertion_bw"].astype(str).str.strip() == ""
        samples.loc[needs_ins, "insertion_bw"] = samples.loc[needs_ins, "atac_bw"]

    if (samples["positions_xls"].astype(str).str.strip() == "").any():
        bad = samples.loc[samples["positions_xls"].astype(str).str.strip() == "", "sample"].tolist()
        raise PipelineError(f"positions_xls is still missing after setup for samples: {bad}")

    samples.to_csv(positions_sample_sheet, sep="\t", index=False)

    if "array" in steps:
        ensure_dir(array_root)
        for row in samples.itertuples(index=False):
            sample = row.sample
            sample_out = array_root / sample
            ensure_dir(sample_out)
            cmd = [
                args.python, str(call_array_script),
                "-i", row.positions_xls,
                "-O", str(sample_out),
                "--method", args.array_method,
                "--k", str(args.array_k),
                "--peak_mode", args.peak_mode,
                "--min_overlap_frac", str(args.min_overlap_frac),
                "--chrom_allow_regex", args.chrom_allow_regex,
            ]
            if args.peak_mode == "sicer2":
                if not row.treatment:
                    raise PipelineError(f"Sample {sample} is missing treatment but --peak_mode sicer2 was requested")
                species = row.species or ""
                if not species:
                    raise PipelineError(f"Sample {sample} is missing species but --peak_mode sicer2 was requested")
                cmd.extend(["--treatment", row.treatment, "--species", species])
                if row.control:
                    cmd.extend(["--control", row.control])
            else:
                peak_bed = row.peak_bed or ""
                if not peak_bed:
                    raise PipelineError(f"Sample {sample} is missing peak_bed but --peak_mode existing was requested")
                cmd.extend(["--peak_bed", peak_bed])
            cmd.extend(shlex_list(args.array_extra_args))
            run_cmd(cmd, log_dir / f"array.{sample}.log")

    cluster_df = build_cluster_sample_sheet(samples, array_root, cluster_sample_sheet)

    region_stats_df = cluster_df[["cell_type", "stats_tsv"]].rename(columns={"cell_type": "sample"})
    region_stats_df.to_csv(region_stats_sheet, sep="\t", index=False)

    peak_df = cluster_df[["cell_type", "peak_bed"]].rename(columns={"cell_type": "sample"})
    peak_df.to_csv(peak_sheet, sep="\t", index=False)

    if "cluster" in steps:
        ensure_dir(cluster_root)
        cmd = [
            args.rscript, str(cluster_script),
            "--sample_sheet", str(cluster_sample_sheet),
            "--out_dir", str(cluster_root),
            "--ks", args.ks,
            "--best_k_source", args.best_k_source,
        ]
        if args.sample_order:
            cmd.extend(["--sample_order", args.sample_order])
        cmd.extend(shlex_list(args.cluster_extra_args))
        run_cmd(cmd, log_dir / "cluster.log")

    best_k, cluster_tsv = locate_best_k_cluster_tsv(cluster_root)

    tss_anno_tsv = args.tss_anno_tsv
    if "describe" in steps:
        if not resolved_gtf:
            raise PipelineError(
                "--gtf is required when running the describe step, unless a single shared annotation_gtf is present in sample_sheet"
            )
        ensure_dir(describe_root)
        cmd = [
            args.rscript, str(describe_script),
            "--cluster_tsv", cluster_tsv,
            "--out_dir", str(describe_root),
            "--region_stats_sheet", str(region_stats_sheet),
            "--gtf", resolved_gtf,
            "--use_transcript_tss", bool_str(args.describe_use_transcript_tss),
        ]
        if args.sample_order:
            cmd.extend(["--sample_order", args.sample_order])
        if resolved_ref_fa:
            cmd.extend(["--ref_fa", resolved_ref_fa])
        if peak_sheet.exists():
            cmd.extend(["--peak_sheet", str(peak_sheet)])
        cmd.extend(shlex_list(args.describe_extra_args))
        run_cmd(cmd, log_dir / "describe.log")
        tss_anno_tsv = infer_tss_output(describe_root, use_transcript_tss=args.describe_use_transcript_tss)

    if "motif" in steps:
        if not resolved_motif_annotation:
            raise PipelineError(
                "--motif_annotation_file is required when running the motif step, unless it can be inherited "
                "from sample_sheet.motif_annotation_file or sample_sheet.annotation_gtf"
            )
        if not resolved_chrom_sizes:
            raise PipelineError(
                "--chrom_sizes is required when running the motif step, unless a single shared chrom_sizes is present in sample_sheet"
            )
        if not resolved_genome:
            raise PipelineError(
                "--genome is required when running the motif step, unless a single shared genome or species value is present in sample_sheet"
            )
        ensure_dir(motif_root)
        cmd = [
            args.python, str(motif_script),
            "--cluster_run_dir", str(cluster_root),
            "-O", str(motif_root),
            "--annotation_file", resolved_motif_annotation,
            "--chrom_sizes", resolved_chrom_sizes,
            "--genome", resolved_genome,
            "--mode", args.motif_mode,
            "--target_clusters", args.target_clusters,
            "--tss_near_bp", str(args.tss_near_bp),
        ]

        if args.homer_auto_bg:
            cmd.append("--homer_auto_bg")
        if args.build_only:
            cmd.append("--build_only")

        if args.motif_mode == "internal_linker":
            dmap = {row["cell_type"]: row["dpos_xls"] for _, row in cluster_df.iterrows()}
            dpos_map_str = ";".join(f"{k}={v}" for k, v in dmap.items())
            cmd.extend(["--dpos_map", dpos_map_str])
        else:
            imap = {row["cell_type"]: row["insertion_bw"] for _, row in cluster_df.iterrows()}
            ins_bw_map_str = ";".join(f"{k}={v}" for k, v in imap.items())
            cmd.extend([
                "--ins_bw_map", ins_bw_map_str,
                "--cluster_dir_default", args.cluster_dir_default,
                "--rel_range", args.rel_range,
            ])
            if args.cluster_dir_map:
                cmd.extend(["--cluster_dir_map", args.cluster_dir_map])

        cmd.extend(shlex_list(args.motif_extra_args))
        run_cmd(cmd, log_dir / "motif.log")

    if "score" in steps:
        if not args.expr_tsv:
            raise PipelineError("--expr_tsv is required when running the score step")
        if not tss_anno_tsv:
            raise PipelineError(
                "Scoring requires TSS annotation. Provide --tss_anno_tsv, or run the describe step with --gtf."
            )
        ensure_dir(score_root)
        cmd = [
            args.rscript, str(score_script),
            "--cluster_run_dir", str(cluster_root),
            "--motif_run_dir", str(motif_root),
            "--tss_anno_tsv", str(tss_anno_tsv),
            "--expr_tsv", args.expr_tsv,
            "--out_dir", str(score_root),
            "--promoter_flag_col", args.promoter_flag_col,
            "--restrict_clusters", args.restrict_clusters,
            "--gene_id_candidates", args.gene_id_candidates,
            "--motif_q_cutoff", str(args.motif_q_cutoff),
            "--min_group_n", str(args.min_group_n),
            "--top_n_plot", str(args.top_n_plot),
        ]
        if args.sample_order:
            cmd.extend(["--sample_order", args.sample_order])
        if args.expr_id_col:
            cmd.extend(["--expr_id_col", args.expr_id_col])
        if args.expr_sample_regex_map:
            cmd.extend(["--expr_sample_regex_map", args.expr_sample_regex_map])
        cmd.extend(shlex_list(args.score_extra_args))
        run_cmd(cmd, log_dir / "score.log")

    metadata_rows = [
        ("sample_sheet", str(sample_sheet)),
        ("sample_sheet_with_positions", str(positions_sample_sheet)),
        ("out_dir", str(out_dir)),
        ("steps", ",".join(steps)),
        ("nucleosome_root", str(nucleosome_root)),
        ("array_root", str(array_root)),
        ("cluster_root", str(cluster_root)),
        ("describe_root", str(describe_root)),
        ("motif_root", str(motif_root)),
        ("score_root", str(score_root)),
        ("cluster_sample_sheet", str(cluster_sample_sheet)),
        ("best_k", str(best_k)),
        ("cluster_tsv", str(cluster_tsv)),
        ("resolved_gtf", resolved_gtf),
        ("resolved_ref_fa", resolved_ref_fa),
        ("resolved_motif_annotation", resolved_motif_annotation),
        ("resolved_chrom_sizes", resolved_chrom_sizes),
        ("resolved_genome", resolved_genome),
        ("tss_anno_tsv", str(tss_anno_tsv) if tss_anno_tsv else ""),
    ]
    write_pipeline_metadata(out_dir / "pipeline_run_metadata.txt", metadata_rows)
    eprint(f"[INFO] Pipeline finished. Metadata: {out_dir / 'pipeline_run_metadata.txt'}")


if __name__ == "__main__":
    try:
        main()
    except PipelineError as e:
        eprint(f"[ERROR] {e}")
        sys.exit(1)
