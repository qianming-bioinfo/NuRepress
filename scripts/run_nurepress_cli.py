#!/usr/bin/env python3
import argparse
import glob
import shlex
import subprocess
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


class PipelineError(RuntimeError):
    pass


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def split_csv(x: Optional[str]) -> List[str]:
    if x is None or str(x).strip() == "":
        return []
    return [i.strip() for i in str(x).split(",") if i.strip()]


def shlex_list(x: Optional[str]) -> List[str]:
    if x is None or str(x).strip() == "":
        return []
    return shlex.split(str(x))


def parse_sample_map(x: Optional[str], arg_name: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if x is None or str(x).strip() == "":
        return out

    for item in str(x).split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise PipelineError(
                f"Invalid format for {arg_name}: '{item}'. "
                f"Expected: sample1=/path/a;sample2=/path/b"
            )
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()
        if not k:
            raise PipelineError(f"Empty sample name in {arg_name}")
        if not v:
            raise PipelineError(f"Empty value for sample '{k}' in {arg_name}")
        if k in out:
            raise PipelineError(f"Duplicated sample '{k}' in {arg_name}")
        out[k] = v
    return out


def materialize_sample_field(
    sample_list: List[str],
    mapping: Dict[str, str],
    arg_name: str,
    required: bool = False,
    default: str = "",
) -> List[str]:
    unknown = sorted(set(mapping.keys()) - set(sample_list))
    if unknown:
        raise PipelineError(
            f"{arg_name} contains unknown samples not listed in --samples: {unknown}"
        )

    out = {s: default for s in sample_list}
    out.update(mapping)

    if required:
        missing = [s for s in sample_list if str(out[s]).strip() == ""]
        if missing:
            raise PipelineError(
                f"{arg_name} is required for all samples. Missing in: {missing}"
            )

    return [str(out[s]).strip() for s in sample_list]


def repo_script_path(base_dir, filename):
    p = (Path(base_dir) / filename).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Bundled script not found: {p}")
    return p


def run_cmd(cmd: List[str], log_path: Path, cwd: Optional[Path] = None):
    ensure_dir(log_path.parent)
    eprint("[RUN]", " ".join(shlex.quote(x) for x in cmd))
    with open(log_path, "w") as log:
        log.write("COMMAND=\n")
        log.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
        proc = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=str(cwd) if cwd else None,
        )
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


def resolve_shared_value(
    cli_value: Optional[str],
    samples: pd.DataFrame,
    column: str,
    fallback_columns: Optional[List[str]] = None,
) -> str:
    if cli_value and str(cli_value).strip():
        return str(cli_value).strip()
    tried = [column] + (fallback_columns or [])
    for col in tried:
        val = first_unique_nonempty(samples, col)
        if val:
            return val
    return ""


def standard_paths(out_dir: Path, nucleosome_out_subdir: str) -> Dict[str, Path]:
    intermediate_dir = out_dir / "intermediate"
    return {
        "out_dir": out_dir,
        "log_dir": out_dir / "logs",
        "nucleosome_root": out_dir / nucleosome_out_subdir,
        "array_root": out_dir / "01_array_call",
        "cluster_root": out_dir / "02_cluster",
        "describe_root": out_dir / "03_describe",
        "motif_root": out_dir / "04_motif",
        "score_root": out_dir / "05_score",
        "intermediate_dir": intermediate_dir,
        "positions_sample_sheet": intermediate_dir / "sample_sheet.with_positions.tsv",
        "cluster_sample_sheet": intermediate_dir / "cluster_sample_sheet.tsv",
        "region_stats_sheet": intermediate_dir / "region_stats_sheet.tsv",
        "peak_sheet": intermediate_dir / "peak_sheet.tsv",
        "pipeline_metadata": out_dir / "pipeline_run_metadata.txt",
    }


def ensure_standard_dirs(paths: Dict[str, Path]):
    for k in ["out_dir", "log_dir", "intermediate_dir"]:
        ensure_dir(paths[k])


def build_samples_df(
    sample_list: List[str],
    atac_bw_map: Optional[str] = "",
    positions_map: Optional[str] = "",
    treatment_map: Optional[str] = "",
    control_map: Optional[str] = "",
    peak_bed_map: Optional[str] = "",
    dpos_map: Optional[str] = "",
    insertion_bw_map: Optional[str] = "",
    genome: Optional[str] = "",
    genome_fasta: Optional[str] = "",
    chrom_sizes: Optional[str] = "",
    annotation_gtf: Optional[str] = "",
    require_atac_bw: bool = False,
    require_positions: bool = False,
    require_treatment: bool = False,
    require_peak_bed: bool = False,
    require_genome: bool = False,
) -> pd.DataFrame:
    if not sample_list:
        raise PipelineError("--samples cannot be empty")

    if len(set(sample_list)) != len(sample_list):
        dup = [x for x in sample_list if sample_list.count(x) > 1]
        raise PipelineError(f"Duplicated sample names in --samples: {sorted(set(dup))}")

    atac_bw = parse_sample_map(atac_bw_map, "--atac_bw_map")
    positions = parse_sample_map(positions_map, "--positions_map")
    treatment = parse_sample_map(treatment_map, "--treatment_map")
    control = parse_sample_map(control_map, "--control_map")
    peak_bed = parse_sample_map(peak_bed_map, "--peak_bed_map")
    dpos = parse_sample_map(dpos_map, "--dpos_map")
    insertion_bw = parse_sample_map(insertion_bw_map, "--insertion_bw_map")

    genome = str(genome or "").strip()
    genome_fasta = str(genome_fasta or "").strip()
    chrom_sizes = str(chrom_sizes or "").strip()
    annotation_gtf = str(annotation_gtf or "").strip()

    if require_genome and not genome:
        raise PipelineError("--genome is required for the selected command")

    df = pd.DataFrame({
        "sample": sample_list,
        "atac_bw": materialize_sample_field(
            sample_list, atac_bw, "--atac_bw_map", required=require_atac_bw
        ),
        "positions_xls": materialize_sample_field(
            sample_list, positions, "--positions_map", required=require_positions
        ),
        "treatment": materialize_sample_field(
            sample_list, treatment, "--treatment_map", required=require_treatment
        ),
        "control": materialize_sample_field(
            sample_list, control, "--control_map", required=False
        ),
        "peak_bed": materialize_sample_field(
            sample_list, peak_bed, "--peak_bed_map", required=require_peak_bed
        ),
        "genome": [genome] * len(sample_list),
        "species": [genome] * len(sample_list),
        "dpos_xls": materialize_sample_field(
            sample_list, dpos, "--dpos_map", required=False
        ),
        "insertion_bw": materialize_sample_field(
            sample_list, insertion_bw, "--insertion_bw_map", required=False
        ),
        "genome_fasta": [genome_fasta] * len(sample_list),
        "chrom_sizes": [chrom_sizes] * len(sample_list),
        "annotation_gtf": [annotation_gtf] * len(sample_list),
    })
    return df


def _normalize_token(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", s).lower()

def discover_danpos_position_file(nucleosome_out_dir: Path, sample: str) -> str:
    pooled = nucleosome_out_dir / "pooled"
    if not pooled.exists():
        raise PipelineError(f"DANPOS pooled directory not found: {pooled}")

    sample_norm = _normalize_token(sample)

    ref_adjust_hits: List[Path] = []
    plain_hits: List[Path] = []

    for fp in pooled.iterdir():
        if not fp.is_file():
            continue

        name = fp.name
        name_norm = _normalize_token(name)

        if sample_norm not in name_norm:
            continue

        if name.endswith(".positions.ref_adjust.xls"):
            ref_adjust_hits.append(fp)
        elif name.endswith(".positions.xls"):
            plain_hits.append(fp)

    hits = sorted(ref_adjust_hits) + sorted(plain_hits)

    if not hits:
        raise PipelineError(
            f"Cannot find DANPOS position output for sample '{sample}' under: {pooled}. "
            f"Available position files: "
            f"{', '.join(sorted(p.name for p in pooled.iterdir() if p.is_file() and 'positions' in p.name))}"
        )

    return str(hits[0])


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
            "sample": sample,
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
        raise PipelineError(f"Missing cluster_tables under: {cluster_dir}")

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


def write_positions_sheet(samples: pd.DataFrame, path: Path):
    ensure_dir(path.parent)
    samples.to_csv(path, sep="\t", index=False)


def require_file(path: Path, msg: str):
    if not path.exists():
        raise PipelineError(msg)


def run_nucleosome_step(args, samples: pd.DataFrame, paths: Dict[str, Path]):
    missing_treat = samples["treatment"].astype(str).str.strip() == ""
    if missing_treat.any():
        bad = samples.loc[missing_treat, "sample"].tolist()
        raise PipelineError(
            f"The nucleosome step requires treatment for all samples. Missing in: {bad}"
        )

    control_series = samples["control"].astype(str).str.strip()
    has_control = control_series != ""

    if has_control.any() and not has_control.all():
        bad = samples.loc[~has_control, "sample"].tolist()
        raise PipelineError(
            "For the nucleosome step, control must be provided for either all samples or none. "
            f"Missing control in: {bad}"
        )

    input_pairs = ",".join(f"{r.sample}:{r.treatment}" for r in samples.itertuples(index=False))

    ensure_dir(paths["nucleosome_root"])
    cmd = [
        args.bash,
        str(args.call_nucleosomes_script),
        "-i",
        input_pairs,
        "-o",
        str(paths["nucleosome_root"]),
        "-x",
        args.nucleosome_command,
        "-p",
        args.danpos_py,
        "-t",
        str(args.nucleosome_threads),
    ]

    if has_control.all():
        bg_pairs = ",".join(f"{r.sample}:{r.control}" for r in samples.itertuples(index=False))
        cmd.extend(["-b", bg_pairs])
        eprint("[INFO] Nucleosome step will run with background/control.")
    else:
        eprint("[INFO] Nucleosome step will run without background/control.")

    if args.skip_nucleosome_index:
        cmd.append("--skip-index")

    extra = shlex_list(args.nucleosome_extra_args)
    if extra:
        cmd.append("--")
        cmd.extend(extra)

    run_cmd(cmd, paths["log_dir"] / "nucleosome.log")


def run_array_step(args, samples: pd.DataFrame, paths: Dict[str, Path], infer_from_nucleosome: bool = True) -> pd.DataFrame:
    df = samples.copy()

    if infer_from_nucleosome and (df["positions_xls"].astype(str).str.strip() == "").any():
        df = enrich_positions_from_nucleosome(df, paths["nucleosome_root"])

    if (df["positions_xls"].astype(str).str.strip() == "").any():
        bad = df.loc[df["positions_xls"].astype(str).str.strip() == "", "sample"].tolist()
        raise PipelineError(f"positions_xls is missing for samples: {bad}")

    needs_dpos = df["dpos_xls"].astype(str).str.strip() == ""
    df.loc[needs_dpos, "dpos_xls"] = df.loc[needs_dpos, "positions_xls"]

    needs_ins = df["insertion_bw"].astype(str).str.strip() == ""
    df.loc[needs_ins, "insertion_bw"] = df.loc[needs_ins, "atac_bw"]

    write_positions_sheet(df, paths["positions_sample_sheet"])

    ensure_dir(paths["array_root"])
    for row in df.itertuples(index=False):
        sample = row.sample
        sample_out = paths["array_root"] / sample
        ensure_dir(sample_out)

        cmd = [
            args.python,
            str(args.call_array_script),
            "-i",
            row.positions_xls,
            "-O",
            str(sample_out),
            "--method",
            args.array_method,
            "--k",
            str(args.array_k),
            "--peak_mode",
            args.peak_mode,
            "--min_overlap_frac",
            str(args.min_overlap_frac),
            "--chrom_allow_regex",
            args.chrom_allow_regex,
        ]

        if args.peak_mode == "sicer2":
            if not row.treatment:
                raise PipelineError(f"Sample {sample} is missing treatment but --peak_mode sicer2 was requested")
            genome = row.genome or ""
            if not genome:
                raise PipelineError(f"Sample {sample} is missing genome but --peak_mode sicer2 was requested")
            cmd.extend(["--treatment", row.treatment, "--species", genome])
            if row.control:
                cmd.extend(["--control", row.control])
        else:
            peak_bed = row.peak_bed or ""
            if not peak_bed:
                raise PipelineError(f"Sample {sample} is missing peak_bed but --peak_mode existing was requested")
            cmd.extend(["--peak_bed", peak_bed])

        cmd.extend(shlex_list(args.array_extra_args))
        run_cmd(cmd, paths["log_dir"] / f"array.{sample}.log")

    return df


def prepare_cluster_bundle(samples: pd.DataFrame, paths: Dict[str, Path]) -> pd.DataFrame:
    cluster_df = build_cluster_sample_sheet(samples, paths["array_root"], paths["cluster_sample_sheet"])

    region_stats_df = cluster_df[["sample", "stats_tsv"]].copy()
    region_stats_df.to_csv(paths["region_stats_sheet"], sep="\t", index=False)

    peak_df = cluster_df[["sample", "peak_bed"]].copy()
    peak_df.to_csv(paths["peak_sheet"], sep="\t", index=False)

    return cluster_df


def run_cluster_step(args, paths: Dict[str, Path]):
    require_file(
        paths["cluster_sample_sheet"],
        "Missing cluster_sample_sheet.tsv. Please run the array/cluster preparation first under the same --out_dir.",
    )
    ensure_dir(paths["cluster_root"])
    cmd = [
        args.rscript,
        str(args.cluster_script),
        "--sample_sheet",
        str(paths["cluster_sample_sheet"]),
        "--out_dir",
        str(paths["cluster_root"]),
        "--ks",
        args.ks,
        "--best_k_source",
        args.best_k_source,
    ]
    if getattr(args, "sample_order", None):
        cmd.extend(["--sample_order", args.sample_order])
    cmd.extend(shlex_list(args.cluster_extra_args))
    run_cmd(cmd, paths["log_dir"] / "cluster.log")


def run_describe_step(args, paths: Dict[str, Path], resolved_gtf: str, resolved_ref_fa: str) -> str:
    require_file(
        paths["region_stats_sheet"],
        "Missing region_stats_sheet.tsv. Please run cluster subcommand first under the same --out_dir.",
    )

    best_k, cluster_tsv = locate_best_k_cluster_tsv(paths["cluster_root"])
    _ = best_k

    ensure_dir(paths["describe_root"])
    cmd = [
        args.rscript,
        str(args.describe_script),
        "--cluster_tsv",
        cluster_tsv,
        "--out_dir",
        str(paths["describe_root"]),
        "--region_stats_sheet",
        str(paths["region_stats_sheet"]),
        "--gtf",
        resolved_gtf,
        "--use_transcript_tss",
        bool_str(args.describe_use_transcript_tss),
    ]
    if getattr(args, "sample_order", None):
        cmd.extend(["--sample_order", args.sample_order])
    if resolved_ref_fa:
        cmd.extend(["--ref_fa", resolved_ref_fa])
    if paths["peak_sheet"].exists():
        cmd.extend(["--peak_sheet", str(paths["peak_sheet"])])
    cmd.extend(shlex_list(args.describe_extra_args))
    run_cmd(cmd, paths["log_dir"] / "describe.log")

    return infer_tss_output(
        paths["describe_root"],
        use_transcript_tss=args.describe_use_transcript_tss,
    )


def run_motif_step(args, paths: Dict[str, Path], resolved_annotation_gtf: str, resolved_chrom_sizes: str, resolved_genome: str):
    require_file(
        paths["cluster_sample_sheet"],
        "Missing cluster_sample_sheet.tsv. Please run cluster subcommand first under the same --out_dir.",
    )
    cluster_df = pd.read_csv(paths["cluster_sample_sheet"], sep="\t", dtype=str).fillna("")

    ensure_dir(paths["motif_root"])
    cmd = [
        args.python,
        str(args.motif_script),
        "--cluster_run_dir",
        str(paths["cluster_root"]),
        "-O",
        str(paths["motif_root"]),
        "--annotation_file",
        resolved_annotation_gtf,
        "--chrom_sizes",
        resolved_chrom_sizes,
        "--genome",
        resolved_genome,
        "--mode",
        args.motif_mode,
        "--target_clusters",
        args.target_clusters,
        "--tss_near_bp",
        str(args.tss_near_bp),
    ]

    if args.homer_auto_bg:
        cmd.append("--homer_auto_bg")
    if args.build_only:
        cmd.append("--build_only")

    if args.motif_mode == "internal_linker":
        dmap = {row["sample"]: row["dpos_xls"] for _, row in cluster_df.iterrows()}
        dpos_map_str = ";".join(f"{k}={v}" for k, v in dmap.items())
        cmd.extend(["--dpos_map", dpos_map_str])
    else:
        imap = {row["sample"]: row["insertion_bw"] for _, row in cluster_df.iterrows()}
        ins_bw_map_str = ";".join(f"{k}={v}" for k, v in imap.items())
        cmd.extend([
            "--ins_bw_map",
            ins_bw_map_str,
            "--cluster_dir_default",
            args.cluster_dir_default,
            "--rel_range",
            args.rel_range,
        ])
        if args.cluster_dir_map:
            cmd.extend(["--cluster_dir_map", args.cluster_dir_map])

    cmd.extend(shlex_list(args.motif_extra_args))
    run_cmd(cmd, paths["log_dir"] / "motif.log")


def run_score_step(args, paths: Dict[str, Path], tss_anno_tsv: str):
    ensure_dir(paths["score_root"])
    cmd = [
        args.rscript,
        str(args.score_script),
        "--cluster_run_dir",
        str(paths["cluster_root"]),
        "--motif_run_dir",
        str(paths["motif_root"]),
        "--tss_anno_tsv",
        str(tss_anno_tsv),
        "--expr_tsv",
        args.expr_tsv,
        "--out_dir",
        str(paths["score_root"]),
        "--promoter_flag_col",
        args.promoter_flag_col,
        "--restrict_clusters",
        args.restrict_clusters,
        "--gene_id_candidates",
        args.gene_id_candidates,
        "--motif_q_cutoff",
        str(args.motif_q_cutoff),
        "--min_group_n",
        str(args.min_group_n),
        "--top_n_plot",
        str(args.top_n_plot),
    ]
    if getattr(args, "sample_order", None):
        cmd.extend(["--sample_order", args.sample_order])
    if args.expr_id_col:
        cmd.extend(["--expr_id_col", args.expr_id_col])
    if args.expr_sample_regex_map:
        cmd.extend(["--expr_sample_regex_map", args.expr_sample_regex_map])
    cmd.extend(shlex_list(args.score_extra_args))
    run_cmd(cmd, paths["log_dir"] / "score.log")


def add_common_tool_args(p):
    p.add_argument("--python", default=sys.executable, help="Python executable")
    p.add_argument("--rscript", default="Rscript", help="Rscript executable")
    p.add_argument("--bash", default="bash", help="Bash executable")
    p.add_argument("--call_nucleosomes_script", default=None, help="Path to call_nucleosomes.sh")
    p.add_argument("--call_array_script", default=None, help="Path to call_peak_and_identify_array.py")
    p.add_argument("--cluster_script", default=None, help="Path to cluster_array_subtype.R")
    p.add_argument("--describe_script", default=None, help="Path to describe_array_subtype.R")
    p.add_argument("--motif_script", default=None, help="Path to run_array_motif_enrichment.py")
    p.add_argument("--score_script", default=None, help="Path to score_repressor_candidates.R")


def add_common_io_args(p):
    p.add_argument("--out_dir", required=True, help="Pipeline output directory")
    p.add_argument("--nucleosome_out_subdir", default="00_nucleosome_call", help="Subdirectory for DANPOS wrapper output")


def add_common_sample_args(p):
    p.add_argument("--samples", required=True, help="Comma-separated sample names, e.g. sample1,sample2,sample3")
    p.add_argument("--atac_bw_map", default="", help="sample=/path;sample=/path")
    p.add_argument("--positions_map", default="", help="sample=/path/to/*.positions.xls;...")
    p.add_argument("--treatment_map", default="", help="sample=/path/to/treatment.bam;...")
    p.add_argument("--control_map", default="", help="sample=/path/to/control.bam;... optional")
    p.add_argument("--peak_bed_map", default="", help="sample=/path/to/peak.bed;... used when --peak_mode existing")
    p.add_argument("--dpos_map", default="", help="sample=/path/to/dpos_or_positions.xls;... optional")
    p.add_argument("--insertion_bw_map", default="", help="sample=/path/to/insertion.bw;... optional")


def add_common_reference_args(p):
    p.add_argument("--genome", default=None, help="Shared genome for all samples, e.g. hg38 or mm10")
    p.add_argument("--genome_fasta", default=None, help="Shared reference FASTA")
    p.add_argument("--chrom_sizes", default=None, help="Shared chrom sizes file")
    p.add_argument("--annotation_gtf", default=None, help="Shared annotation file for both describe and motif steps")


def add_sample_order_arg(p):
    p.add_argument("--sample_order", default=None, help="Comma-separated sample order")


def add_nucleosome_args(p):
    p.add_argument("--nucleosome_command", default="dpos", help="DANPOS subcommand passed to call_nucleosomes.sh")
    p.add_argument("--nucleosome_threads", type=int, default=1, help="Threads passed to call_nucleosomes.sh")
    p.add_argument("--danpos_py", default="danpos.py", help="Path to danpos.py passed to call_nucleosomes.sh")
    p.add_argument("--skip_nucleosome_index", action="store_true", default=False, help="Pass --skip-index to call_nucleosomes.sh")
    p.add_argument("--nucleosome_extra_args", default="", help="Extra arguments passed to call_nucleosomes.sh after --")


def add_array_args(p):
    p.add_argument("--array_method", choices=["seedextend", "mergewindow"], default="seedextend")
    p.add_argument("--array_k", type=int, default=3)
    p.add_argument("--peak_mode", choices=["sicer2", "existing"], default="sicer2")
    p.add_argument("--min_overlap_frac", type=float, default=0.30)
    p.add_argument("--chrom_allow_regex", default=r"^chr([0-9]+|X|Y)$")
    p.add_argument("--array_extra_args", default="", help="Extra arguments passed to call_peak_and_identify_array.py")


def add_cluster_args(p):
    p.add_argument("--ks", default="2,3,4,5,6")
    p.add_argument("--best_k_source", default="within_sample")
    p.add_argument("--cluster_extra_args", default="", help="Extra arguments passed to cluster_array_subtype.R")


def add_describe_args(p):
    p.add_argument("--describe_use_transcript_tss", action="store_true", default=True)
    p.add_argument("--describe_gene_tss", dest="describe_use_transcript_tss", action="store_false")
    p.add_argument("--describe_extra_args", default="", help="Extra arguments passed to describe_array_subtype.R")


def add_motif_args(p):
    p.add_argument("--motif_mode", choices=["internal_linker", "edge_outside"], default="internal_linker")
    p.add_argument("--target_clusters", default="ALL")
    p.add_argument("--tss_near_bp", type=int, default=3000)
    p.add_argument("--cluster_dir_default", choices=["neutral", "inside_low", "inside_high"], default="neutral")
    p.add_argument("--cluster_dir_map", default=None, help="Example: C1=inside_low;C2=neutral")
    p.add_argument("--rel_range", default="+275:+475")
    p.add_argument("--homer_auto_bg", action="store_true", default=False)
    p.add_argument("--build_only", action="store_true", default=False)
    p.add_argument("--motif_extra_args", default="", help="Extra arguments passed to run_array_motif_enrichment.py")


def add_score_args(p):
    p.add_argument("--expr_tsv", default=None, help="Expression matrix TSV for scoring step")
    p.add_argument("--tss_anno_tsv", default=None, help="Optional existing TSS annotation TSV. Overrides describe-step output when provided")
    p.add_argument("--promoter_flag_col", default="overlap_TSSwin_2000")
    p.add_argument("--restrict_clusters", default="ALL")
    p.add_argument("--expr_id_col", default=None)
    p.add_argument("--expr_sample_regex_map", default=None)
    p.add_argument("--gene_id_candidates", default="nearest_TSS_gene_name,gene_symbol_mapped,nearest_TSS_id")
    p.add_argument("--motif_q_cutoff", type=float, default=0.05)
    p.add_argument("--min_group_n", type=int, default=5)
    p.add_argument("--top_n_plot", type=int, default=20)
    p.add_argument("--score_extra_args", default="", help="Extra arguments passed to score_repressor_candidates.R")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "NuRepress pipeline runner with two levels: "
            "a top-level 'run' controller and per-step subcommands."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    sub = ap.add_subparsers(dest="subcommand", required=True)

    p_run = sub.add_parser("run", help="Run multiple steps under one controller")
    add_common_io_args(p_run)
    add_common_tool_args(p_run)
    add_common_sample_args(p_run)
    add_common_reference_args(p_run)
    add_nucleosome_args(p_run)
    add_array_args(p_run)
    add_cluster_args(p_run)
    add_describe_args(p_run)
    add_motif_args(p_run)
    add_score_args(p_run)
    add_sample_order_arg(p_run)
    p_run.add_argument(
        "--steps",
        default="nucleosome,array,cluster,describe,motif,score",
        help="Comma-separated steps to run. Choices: nucleosome,array,cluster,describe,motif,score",
    )

    p_nuc = sub.add_parser("nucleosome", help="Run only DANPOS wrapper step")
    add_common_io_args(p_nuc)
    add_common_tool_args(p_nuc)
    add_common_sample_args(p_nuc)
    add_nucleosome_args(p_nuc)

    p_array = sub.add_parser("array", help="Run only array calling step")
    add_common_io_args(p_array)
    add_common_tool_args(p_array)
    add_common_sample_args(p_array)
    add_common_reference_args(p_array)
    add_array_args(p_array)

    p_cluster = sub.add_parser("cluster", help="Prepare cluster bundle and run cluster step")
    add_common_io_args(p_cluster)
    add_common_tool_args(p_cluster)
    add_common_sample_args(p_cluster)
    add_common_reference_args(p_cluster)
    add_cluster_args(p_cluster)
    add_sample_order_arg(p_cluster)

    p_describe = sub.add_parser("describe", help="Run describe/TSS annotation step")
    add_common_io_args(p_describe)
    add_common_tool_args(p_describe)
    add_common_reference_args(p_describe)
    add_describe_args(p_describe)
    add_sample_order_arg(p_describe)

    p_motif = sub.add_parser("motif", help="Run motif enrichment step")
    add_common_io_args(p_motif)
    add_common_tool_args(p_motif)
    add_common_reference_args(p_motif)
    add_motif_args(p_motif)

    p_score = sub.add_parser("score", help="Run scoring step")
    add_common_io_args(p_score)
    add_common_tool_args(p_score)
    add_score_args(p_score)
    add_sample_order_arg(p_score)

    return ap


def resolve_tool_paths(args):
    script_base = Path(__file__).resolve().parent

    if not getattr(args, "call_nucleosomes_script", None):
        args.call_nucleosomes_script = repo_script_path(script_base / "00_nucleosome", "call_nucleosomes.sh")
    else:
        args.call_nucleosomes_script = Path(args.call_nucleosomes_script).resolve()

    if not getattr(args, "call_array_script", None):
        args.call_array_script = repo_script_path(script_base / "01_array", "call_peak_and_identify_array.py")
    else:
        args.call_array_script = Path(args.call_array_script).resolve()

    if not getattr(args, "cluster_script", None):
        args.cluster_script = repo_script_path(script_base / "02_cluster", "cluster_array_subtype.R")
    else:
        args.cluster_script = Path(args.cluster_script).resolve()

    if not getattr(args, "describe_script", None):
        args.describe_script = repo_script_path(script_base / "02_cluster", "describe_array_subtype.R")
    else:
        args.describe_script = Path(args.describe_script).resolve()

    if not getattr(args, "motif_script", None):
        args.motif_script = repo_script_path(script_base / "03_motif", "run_array_motif_enrichment.py")
    else:
        args.motif_script = Path(args.motif_script).resolve()

    if not getattr(args, "score_script", None):
        args.score_script = repo_script_path(script_base / "04_score", "score_repressor_candidates.R")
    else:
        args.score_script = Path(args.score_script).resolve()


def handle_run(args):
    steps = [s.lower() for s in split_csv(args.steps)]
    allowed_steps = {"nucleosome", "array", "cluster", "describe", "motif", "score"}
    unknown_steps = [s for s in steps if s not in allowed_steps]
    if unknown_steps:
        raise PipelineError(f"Unknown steps: {unknown_steps}")

    sample_list = split_csv(args.samples)
    need_atac = any(s in steps for s in ["array", "cluster"])
    need_treatment = ("nucleosome" in steps) or ("array" in steps and args.peak_mode == "sicer2")
    need_peak_bed = ("array" in steps and args.peak_mode == "existing")
    need_genome = ("array" in steps and args.peak_mode == "sicer2") or ("motif" in steps)

    samples = build_samples_df(
        sample_list=sample_list,
        atac_bw_map=args.atac_bw_map,
        positions_map=args.positions_map,
        treatment_map=args.treatment_map,
        control_map=args.control_map,
        peak_bed_map=args.peak_bed_map,
        dpos_map=args.dpos_map,
        insertion_bw_map=args.insertion_bw_map,
        genome=args.genome,
        genome_fasta=args.genome_fasta,
        chrom_sizes=args.chrom_sizes,
        annotation_gtf=args.annotation_gtf,
        require_atac_bw=need_atac,
        require_positions=False,
        require_treatment=need_treatment,
        require_peak_bed=need_peak_bed,
        require_genome=need_genome,
    )

    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)

    if "nucleosome" in steps:
        run_nucleosome_step(args, samples, paths)

    if "array" in steps:
        samples = run_array_step(args, samples, paths, infer_from_nucleosome=True)

    if "cluster" in steps:
        if (samples["positions_xls"].astype(str).str.strip() == "").any():
            samples = enrich_positions_from_nucleosome(samples, paths["nucleosome_root"])
            write_positions_sheet(samples, paths["positions_sample_sheet"])
        prepare_cluster_bundle(samples, paths)
        run_cluster_step(args, paths)

    tss_anno_tsv = args.tss_anno_tsv
    if "describe" in steps:
        resolved_gtf = resolve_shared_value(args.annotation_gtf, samples, "annotation_gtf")
        resolved_ref_fa = resolve_shared_value(args.genome_fasta, samples, "genome_fasta")
        if not resolved_gtf:
            raise PipelineError("--annotation_gtf is required when running describe")
        tss_anno_tsv = run_describe_step(args, paths, resolved_gtf, resolved_ref_fa)

    if "motif" in steps:
        resolved_annotation_gtf = resolve_shared_value(args.annotation_gtf, samples, "annotation_gtf")
        resolved_chrom_sizes = resolve_shared_value(args.chrom_sizes, samples, "chrom_sizes")
        resolved_genome = resolve_shared_value(args.genome, samples, "genome")
        if not resolved_annotation_gtf:
            raise PipelineError("--annotation_gtf is required when running motif")
        if not resolved_chrom_sizes:
            raise PipelineError("--chrom_sizes is required when running motif")
        if not resolved_genome:
            raise PipelineError("--genome is required when running motif")
        run_motif_step(args, paths, resolved_annotation_gtf, resolved_chrom_sizes, resolved_genome)

    if "score" in steps:
        if not args.expr_tsv:
            raise PipelineError("--expr_tsv is required when running score")
        if not tss_anno_tsv:
            raise PipelineError("Scoring requires TSS annotation. Provide --tss_anno_tsv or run describe first.")
        run_score_step(args, paths, tss_anno_tsv)

    metadata_rows = [
        ("subcommand", "run"),
        ("samples", ",".join(samples["sample"].tolist())),
        ("steps", ",".join(steps)),
        ("out_dir", str(paths["out_dir"])),
        ("nucleosome_root", str(paths["nucleosome_root"])),
        ("array_root", str(paths["array_root"])),
        ("cluster_root", str(paths["cluster_root"])),
        ("describe_root", str(paths["describe_root"])),
        ("motif_root", str(paths["motif_root"])),
        ("score_root", str(paths["score_root"])),
        ("positions_sample_sheet", str(paths["positions_sample_sheet"]) if paths["positions_sample_sheet"].exists() else ""),
        ("cluster_sample_sheet", str(paths["cluster_sample_sheet"]) if paths["cluster_sample_sheet"].exists() else ""),
        ("tss_anno_tsv", str(tss_anno_tsv) if tss_anno_tsv else ""),
    ]
    write_pipeline_metadata(paths["pipeline_metadata"], metadata_rows)
    eprint(f"[INFO] Done. Metadata: {paths['pipeline_metadata']}")


def handle_nucleosome(args):
    sample_list = split_csv(args.samples)
    samples = build_samples_df(
        sample_list=sample_list,
        treatment_map=args.treatment_map,
        control_map=args.control_map,
        require_treatment=True,
    )
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)
    run_nucleosome_step(args, samples, paths)


def handle_array(args):
    sample_list = split_csv(args.samples)
    samples = build_samples_df(
        sample_list=sample_list,
        atac_bw_map=args.atac_bw_map,
        positions_map=args.positions_map,
        treatment_map=args.treatment_map,
        control_map=args.control_map,
        peak_bed_map=args.peak_bed_map,
        dpos_map=args.dpos_map,
        insertion_bw_map=args.insertion_bw_map,
        genome=args.genome,
        require_atac_bw=False,
        require_positions=False,
        require_treatment=(args.peak_mode == "sicer2"),
        require_peak_bed=(args.peak_mode == "existing"),
        require_genome=(args.peak_mode == "sicer2"),
    )
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)
    run_array_step(args, samples, paths, infer_from_nucleosome=True)


def handle_cluster(args):
    sample_list = split_csv(args.samples)
    samples = build_samples_df(
        sample_list=sample_list,
        atac_bw_map=args.atac_bw_map,
        positions_map=args.positions_map,
        dpos_map=args.dpos_map,
        insertion_bw_map=args.insertion_bw_map,
        require_atac_bw=True,
        require_positions=False,
    )
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)

    if (samples["positions_xls"].astype(str).str.strip() == "").any():
        samples = enrich_positions_from_nucleosome(samples, paths["nucleosome_root"])
    needs_dpos = samples["dpos_xls"].astype(str).str.strip() == ""
    samples.loc[needs_dpos, "dpos_xls"] = samples.loc[needs_dpos, "positions_xls"]
    needs_ins = samples["insertion_bw"].astype(str).str.strip() == ""
    samples.loc[needs_ins, "insertion_bw"] = samples.loc[needs_ins, "atac_bw"]

    write_positions_sheet(samples, paths["positions_sample_sheet"])
    prepare_cluster_bundle(samples, paths)
    run_cluster_step(args, paths)


def handle_describe(args):
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)

    if not args.annotation_gtf:
        raise PipelineError("--annotation_gtf is required for describe subcommand")

    tss_anno_tsv = run_describe_step(
        args,
        paths,
        resolved_gtf=str(args.annotation_gtf),
        resolved_ref_fa=str(args.genome_fasta or ""),
    )
    eprint(f"[INFO] TSS annotation: {tss_anno_tsv}")


def handle_motif(args):
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)

    if not args.annotation_gtf:
        raise PipelineError("--annotation_gtf is required for motif subcommand")
    if not args.chrom_sizes:
        raise PipelineError("--chrom_sizes is required for motif subcommand")
    if not args.genome:
        raise PipelineError("--genome is required for motif subcommand")

    run_motif_step(
        args,
        paths,
        resolved_annotation_gtf=str(args.annotation_gtf),
        resolved_chrom_sizes=str(args.chrom_sizes),
        resolved_genome=str(args.genome),
    )


def handle_score(args):
    out_dir = Path(args.out_dir).resolve()
    paths = standard_paths(out_dir, args.nucleosome_out_subdir)
    ensure_standard_dirs(paths)
    resolve_tool_paths(args)

    if not args.expr_tsv:
        raise PipelineError("--expr_tsv is required for score subcommand")

    tss_anno_tsv = args.tss_anno_tsv
    if not tss_anno_tsv:
        inferred = paths["describe_root"] / "03_TSS_annotation" / "region_TSS_association.detail.transcriptTSS.tsv"
        if inferred.exists():
            tss_anno_tsv = str(inferred)
        else:
            raise PipelineError(
                "--tss_anno_tsv is required for score subcommand unless describe output already exists"
            )

    run_score_step(args, paths, tss_anno_tsv)


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.subcommand == "run":
        handle_run(args)
    elif args.subcommand == "nucleosome":
        handle_nucleosome(args)
    elif args.subcommand == "array":
        handle_array(args)
    elif args.subcommand == "cluster":
        handle_cluster(args)
    elif args.subcommand == "describe":
        handle_describe(args)
    elif args.subcommand == "motif":
        handle_motif(args)
    elif args.subcommand == "score":
        handle_score(args)
    else:
        raise PipelineError(f"Unsupported subcommand: {args.subcommand}")


if __name__ == "__main__":
    try:
        main()
    except PipelineError as e:
        eprint(f"[ERROR] {e}")
        sys.exit(1)