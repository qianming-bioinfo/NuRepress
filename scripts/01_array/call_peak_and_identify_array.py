#!/usr/bin/env python3
import os
import sys
import re
import math
import glob
import shlex
import shutil
import argparse
import subprocess
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def read_ref_adjust_table(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
        if df.shape[1] >= 2:
            return df
    except Exception:
        pass

    try:
        df = pd.read_csv(path, sep=r"\s+", comment="#", low_memory=False, engine="python")
        if df.shape[1] >= 2:
            return df
    except Exception:
        pass

    try:
        df = pd.read_excel(path)
        return df
    except Exception as e:
        raise RuntimeError(f"Failed to read input file as TSV/whitespace/Excel: {path}\nLast error: {e}")


def basename_before_first_dot(path: str) -> str:
    return os.path.basename(path).split(".", 1)[0]


def fmt_num(x) -> str:
    if x is None:
        return "NA"
    if isinstance(x, (int, np.integer)):
        s = str(int(x))
    else:
        s = f"{float(x):g}"
    return s.replace(".", "p")


def clean_positions_df(df: pd.DataFrame, chrom_allow_regex: str, keep_all_centers: bool) -> pd.DataFrame:
    required = ["chr", "start", "end", "smt_pos", "smt_value", "fuzziness_score"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}\nAvailable columns: {list(df.columns)}")

    out = df[required].copy()
    out = out.rename(columns={
        "smt_pos": "center",
        "smt_value": "signal",
        "fuzziness_score": "fuzz",
    })

    for c in ["start", "end", "center", "signal", "fuzz"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    out = out.replace([np.inf, -np.inf], np.nan).dropna().copy()
    out["start"] = out["start"].astype(np.int64)
    out["end"] = out["end"].astype(np.int64)
    out["center"] = out["center"].astype(np.int64)
    out["signal"] = out["signal"].astype(float)
    out["fuzz"] = out["fuzz"].astype(float)

    if chrom_allow_regex.strip() != "":
        pat = re.compile(chrom_allow_regex)
        out = out[out["chr"].astype(str).map(lambda x: pat.match(x) is not None)].copy()

    if not keep_all_centers:
        out = (out.sort_values(["chr", "center", "signal"], ascending=[True, True, False])
                 .drop_duplicates(subset=["chr", "center"], keep="first")
                 .reset_index(drop=True))

    out = out.sort_values(["chr", "center"]).reset_index(drop=True)
    return out


def merge_intervals(df_regions: pd.DataFrame, allow_gap_bp: int = 0) -> pd.DataFrame:
    if df_regions.empty:
        return df_regions

    df = df_regions.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    merged = []
    cur_chr = None
    cur_s = None
    cur_e = None

    for row in df.itertuples(index=False):
        ch = row.chr
        s = int(row.start)
        e = int(row.end)

        if cur_chr is None:
            cur_chr, cur_s, cur_e = ch, s, e
            continue

        if ch != cur_chr or s > (cur_e + allow_gap_bp):
            merged.append((cur_chr, cur_s, cur_e))
            cur_chr, cur_s, cur_e = ch, s, e
        else:
            if e > cur_e:
                cur_e = e

    merged.append((cur_chr, cur_s, cur_e))
    return pd.DataFrame(merged, columns=["chr", "start", "end"])


def region_stats_for_center_span(
    chrom: str,
    centers: np.ndarray,
    signals: np.ndarray,
    fuzzes: np.ndarray,
    regions_chr: pd.DataFrame,
    spacing_min: float,
    spacing_max: float,
) -> pd.DataFrame:
    rows = []
    for r in regions_chr.itertuples(index=False):
        rs = int(r.start)
        re_ = int(r.end)
        L = int(np.searchsorted(centers, rs, side="left"))
        R = int(np.searchsorted(centers, re_, side="right")) - 1

        if L < 0:
            L = 0
        if R >= centers.size:
            R = centers.size - 1

        if R < L:
            rows.append({
                "chr": chrom,
                "region_start": rs,
                "region_end": re_,
                "length_bp": int(re_ - rs),
                "n_nuc": 0,
                "mean_spacing": np.nan,
                "cv": np.nan,
                "good_ratio": np.nan,
                "mean_signal": np.nan,
                "p20_signal": np.nan,
                "min_signal": np.nan,
                "median_fuzz": np.nan,
                "p80_fuzz": np.nan,
                "score": np.nan,
                "L": L,
                "R": R,
            })
            continue

        seg_centers = centers[L:R + 1]
        seg_signal = signals[L:R + 1]
        seg_fuzz = fuzzes[L:R + 1]
        n_nuc = int(seg_centers.size)

        if n_nuc >= 2:
            gaps = np.diff(seg_centers).astype(float)
            mean_spacing = float(np.mean(gaps))
            sd = float(np.std(gaps, ddof=1)) if gaps.size >= 2 else 0.0
            cv = float(sd / mean_spacing) if mean_spacing > 0 else np.inf
            good_ratio = float(np.mean((gaps >= spacing_min) & (gaps <= spacing_max)))
        else:
            mean_spacing = np.nan
            cv = np.nan
            good_ratio = np.nan

        mean_signal = float(np.mean(seg_signal)) if n_nuc > 0 else np.nan
        p20_signal = float(np.quantile(seg_signal, 0.20)) if n_nuc > 0 else np.nan
        min_signal = float(np.min(seg_signal)) if n_nuc > 0 else np.nan
        median_fuzz = float(np.median(seg_fuzz)) if n_nuc > 0 else np.nan
        p80_fuzz = float(np.quantile(seg_fuzz, 0.80)) if n_nuc > 0 else np.nan
        length_bp = int(re_ - rs)

        if np.isfinite(cv) and np.isfinite(good_ratio) and np.isfinite(mean_signal) and np.isfinite(median_fuzz):
            score = float((good_ratio * 100.0) - (cv * 50.0) + np.log1p(mean_signal) - 0.2 * np.log1p(median_fuzz))
        else:
            score = np.nan

        rows.append({
            "chr": chrom,
            "region_start": rs,
            "region_end": re_,
            "length_bp": length_bp,
            "n_nuc": n_nuc,
            "mean_spacing": mean_spacing,
            "cv": cv,
            "good_ratio": good_ratio,
            "mean_signal": mean_signal,
            "p20_signal": p20_signal,
            "min_signal": min_signal,
            "median_fuzz": median_fuzz,
            "p80_fuzz": p80_fuzz,
            "score": score,
            "L": int(L),
            "R": int(R),
        })

    return pd.DataFrame(rows)


def calc_windows_for_chr(df_chr: pd.DataFrame, k: int) -> Optional[pd.DataFrame]:
    df_chr = df_chr.sort_values("center").reset_index(drop=True)
    centers = df_chr["center"].to_numpy(dtype=np.int64)
    fuzz = df_chr["fuzz"].to_numpy(dtype=float)
    signal = df_chr["signal"].to_numpy(dtype=float)

    n = centers.shape[0]
    if n < k:
        return None

    d = np.diff(centers).astype(float)
    S = np.concatenate(([0.0], np.cumsum(d)))
    S2 = np.concatenate(([0.0], np.cumsum(d * d)))

    m = n - k + 1
    idx = np.arange(m, dtype=np.int64)
    a = idx
    b = idx + k - 2

    denom = float(k - 1)
    sum_d = S[b + 1] - S[a]
    sum2_d = S2[b + 1] - S2[a]

    mean_spacing = sum_d / denom
    var_spacing = (sum2_d / denom) - mean_spacing ** 2
    var_spacing[var_spacing < 0] = 0.0
    sd_spacing = np.sqrt(var_spacing)
    cv_spacing = np.where(mean_spacing > 0, sd_spacing / mean_spacing, np.nan)

    start = centers[0:m]
    end = centers[idx + k - 1]

    sig_s = pd.Series(signal)
    fuzz_s = pd.Series(fuzz)
    mean_signal = sig_s.rolling(window=k).mean().to_numpy()[k - 1:]
    min_signal = sig_s.rolling(window=k).min().to_numpy()[k - 1:]
    median_fuzz = fuzz_s.rolling(window=k).median().to_numpy()[k - 1:]

    return pd.DataFrame({
        "chr": df_chr["chr"].iloc[0],
        "start": start,
        "end": end,
        "n_nuc": k,
        "mean_spacing": mean_spacing,
        "sd_spacing": sd_spacing,
        "cv_spacing": cv_spacing,
        "mean_signal": mean_signal,
        "min_signal": min_signal,
        "median_fuzz": median_fuzz,
    })


def run_mergewindow(df: pd.DataFrame, args, sample_prefix: str):
    all_w = []
    chrom_groups = list(df.groupby("chr", sort=False))
    for chrom, df_chr in tqdm(chrom_groups, total=len(chrom_groups), desc=f"Sliding windows (k={args.k})"):
        w = calc_windows_for_chr(df_chr, k=args.k)
        if w is not None and not w.empty:
            all_w.append(w)

    if not all_w:
        raise RuntimeError("No windows generated. Try lowering --k.")

    windows = pd.concat(all_w, ignore_index=True)
    print(f"[INFO] Total windows: {windows.shape[0]}", file=sys.stderr)

    good = windows[
        (windows["mean_spacing"] >= args.spacing_min)
        & (windows["mean_spacing"] <= args.spacing_max)
        & (windows["cv_spacing"] <= args.cv_cutoff)
    ].copy()
    print(f"[INFO] After spacing+CV: {good.shape[0]}", file=sys.stderr)
    if good.empty:
        raise RuntimeError("No windows pass spacing+CV. Relax thresholds.")

    mean_thr = good["mean_signal"].quantile(args.signal_mean_quantile)
    min_thr = good["min_signal"].quantile(args.signal_min_quantile)
    good = good[(good["mean_signal"] >= mean_thr) & (good["min_signal"] >= min_thr)].copy()
    print(f"[INFO] After signal: {good.shape[0]} (mean_thr={mean_thr:.4g}, min_thr={min_thr:.4g})", file=sys.stderr)
    if good.empty:
        raise RuntimeError("No windows pass signal filters. Lower quantiles.")

    fuzz_thr = good["median_fuzz"].quantile(args.fuzz_quantile)
    good = good[good["median_fuzz"] <= fuzz_thr].copy()
    print(f"[INFO] After fuzz: {good.shape[0]} (fuzz_thr={fuzz_thr:.4g})", file=sys.stderr)
    if good.empty:
        raise RuntimeError("No windows pass fuzz filter. Increase --fuzz_quantile.")

    merged = merge_intervals(good[["chr", "start", "end"]].copy(), allow_gap_bp=args.merge_allow_gap_bp)
    print(f"[INFO] Merged regions: {merged.shape[0]}", file=sys.stderr)
    if args.min_region_bp > 0:
        merged = merged[(merged["end"] - merged["start"] + 1) >= args.min_region_bp].copy()
        print(f"[INFO] After min_region_bp={args.min_region_bp}: {merged.shape[0]}", file=sys.stderr)
    if merged.empty:
        raise RuntimeError("No merged regions left.")

    nuc_by_chr = {}
    for chrom, df_chr in df.groupby("chr", sort=False):
        df_chr = df_chr.sort_values("center").reset_index(drop=True)
        nuc_by_chr[chrom] = (
            df_chr["center"].to_numpy(np.int64),
            df_chr["signal"].to_numpy(float),
            df_chr["fuzz"].to_numpy(float),
        )

    stats_all = []
    merged_groups = list(merged.groupby("chr", sort=False))
    for chrom, reg_chr in tqdm(merged_groups, total=len(merged_groups), desc="Region stats"):
        centers, signals, fuzzes = nuc_by_chr[chrom]
        st = region_stats_for_center_span(
            chrom=chrom,
            centers=centers,
            signals=signals,
            fuzzes=fuzzes,
            regions_chr=reg_chr,
            spacing_min=args.spacing_min,
            spacing_max=args.spacing_max,
        )
        stats_all.append(st)

    stats_df = pd.concat(stats_all, ignore_index=True)
    stats_df = stats_df.sort_values(["chr", "region_start", "region_end"]).reset_index(drop=True)
    stats_df.insert(0, "region_id", [f"{sample_prefix}_region_{i + 1}" for i in range(stats_df.shape[0])])

    tag = (
        f"{sample_prefix}.mergewindow"
        f".k{args.k}"
        f".sp{fmt_num(args.spacing_min)}-{fmt_num(args.spacing_max)}"
        f".cv{fmt_num(args.cv_cutoff)}"
        f".meanQ{fmt_num(args.signal_mean_quantile)}"
        f".minQ{fmt_num(args.signal_min_quantile)}"
        f".fuzzQ{fmt_num(args.fuzz_quantile)}"
        f".gap{fmt_num(args.merge_allow_gap_bp)}"
        f".minLen{fmt_num(args.min_region_bp)}"
    )
    return stats_df, tag, True


def segment_metrics(pos, starts, ends, signals, fuzzes, L, R, args) -> Optional[dict]:
    if L < 0 or R >= len(pos) or L >= R:
        return None

    seg_pos = pos[L:R + 1]
    seg_sig = signals[L:R + 1]
    seg_fuz = fuzzes[L:R + 1]
    seg_st = starts[L:R + 1]
    seg_en = ends[L:R + 1]
    gaps = np.diff(seg_pos)
    n_nuc = len(seg_pos)
    if len(gaps) == 0:
        return None

    mean_spacing = float(np.mean(gaps))
    sd = float(np.std(gaps, ddof=1)) if len(gaps) >= 2 else 0.0
    cv = float(sd / mean_spacing) if mean_spacing > 0 else np.inf
    good_ratio = float(np.mean((gaps >= args.spacing_min) & (gaps <= args.spacing_max)))

    mean_signal = float(np.mean(seg_sig))
    p20_signal = float(np.quantile(seg_sig, 0.20))
    min_signal = float(np.min(seg_sig))
    median_fuzz = float(np.median(seg_fuz))
    p80_fuzz = float(np.quantile(seg_fuz, 0.80))
    region_start = int(np.min(seg_st))
    region_end = int(np.max(seg_en))
    length_bp = int(region_end - region_start)
    score = float((good_ratio * 100.0) - (cv * 50.0) + math.log1p(mean_signal) - 0.2 * math.log1p(median_fuzz))

    return {
        "L": int(L), "R": int(R),
        "n_nuc": int(n_nuc),
        "region_start": region_start,
        "region_end": region_end,
        "length_bp": length_bp,
        "mean_spacing": mean_spacing,
        "cv": cv,
        "good_ratio": good_ratio,
        "mean_signal": mean_signal,
        "p20_signal": p20_signal,
        "min_signal": min_signal,
        "median_fuzz": median_fuzz,
        "p80_fuzz": p80_fuzz,
        "score": score,
    }


def passes_basic_extend(metrics: Optional[dict], args) -> bool:
    if metrics is None:
        return False
    if metrics["good_ratio"] < args.min_good_ratio_extend:
        return False
    if metrics["cv"] > (args.cv_cutoff * 1.5):
        return False
    return True


def passes_basic_final(metrics: Optional[dict], args) -> bool:
    if metrics is None:
        return False
    if metrics["n_nuc"] < args.min_nuc_final:
        return False
    if metrics["good_ratio"] < args.min_good_ratio_final:
        return False
    if metrics["cv"] > args.cv_cutoff:
        return False
    if not (args.spacing_min <= metrics["mean_spacing"] <= args.spacing_max):
        return False
    if metrics["length_bp"] < args.min_region_bp:
        return False
    return True


def find_seed_runs(pos: np.ndarray, args):
    n = len(pos)
    if n < args.k:
        return []
    gaps = np.diff(pos)
    good_gap = (gaps >= args.spacing_min) & (gaps <= args.spacing_max)

    if args.k < 3:
        raise ValueError("seedextend method requires --k >= 3")

    seed_mask = np.zeros(n - args.k + 1, dtype=bool)
    for i in range(n - args.k + 1):
        seg = gaps[i:i + args.k - 1]
        if seg.size != args.k - 1:
            continue
        if not np.all(good_gap[i:i + args.k - 1]):
            continue
        if float(np.max(seg) - np.min(seg)) <= args.seed_delta_bp:
            seed_mask[i] = True

    runs = []
    i = 0
    while i < len(seed_mask):
        if not seed_mask[i]:
            i += 1
            continue
        s = i
        while i + 1 < len(seed_mask) and seed_mask[i + 1]:
            i += 1
        e = i
        L = s
        R = e + args.k - 1
        runs.append((L, R))
        i += 1
    return runs


def extend_segment(pos, starts, ends, signals, fuzzes, L0, R0, args) -> Optional[dict]:
    n = len(pos)
    L, R = max(0, int(L0)), min(n - 1, int(R0))
    if R - L + 1 < args.k:
        return None

    fail_left = 0
    fail_right = 0
    cur = segment_metrics(pos, starts, ends, signals, fuzzes, L, R, args)
    if cur is None:
        return None

    while True:
        can_left = (L > 0) and (fail_left < args.extend_max_fail_per_side)
        can_right = (R < n - 1) and (fail_right < args.extend_max_fail_per_side)
        if not can_left and not can_right:
            break

        best_side = None
        best_metrics = None

        if can_left:
            m_left = segment_metrics(pos, starts, ends, signals, fuzzes, L - 1, R, args)
            if passes_basic_extend(m_left, args):
                best_side = "L"
                best_metrics = m_left
            else:
                fail_left += 1

        if can_right:
            m_right = segment_metrics(pos, starts, ends, signals, fuzzes, L, R + 1, args)
            if passes_basic_extend(m_right, args):
                if best_metrics is None:
                    best_side = "R"
                    best_metrics = m_right
                else:
                    key_best = (best_metrics["good_ratio"], -best_metrics["cv"], best_metrics["mean_signal"], -best_metrics["median_fuzz"])
                    key_new = (m_right["good_ratio"], -m_right["cv"], m_right["mean_signal"], -m_right["median_fuzz"])
                    if key_new > key_best:
                        best_side = "R"
                        best_metrics = m_right
            else:
                fail_right += 1

        if best_side is None:
            if fail_left >= args.extend_max_fail_per_side and fail_right >= args.extend_max_fail_per_side:
                break
            continue

        if best_side == "L":
            L -= 1
            fail_left = 0
        else:
            R += 1
            fail_right = 0
        cur = best_metrics

    return segment_metrics(pos, starts, ends, signals, fuzzes, L, R, args)


def dedup_segments(seg_df: pd.DataFrame) -> pd.DataFrame:
    seg_df = seg_df.sort_values(["chr", "region_start", "region_end", "score"], ascending=[True, True, True, False]).copy()
    seg_df = seg_df.drop_duplicates(subset=["chr", "region_start", "region_end"], keep="first").copy()
    return seg_df.reset_index(drop=True)


def merge_segments_by_coord(seg_df: pd.DataFrame, args) -> pd.DataFrame:
    if seg_df.empty:
        return seg_df

    merged_rows = []
    for chrom, sub in seg_df.groupby("chr", sort=True):
        sub = sub.sort_values("region_start").copy()
        cur = sub.iloc[0].to_dict()

        for _, row in sub.iloc[1:].iterrows():
            gap = int(row["region_start"] - cur["region_end"])
            if gap <= args.merge_allow_gap_bp:
                new_start = int(min(cur["region_start"], row["region_start"]))
                new_end = int(max(cur["region_end"], row["region_end"]))
                n1 = int(cur["n_nuc"])
                n2 = int(row["n_nuc"])
                w_gap1 = max(n1 - 1, 1)
                w_gap2 = max(n2 - 1, 1)

                merged = dict(cur)
                merged["region_start"] = new_start
                merged["region_end"] = new_end
                merged["length_bp"] = int(new_end - new_start)
                merged["n_nuc"] = n1 + n2
                merged["mean_spacing"] = float((cur["mean_spacing"] * w_gap1 + row["mean_spacing"] * w_gap2) / (w_gap1 + w_gap2))
                merged["cv"] = float(max(cur["cv"], row["cv"]))
                merged["good_ratio"] = float(min(cur["good_ratio"], row["good_ratio"]))
                merged["mean_signal"] = float((cur["mean_signal"] * n1 + row["mean_signal"] * n2) / (n1 + n2))
                merged["p20_signal"] = float(min(cur["p20_signal"], row["p20_signal"]))
                merged["min_signal"] = float(min(cur["min_signal"], row["min_signal"]))
                merged["median_fuzz"] = float(max(cur["median_fuzz"], row["median_fuzz"]))
                merged["p80_fuzz"] = float(max(cur["p80_fuzz"], row["p80_fuzz"]))
                merged["score"] = float(min(cur["score"], row["score"]))

                if passes_basic_final(merged, args):
                    cur = merged
                else:
                    merged_rows.append(cur)
                    cur = row.to_dict()
            else:
                merged_rows.append(cur)
                cur = row.to_dict()

        merged_rows.append(cur)

    return pd.DataFrame(merged_rows).reset_index(drop=True)


def run_seedextend(df: pd.DataFrame, args, sample_prefix: str):
    all_segments = []
    chrom_groups = list(df.groupby("chr", sort=True))
    for chrom, sub in tqdm(chrom_groups, total=len(chrom_groups), desc="Seed+extend per chromosome"):
        sub = sub.sort_values("center").copy()
        pos = sub["center"].to_numpy(dtype=float)
        starts = sub["start"].to_numpy(dtype=float)
        ends = sub["end"].to_numpy(dtype=float)
        signals = sub["signal"].to_numpy(dtype=float)
        fuzzes = sub["fuzz"].to_numpy(dtype=float)

        if len(pos) < args.k:
            continue

        seed_runs = find_seed_runs(pos, args)
        if not seed_runs:
            continue

        for (L0, R0) in tqdm(seed_runs, desc=f"{chrom}: extend seeds", leave=False):
            met = extend_segment(pos, starts, ends, signals, fuzzes, L0, R0, args)
            if met is not None:
                all_segments.append({"chr": chrom, **met})

    seg_df = pd.DataFrame(all_segments)
    print(f"[INFO] Segments from seed+extend (before dedup): {seg_df.shape}", file=sys.stderr)
    if seg_df.empty:
        raise RuntimeError("No segments found. Consider loosening seed/extend thresholds.")

    seg_df = dedup_segments(seg_df)
    print(f"[INFO] After dedup: {seg_df.shape}", file=sys.stderr)

    seg_basic = seg_df[seg_df.apply(lambda x: passes_basic_final(x, args), axis=1)].copy()
    print(f"[INFO] After basic final filters: {seg_basic.shape}", file=sys.stderr)
    if seg_basic.empty:
        raise RuntimeError("No segments remain after basic final filters.")

    seg_merged = merge_segments_by_coord(seg_basic, args)
    print(f"[INFO] After merge: {seg_merged.shape}", file=sys.stderr)

    seg_merged = seg_merged[seg_merged.apply(lambda x: passes_basic_final(x, args), axis=1)].copy().reset_index(drop=True)
    print(f"[INFO] After re-check final filters: {seg_merged.shape}", file=sys.stderr)
    if seg_merged.empty:
        raise RuntimeError("No segments remain after merge and final re-check.")

    if args.min_mean_signal_hard is not None:
        seg_merged = seg_merged[seg_merged["mean_signal"] >= float(args.min_mean_signal_hard)].copy()
    if args.min_p20_signal_hard is not None:
        seg_merged = seg_merged[seg_merged["p20_signal"] >= float(args.min_p20_signal_hard)].copy()
    if args.max_median_fuzz_hard is not None:
        seg_merged = seg_merged[seg_merged["median_fuzz"] <= float(args.max_median_fuzz_hard)].copy()
    seg_merged = seg_merged.reset_index(drop=True)
    print(f"[INFO] After hard cutoffs: {seg_merged.shape}", file=sys.stderr)
    if seg_merged.empty:
        raise RuntimeError("All segments filtered out by hard cutoffs.")

    thr_mean_sig = float(seg_merged["mean_signal"].quantile(args.signal_mean_quantile))
    thr_p20_sig = float(seg_merged["p20_signal"].quantile(args.signal_p20_quantile))
    thr_med_fuzz = float(seg_merged["median_fuzz"].quantile(args.fuzz_quantile))

    seg_final = seg_merged[
        (seg_merged["mean_signal"] >= thr_mean_sig)
        & (seg_merged["p20_signal"] >= thr_p20_sig)
        & (seg_merged["median_fuzz"] <= thr_med_fuzz)
    ].copy().reset_index(drop=True)
    print(
        f"[INFO] Quantile thresholds: mean_signal>={thr_mean_sig:.4g}, p20_signal>={thr_p20_sig:.4g}, median_fuzz<={thr_med_fuzz:.4g}",
        file=sys.stderr,
    )
    print(f"[INFO] After quantile filters: {seg_final.shape}", file=sys.stderr)
    if seg_final.empty:
        raise RuntimeError("All segments filtered out at the final stage. Try relaxing quantile thresholds.")

    seg_final = seg_final.sort_values(["chr", "region_start", "region_end"]).reset_index(drop=True)
    seg_final.insert(0, "region_id", [f"{sample_prefix}_region_{i + 1}" for i in range(seg_final.shape[0])])

    tag = (
        f"{sample_prefix}.seedextend"
        f".k{args.k}"
        f".sp{fmt_num(args.spacing_min)}-{fmt_num(args.spacing_max)}"
        f".cv{fmt_num(args.cv_cutoff)}"
        f".sd{fmt_num(args.seed_delta_bp)}"
        f".extFail{fmt_num(args.extend_max_fail_per_side)}"
        f".grE{fmt_num(args.min_good_ratio_extend)}"
        f".grF{fmt_num(args.min_good_ratio_final)}"
        f".minN{fmt_num(args.min_nuc_final)}"
        f".gap{fmt_num(args.merge_allow_gap_bp)}"
        f".minBP{fmt_num(args.min_region_bp)}"
        f".meanQ{fmt_num(args.signal_mean_quantile)}"
        f".p20Q{fmt_num(args.signal_p20_quantile)}"
        f".fuzzQ{fmt_num(args.fuzz_quantile)}"
    )
    return seg_final, tag, False


def write_outputs(stats_df: pd.DataFrame, out_dir: str, tag: str, center_span_bed: bool) -> Tuple[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    bed_path = os.path.join(out_dir, f"{tag}.bed")
    stats_path = os.path.join(out_dir, f"{tag}.stats.tsv")

    bed_df = stats_df[["chr", "region_start", "region_end", "region_id"]].copy()
    if center_span_bed:
        bed_df["region_start"] = bed_df["region_start"].astype(int) - 1
        bed_df.loc[bed_df["region_start"] < 0, "region_start"] = 0

    bed_df["region_start"] = bed_df["region_start"].astype(int)
    bed_df["region_end"] = bed_df["region_end"].astype(int)
    bed_df.to_csv(bed_path, sep="\t", header=False, index=False)

    cols = [
        "region_id", "chr", "region_start", "region_end", "length_bp",
        "n_nuc", "mean_spacing", "cv", "good_ratio",
        "mean_signal", "p20_signal", "min_signal",
        "median_fuzz", "p80_fuzz", "score", "L", "R",
    ]
    stats_df[cols].to_csv(stats_path, sep="\t", index=False)
    print(f"[INFO] Saved BED: {bed_path}", file=sys.stderr)
    print(f"[INFO] Saved stats: {stats_path}", file=sys.stderr)
    return bed_path, stats_path



def read_bed3_ignore_extra(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    if df.shape[1] < 3:
        raise ValueError(f"BED must have >=3 columns: {path}")
    df = df.iloc[:, :3].copy()
    df.columns = ["chr", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"]).copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df = df[df["end"] > df["start"]].copy().reset_index(drop=True)
    df["row_idx"] = np.arange(df.shape[0], dtype=np.int64)
    return df


def read_peak_bed3(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, comment="#", engine="python", dtype={0: str})
    if df.shape[1] < 3:
        raise ValueError(f"Peak file must have >=3 columns: {path}")
    df = df.iloc[:, :3].copy()
    df.columns = ["chr", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"]).copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df = df[df["end"] > df["start"]].copy().reset_index(drop=True)
    return df


def union_overlap_len_with_intervals(a_start: int, a_end: int, b_starts: np.ndarray, b_ends: np.ndarray) -> int:
    s = np.maximum(b_starts, a_start)
    e = np.minimum(b_ends, a_end)
    m = e > s
    if not np.any(m):
        return 0
    s = s[m]
    e = e[m]
    order = np.argsort(s, kind="mergesort")
    s = s[order]
    e = e[order]
    total = 0
    cur_s = int(s[0])
    cur_e = int(e[0])
    for i in range(1, len(s)):
        si = int(s[i])
        ei = int(e[i])
        if si > cur_e:
            total += (cur_e - cur_s)
            cur_s, cur_e = si, ei
        elif ei > cur_e:
            cur_e = ei
    total += (cur_e - cur_s)
    return int(total)


def filter_phased_by_peak_overlap(phased: pd.DataFrame, peaks: pd.DataFrame, min_frac: float = 0.30) -> pd.DataFrame:
    out_rows = []
    phased_groups = list(phased.groupby("chr", sort=False))
    peak_by_chr = {c: sub.sort_values("start").reset_index(drop=True) for c, sub in peaks.groupby("chr", sort=False)}

    for chrom, ph_chr in tqdm(phased_groups, desc="Overlap per chromosome", total=len(phased_groups)):
        ph_chr = ph_chr.sort_values("start").reset_index(drop=True)
        pk_chr = peak_by_chr.get(chrom, None)

        if pk_chr is None or pk_chr.empty:
            for row in ph_chr.itertuples(index=False):
                a_len = int(row.end - row.start)
                out_rows.append((chrom, int(row.start), int(row.end), int(row.row_idx), 0, a_len, 0.0, False))
            continue

        pk_st = pk_chr["start"].to_numpy(np.int64)
        pk_en = pk_chr["end"].to_numpy(np.int64)

        for row in ph_chr.itertuples(index=False):
            a_s = int(row.start)
            a_e = int(row.end)
            a_len = a_e - a_s
            ridx = int(row.row_idx)

            if a_len <= 0:
                out_rows.append((chrom, a_s, a_e, ridx, 0, max(a_len, 1), 0.0, False))
                continue

            j = np.searchsorted(pk_st, a_e, side="left")
            if j <= 0:
                out_rows.append((chrom, a_s, a_e, ridx, 0, a_len, 0.0, False))
                continue

            cand_st = pk_st[:j]
            cand_en = pk_en[:j]
            m = cand_en > a_s
            if not np.any(m):
                out_rows.append((chrom, a_s, a_e, ridx, 0, a_len, 0.0, False))
                continue

            ov = union_overlap_len_with_intervals(a_s, a_e, cand_st[m], cand_en[m])
            frac = float(ov / a_len)
            keep = frac >= float(min_frac)
            out_rows.append((chrom, a_s, a_e, ridx, int(ov), int(a_len), frac, keep))

    return pd.DataFrame(
        out_rows,
        columns=["chr", "start", "end", "row_idx", "overlap_bp", "phased_len", "overlap_frac", "keep"],
    )


def locate_peak_bed(peak_dir: str) -> str:
    all_files = sorted(glob.glob(os.path.join(peak_dir, "**", "*"), recursive=True))
    all_files = [x for x in all_files if os.path.isfile(x)]

    if not all_files:
        raise FileNotFoundError(f"No files found under SICER output directory: {peak_dir}")

    bed_files = [x for x in all_files if x.lower().endswith(".bed")]
    scoreisland_files = [x for x in all_files if x.lower().endswith(".scoreisland")]

    def bed_score(path: str) -> Tuple[int, int, str]:
        base = os.path.basename(path).lower()
        if base.endswith("island.bed"):
            return (0, len(base), path)
        if "island" in base and "filtered" not in base and "significant" not in base:
            return (1, len(base), path)
        if base.endswith("islandfiltered.bed"):
            return (2, len(base), path)
        return (3, len(base), path)

    if bed_files:
        ranked = sorted(bed_files, key=bed_score)
        return ranked[0]

    if scoreisland_files:
        ranked = sorted(scoreisland_files, key=lambda x: (len(os.path.basename(x)), x))
        return ranked[0]

    raise FileNotFoundError(
        f"No SICER island-like output found under: {peak_dir}\n"
        f"Expected either *.bed (with control) or *.scoreisland (without control)."
    )


def run_sicer2(args, peak_dir: str) -> Tuple[str, str]:
    os.makedirs(peak_dir, exist_ok=True)
    cmd = [
        args.sicer_exe,
        "-t", args.treatment,
        "-s", args.species,
        "-o", peak_dir,
        "-w", str(args.sicer_window_size),
        "-rt", str(args.sicer_redundancy_threshold),
        "-f", str(args.sicer_fragment_size),
        "-egf", str(args.sicer_effective_genome_fraction),
        "-fdr", str(args.sicer_fdr),
        "-g", str(args.sicer_gap_size),
    ]
    if args.control:
        cmd.extend(["-c", args.control])
    else:
        cmd.extend(["-e", str(int(args.sicer_e_value))])
    if args.sicer_cpu is not None:
        cmd.extend(["-cpu", str(args.sicer_cpu)])
    if args.sicer_significant_reads:
        cmd.append("--significant_reads")
    if args.sicer_extra_args:
        cmd.extend(shlex.split(args.sicer_extra_args))

    print("[INFO] Running SICER2:", " ".join(shlex.quote(x) for x in cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    peak_bed = locate_peak_bed(peak_dir)
    print(f"[INFO] Using peak BED: {peak_bed}", file=sys.stderr)
    return peak_bed, " ".join(shlex.quote(x) for x in cmd)


def filter_outputs_by_peaks(
    phased_bed: str,
    phased_stats: str,
    peak_bed: str,
    out_dir: str,
    min_overlap_frac: float,
) -> Tuple[str, str, str]:
    os.makedirs(out_dir, exist_ok=True)
    ph = read_bed3_ignore_extra(phased_bed)
    pk = read_peak_bed3(peak_bed)
    ov_df = filter_phased_by_peak_overlap(ph, pk, min_frac=min_overlap_frac)

    keep_row_idx = set(ov_df.loc[ov_df["keep"], "row_idx"].astype(int).tolist())
    ph_keep = ph[ph["row_idx"].isin(keep_row_idx)].copy()

    st = pd.read_csv(phased_stats, sep="\t")
    st["_row_idx"] = np.arange(st.shape[0], dtype=np.int64)
    st_keep = st[st["_row_idx"].isin(keep_row_idx)].drop(columns=["_row_idx"]).copy()

    base = os.path.basename(phased_bed).replace(".bed", "")
    tag = f"peakOvFrac{str(min_overlap_frac).replace('.', 'p')}"
    out_bed = os.path.join(out_dir, f"{base}.{tag}.bed")
    out_stats = os.path.join(out_dir, f"{base}.{tag}.stats.tsv")
    out_overlap = os.path.join(out_dir, f"{base}.{tag}.overlap.tsv")

    ph_keep[["chr", "start", "end"]].to_csv(out_bed, sep="\t", header=False, index=False)
    st_keep.to_csv(out_stats, sep="\t", index=False)
    ov_df.to_csv(out_overlap, sep="\t", index=False)

    print(f"[INFO] Saved peak-filtered BED: {out_bed}", file=sys.stderr)
    print(f"[INFO] Saved peak-filtered stats: {out_stats}", file=sys.stderr)
    print(f"[INFO] Saved overlap report: {out_overlap}", file=sys.stderr)
    return out_bed, out_stats, out_overlap


def write_run_metadata(
    out_dir: str,
    args,
    raw_bed: str,
    raw_stats: str,
    peak_bed: str,
    filtered_bed: str,
    filtered_stats: str,
    overlap_tsv: str,
    sicer_cmd: str,
) -> str:
    meta_path = os.path.join(out_dir, "run_metadata.txt")
    lines = [
        f"input={args.input}",
        f"method={args.method}",
        f"output_dir={args.output_dir}",
        f"raw_bed={raw_bed}",
        f"raw_stats={raw_stats}",
        f"peak_mode={args.peak_mode}",
        f"peak_bed={peak_bed}",
        f"filtered_bed={filtered_bed}",
        f"filtered_stats={filtered_stats}",
        f"overlap_tsv={overlap_tsv}",
        f"min_overlap_frac={args.min_overlap_frac}",
        f"sicer_command={sicer_cmd if sicer_cmd else 'NA'}",
    ]
    with open(meta_path, "w") as f:
        for line in lines:
            f.write(line + "\n")
    return meta_path

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Run SICER2 broad peak calling, identify well-phased nucleosome arrays, and optionally keep only arrays overlapping SICER peaks."
    )
    ap.add_argument("-i", "--input", required=True, help="DANPOS output: *positions*.xls")
    ap.add_argument("-O", "--output_dir", required=True, help="Output directory")

    ap.add_argument("--method", choices=["seedextend", "mergewindow"], default="seedextend", help="Array-identification method (default: seedextend)")
    ap.add_argument("--chrom_allow_regex", type=str, default=r"^chr([0-9]+|X|Y)$", help="Keep chromosomes matching this regex. Use empty string to disable.")
    ap.add_argument("--keep_all_centers", action="store_true", help="Do not deduplicate duplicate centers within chromosome.")
    ap.add_argument("--k", type=int, default=3, help="Window/seed size in nucleosomes for both methods (default: 3)")
    ap.add_argument("--spacing_min", type=float, default=160.0, help="Minimum spacing bp")
    ap.add_argument("--spacing_max", type=float, default=220.0, help="Maximum spacing bp")
    ap.add_argument("--cv_cutoff", type=float, default=0.15, help="Maximum spacing CV")
    ap.add_argument("--merge_allow_gap_bp", type=int, default=100, help="Merge neighboring intervals if gap <= this value")
    ap.add_argument("--min_region_bp", type=int, default=0, help="Minimum retained region length in bp")
    ap.add_argument("--signal_mean_quantile", type=float, default=0.5, help="Quantile cutoff for mean signal")
    ap.add_argument("--signal_min_quantile", type=float, default=0.5, help="Quantile cutoff for minimum signal (mergewindow only)")
    ap.add_argument("--signal_p20_quantile", type=float, default=0.5, help="Quantile cutoff for p20 signal (seedextend only)")
    ap.add_argument("--fuzz_quantile", type=float, default=None, help="Quantile cutoff for fuzz; default is 0.3 for seedextend and 0.1 for mergewindow")
    ap.add_argument("--seed_delta_bp", type=float, default=60.0, help="Maximum allowed max-min gap difference within a seed (seedextend only)")
    ap.add_argument("--extend_max_fail_per_side", type=int, default=3, help="Stop extension after this many consecutive failures on a side (seedextend only)")
    ap.add_argument("--min_good_ratio_extend", type=float, default=0.60, help="Minimum good-ratio during extension (seedextend only)")
    ap.add_argument("--min_good_ratio_final", type=float, default=0.60, help="Minimum good-ratio for final regions (seedextend only)")
    ap.add_argument("--min_nuc_final", type=int, default=3, help="Minimum nucleosome count for final regions (seedextend only)")
    ap.add_argument("--min_mean_signal_hard", type=float, default=None, help="Optional hard cutoff on mean signal (seedextend only)")
    ap.add_argument("--min_p20_signal_hard", type=float, default=None, help="Optional hard cutoff on p20 signal (seedextend only)")
    ap.add_argument("--max_median_fuzz_hard", type=float, default=None, help="Optional hard cutoff on median fuzz (seedextend only)")

    ap.add_argument("--peak_mode", choices=["sicer2", "existing"], default="sicer2", help="How to obtain peak BED (default: sicer2)")
    ap.add_argument("--peak_bed", type=str, default=None, help="Existing peak BED to use when --peak_mode existing")
    ap.add_argument("--min_overlap_frac", type=float, default=0.30, help="Minimum fraction of an array region that must be covered by peak BED (default: 0.30)")

    ap.add_argument("--sicer_exe", type=str, default="sicer", help="SICER2 executable (default: sicer)")
    ap.add_argument("--treatment", type=str, default=None, help="Treatment BED/BAM for SICER2")
    ap.add_argument("--control", type=str, default=None, help="Control BED/BAM for SICER2")
    ap.add_argument("--species", type=str, default=None, help="SICER2 species key, e.g. hg38 or mm10")
    ap.add_argument("--sicer_window_size", type=int, default=200, help="SICER2 -w / --window_size (default: 200)")
    ap.add_argument("--sicer_gap_size", type=int, default=600, help="SICER2 -g / --gap_size (default: 600)")
    ap.add_argument("--sicer_fdr", type=float, default=0.01, help="SICER2 -fdr / --false_discovery_rate (default: 0.01)")
    ap.add_argument("--sicer_redundancy_threshold", type=int, default=1, help="SICER2 -rt / --redundancy_threshold (default: 1)")
    ap.add_argument("--sicer_fragment_size", type=int, default=150, help="SICER2 -f / --fragment_size (default: 150)")
    ap.add_argument("--sicer_effective_genome_fraction", type=float, default=0.74, help="SICER2 -egf / --effective_genome_faction (default: 0.74)")
    ap.add_argument("--sicer_e_value", type=int, default=1000, help="SICER2 -e / --e_value, used when no control is provided (default: 1000)")
    ap.add_argument("--sicer_cpu", type=int, default=None, help="SICER2 -cpu / --cpu")
    ap.add_argument("--sicer_significant_reads", action="store_true", help="Pass --significant_reads to SICER2")
    ap.add_argument("--sicer_extra_args", type=str, default="", help='Extra SICER2 arguments as a single quoted string, e.g. "--foo 1 --bar 2"')

    return ap


def main():
    args = build_parser().parse_args()
    if args.fuzz_quantile is None:
        args.fuzz_quantile = 0.3 if args.method == "seedextend" else 0.1
    if args.k < 2:
        raise ValueError("--k must be >= 2")
    if args.method == "seedextend" and args.k < 3:
        raise ValueError("seedextend method requires --k >= 3")
    if not (0.0 <= args.min_overlap_frac <= 1.0):
        raise ValueError("--min_overlap_frac must be between 0 and 1")

    if args.peak_mode == "existing":
        if not args.peak_bed:
            raise ValueError("--peak_bed is required when --peak_mode existing")
        if not os.path.exists(args.peak_bed):
            raise FileNotFoundError(f"Peak BED not found: {args.peak_bed}")
    else:
        if not args.treatment:
            raise ValueError("--treatment is required when --peak_mode sicer2")
        if not args.species:
            raise ValueError("--species is required when --peak_mode sicer2")

    os.makedirs(args.output_dir, exist_ok=True)
    sample_prefix = basename_before_first_dot(args.input)

    print(f"[INFO] Reading: {args.input}", file=sys.stderr)
    df_raw = read_ref_adjust_table(args.input)
    df = clean_positions_df(df_raw, chrom_allow_regex=args.chrom_allow_regex, keep_all_centers=args.keep_all_centers)
    print(f"[INFO] Total nucleosomes after clean: {df.shape[0]}", file=sys.stderr)
    if df.empty:
        raise RuntimeError("Empty table after cleaning.")

    if args.method == "mergewindow":
        stats_df, tag, center_span_bed = run_mergewindow(df, args, sample_prefix)
    else:
        stats_df, tag, center_span_bed = run_seedextend(df, args, sample_prefix)

    raw_dir = os.path.join(args.output_dir, "02_array_raw")
    raw_bed, raw_stats = write_outputs(stats_df, raw_dir, tag, center_span_bed=center_span_bed)

    sicer_cmd = ""
    if args.peak_mode == "existing":
        peak_bed = args.peak_bed
    else:
        peak_dir = os.path.join(args.output_dir, "01_sicer")
        peak_bed, sicer_cmd = run_sicer2(args, peak_dir)

    final_dir = os.path.join(args.output_dir, "03_array_peak_filtered")
    filtered_bed, filtered_stats, overlap_tsv = filter_outputs_by_peaks(
        phased_bed=raw_bed,
        phased_stats=raw_stats,
        peak_bed=peak_bed,
        out_dir=final_dir,
        min_overlap_frac=args.min_overlap_frac,
    )

    meta_path = write_run_metadata(
        out_dir=args.output_dir,
        args=args,
        raw_bed=raw_bed,
        raw_stats=raw_stats,
        peak_bed=peak_bed,
        filtered_bed=filtered_bed,
        filtered_stats=filtered_stats,
        overlap_tsv=overlap_tsv,
        sicer_cmd=sicer_cmd,
    )

    print(f"[INFO] Metadata: {meta_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
