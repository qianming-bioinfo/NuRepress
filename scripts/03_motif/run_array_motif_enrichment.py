#!/usr/bin/env python3
import os
import re
import sys
import math
import gzip
import shutil
import tempfile
import subprocess
import glob
from collections import defaultdict

import numpy as np
import pandas as pd


def need_cmd(cmd: str) -> str:
    path = shutil.which(cmd)
    if path is None:
        raise RuntimeError(f"Required command not found in PATH: {cmd}")
    return path


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_csv(x: str):
    return [t.strip() for t in str(x).split(",") if t.strip()]


def parse_map_semicolon(x: str):
    out = {}
    for p in [z.strip() for z in str(x).split(";") if z.strip()]:
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def pick_first_existing(cols, candidates):
    cols_low = {c.lower(): c for c in cols}
    for cand in candidates:
        x = cols_low.get(cand.lower())
        if x is not None:
            return x
    return None


def load_chrom_sizes(path: str):
    cs = pd.read_csv(path, sep=r"\s+", header=None, names=["chr", "len"])
    return dict(zip(cs["chr"].astype(str), cs["len"].astype(int)))


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        last = merged[-1]
        if s <= last[1]:
            if e > last[1]:
                last[1] = e
        else:
            merged.append([s, e])
    return [(int(s), int(e)) for s, e in merged if e > s]


def build_interval_index(intervals_by_chr):
    starts = {}
    ends = {}
    for chrom, ints in intervals_by_chr.items():
        ints2 = merge_intervals(ints)
        intervals_by_chr[chrom] = ints2
        starts[chrom] = np.array([x[0] for x in ints2], dtype=np.int64)
        ends[chrom] = np.array([x[1] for x in ints2], dtype=np.int64)
    return starts, ends


def overlaps_any_merged(chrom, start0, end0, starts_idx, ends_idx):
    arr_s = starts_idx.get(chrom)
    arr_e = ends_idx.get(chrom)
    if arr_s is None or arr_s.size == 0:
        return False
    i = np.searchsorted(arr_s, end0, side="left") - 1
    if i >= 0 and arr_e[i] > start0:
        return True
    return False


def open_text_auto(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_tss_as_windows(annotation_path: str, tss_near_bp: int, chrlen: dict):
    intervals = defaultdict(list)
    path_lower = str(annotation_path).lower()
    is_gtf = path_lower.endswith(".gtf") or path_lower.endswith(".gtf.gz")

    with open_text_auto(annotation_path) as f:
        if not is_gtf:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fs = re.split(r"\t+|\s+", line.rstrip("\n"))
                if len(fs) < 3:
                    continue

                chrom = fs[0]
                if chrom not in chrlen:
                    continue

                try:
                    st0 = int(fs[1])
                    en0 = int(fs[2])
                except Exception:
                    continue

                strand = fs[5] if len(fs) >= 6 else "."

                if strand == "+":
                    tss0 = st0
                elif strand == "-":
                    tss0 = en0 - 1
                else:
                    if (en0 - st0) == 1:
                        tss0 = st0
                    else:
                        tss0 = int((st0 + en0 - 1) // 2)

                s = max(0, int(tss0) - int(tss_near_bp))
                e = min(chrlen[chrom], int(tss0) + int(tss_near_bp) + 1)
                if e > s:
                    intervals[chrom].append((s, e))

        else:
            transcript_rows = []
            gene_rows = []

            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue

                fs = line.rstrip("\n").split("\t")
                if len(fs) < 9:
                    continue

                chrom = fs[0]
                feature = fs[2]
                strand = fs[6]

                if chrom not in chrlen:
                    continue

                try:
                    start1 = int(fs[3])
                    end1 = int(fs[4])
                except Exception:
                    continue

                row = (chrom, feature, start1, end1, strand)

                if feature == "transcript":
                    transcript_rows.append(row)
                elif feature == "gene":
                    gene_rows.append(row)

            use_rows = transcript_rows if len(transcript_rows) > 0 else gene_rows

            for chrom, feature, start1, end1, strand in use_rows:
                if strand == "+":
                    tss0 = start1 - 1
                elif strand == "-":
                    tss0 = end1 - 1
                else:
                    tss0 = int(((start1 + end1) // 2) - 1)

                s = max(0, int(tss0) - int(tss_near_bp))
                e = min(chrlen[chrom], int(tss0) + int(tss_near_bp) + 1)
                if e > s:
                    intervals[chrom].append((s, e))

    starts, ends = build_interval_index(intervals)
    return intervals, starts, ends


def read_danpos_positions(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            low_memory=False,
            usecols=["chr", "smt_pos"],
            dtype={"chr": str},
        )
    except Exception:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            comment="#",
            low_memory=False,
            usecols=["chr", "smt_pos"],
            dtype={"chr": str},
        )
    df["smt_pos"] = pd.to_numeric(df["smt_pos"], errors="coerce")
    df = df.dropna().copy()
    df["smt_pos"] = df["smt_pos"].astype(np.int64)
    df["chr"] = df["chr"].astype(str)
    return df


def smtpos_to_dyad0(smt_pos: np.ndarray, dpos_coord: str) -> np.ndarray:
    x = (dpos_coord or "").strip().lower()
    if x in ("1based", "1based_closed", "onebased"):
        return smt_pos.astype(np.int64) - 1
    if x in ("0based", "0based_halfopen", "zerobased"):
        return smt_pos.astype(np.int64)
    raise RuntimeError(f"Unknown DPOS_COORD={dpos_coord}")


def build_mask_intervals(dfp: pd.DataFrame, core_half_bp: int, chrlen: dict, dpos_coord: str):
    dyad0 = smtpos_to_dyad0(dfp["smt_pos"].to_numpy(dtype=np.int64), dpos_coord)
    chroms = dfp["chr"].to_numpy(dtype=str)
    intervals = defaultdict(list)
    for c, d in zip(chroms, dyad0):
        if c not in chrlen:
            continue
        s = max(0, int(d) - int(core_half_bp))
        e = min(chrlen[c], int(d) + int(core_half_bp) + 1)
        if e > s:
            intervals[c].append((s, e))
    starts, ends = build_interval_index(intervals)
    return intervals, starts, ends


def build_chr_to_dyads(dfp: pd.DataFrame, chrlen: dict, dpos_coord: str):
    dfp2 = dfp[dfp["chr"].isin(chrlen.keys())].copy()
    dyad0 = smtpos_to_dyad0(dfp2["smt_pos"].to_numpy(dtype=np.int64), dpos_coord)
    chroms = dfp2["chr"].to_numpy(dtype=str)
    out = {}
    for c in np.unique(chroms):
        arr = dyad0[chroms == c]
        arr = np.sort(arr.astype(np.int64))
        out[c] = arr
    return out


def subtract_interval_by_merged(span_s, span_e, merged_intervals):
    if span_e <= span_s:
        return []
    out = []
    cur = span_s
    for s, e in merged_intervals:
        if e <= cur:
            continue
        if s >= span_e:
            break
        if s > cur:
            out.append((cur, min(s, span_e)))
        cur = max(cur, e)
        if cur >= span_e:
            break
    if cur < span_e:
        out.append((cur, span_e))
    return [(int(s), int(e)) for s, e in out if e > s]


def linker_segments_for_region(
    chrom,
    start0,
    end1,
    region_id,
    chr_to_dyads,
    mask_intervals_by_chr,
    require_two_dyads=True,
    min_linker_len=20,
    pad_bp=0,
    chrlen=None,
):
    dy = chr_to_dyads.get(chrom)
    if dy is None or dy.size == 0:
        return []

    left_idx = int(np.searchsorted(dy, int(start0), side="left"))
    right_idx = int(np.searchsorted(dy, int(end1), side="left"))
    n_dyads = right_idx - left_idx
    if n_dyads <= 0:
        return []
    if require_two_dyads and n_dyads < 2:
        return []

    first = int(dy[left_idx])
    last = int(dy[right_idx - 1])
    span_s = first
    span_e = last + 1
    if span_e <= span_s:
        return []

    masks = mask_intervals_by_chr.get(chrom, [])
    raw = subtract_interval_by_merged(span_s, span_e, masks)

    if int(pad_bp) > 0:
        padded = []
        L = chrlen[chrom]
        for s, e in raw:
            s2 = max(0, s - int(pad_bp))
            e2 = min(L, e + int(pad_bp))
            if e2 > s2:
                padded.append((s2, e2))
        padded = merge_intervals(padded)
        final = []
        for s, e in padded:
            final.extend(subtract_interval_by_merged(s, e, masks))
    else:
        final = raw

    out = []
    seg_idx = 0
    for s, e in final:
        if (e - s) >= int(min_linker_len):
            seg_idx += 1
            out.append((chrom, int(s), int(e), f"{region_id}|linker{seg_idx:04d}"))
    return out


def nucleosome_segments_for_region(
    chrom,
    start0,
    end1,
    region_id,
    chr_to_dyads,
    core_half_bp=70,
    pad_bp=0,
    chrlen=None,
):
    dy = chr_to_dyads.get(chrom)
    if dy is None or dy.size == 0:
        return []

    left_idx = int(np.searchsorted(dy, int(start0), side="left"))
    right_idx = int(np.searchsorted(dy, int(end1), side="left"))
    n_dyads = right_idx - left_idx
    if n_dyads <= 0:
        return []

    intervals = []
    L = chrlen[chrom]
    for d in dy[left_idx:right_idx]:
        s = max(0, int(d) - int(core_half_bp) - int(pad_bp))
        e = min(L, int(d) + int(core_half_bp) + 1 + int(pad_bp))

        if s < int(start0):
            s = int(start0)
        if e > int(end1):
            e = int(end1)

        if e > s:
            intervals.append((s, e))

    intervals = merge_intervals(intervals)

    out = []
    seg_idx = 0
    for s, e in intervals:
        if e > s:
            seg_idx += 1
            out.append((chrom, int(s), int(e), f"{region_id}|nuc{seg_idx:04d}"))
    return out


def internal_segments_for_region(
    chrom,
    start0,
    end1,
    region_id,
    chr_to_dyads,
    mask_intervals_by_chr,
    internal_region_mode="linker",
    require_two_dyads=True,
    min_linker_len=20,
    linker_pad_bp=0,
    nuc_core_half_bp=70,
    nuc_pad_bp=0,
    chrlen=None,
):
    mode = str(internal_region_mode).strip().lower()
    if mode not in ("linker", "nucleosome", "both"):
        raise RuntimeError(f"Unknown --internal_region_mode: {internal_region_mode}")

    out = []

    if mode in ("linker", "both"):
        out.extend(
            linker_segments_for_region(
                chrom=chrom,
                start0=start0,
                end1=end1,
                region_id=region_id,
                chr_to_dyads=chr_to_dyads,
                mask_intervals_by_chr=mask_intervals_by_chr,
                require_two_dyads=require_two_dyads,
                min_linker_len=min_linker_len,
                pad_bp=linker_pad_bp,
                chrlen=chrlen,
            )
        )

    if mode in ("nucleosome", "both"):
        out.extend(
            nucleosome_segments_for_region(
                chrom=chrom,
                start0=start0,
                end1=end1,
                region_id=region_id,
                chr_to_dyads=chr_to_dyads,
                core_half_bp=nuc_core_half_bp,
                pad_bp=nuc_pad_bp,
                chrlen=chrlen,
            )
        )

    return out


def sample_background_regions(bg_df: pd.DataFrame, tar_df: pd.DataFrame, n_need: int, match_location: bool, seed=1):
    rng = np.random.default_rng(seed)
    if bg_df.shape[0] == 0 or int(n_need) <= 0:
        return bg_df.iloc[0:0].copy()

    if (not match_location) or ("location" not in bg_df.columns) or ("location" not in tar_df.columns):
        if bg_df.shape[0] <= n_need:
            return bg_df.copy()
        take = rng.choice(bg_df.index.to_numpy(), size=n_need, replace=False)
        return bg_df.loc[take].copy()

    loc_counts = tar_df["location"].value_counts(dropna=False)
    loc_frac = loc_counts / loc_counts.sum()

    picked = []
    used = set()

    for loc, frac in loc_frac.items():
        nk = int(round(n_need * float(frac)))
        if nk <= 0:
            continue
        pool = bg_df[bg_df["location"] == loc]
        if pool.shape[0] == 0:
            continue
        nk2 = min(nk, pool.shape[0])
        take = rng.choice(pool.index.to_numpy(), size=nk2, replace=False)
        picked.append(bg_df.loc[take])
        used.update(take.tolist())

    out = pd.concat(picked, axis=0) if picked else bg_df.iloc[0:0].copy()
    if out.shape[0] < n_need:
        remain = n_need - out.shape[0]
        rest = bg_df.drop(index=list(used), errors="ignore")
        if rest.shape[0] > 0:
            take_n = min(remain, rest.shape[0])
            take = rng.choice(rest.index.to_numpy(), size=take_n, replace=False)
            out = pd.concat([out, rest.loc[take]], axis=0)

    return out.copy()


def write_bed(df_bed: pd.DataFrame, out_bed: str):
    if df_bed.shape[0] == 0:
        open(out_bed, "w").close()
        return
    df_bed = df_bed.sort_values(["chr", "start0", "end1", "id"]).copy()
    df_bed.to_csv(out_bed, sep="\t", index=False, header=False)


def clip_1based_inclusive(chr_arr, st1, en1, chrlen):
    out_st1 = np.array(st1, dtype=np.int64, copy=True)
    out_en1 = np.array(en1, dtype=np.int64, copy=True)
    ok = np.ones(len(chr_arr), dtype=bool)
    for i, c in enumerate(chr_arr):
        L = chrlen.get(c)
        if L is None:
            ok[i] = False
            continue
        s = out_st1[i]
        e = out_en1[i]
        if s < 1:
            s = 1
        if e > L:
            e = L
        if s > e:
            ok[i] = False
            continue
        out_st1[i] = s
        out_en1[i] = e
    return ok, out_st1, out_en1


def make_edge_windows(start0, end1, flank_bp, edge_bp):
    summit_L = start0 + 1
    summit_R = end1

    outL_st1 = summit_L - flank_bp
    outL_en1 = summit_L - 1

    inL_st1 = summit_L
    inL_en1 = summit_L + edge_bp - 1

    inR_st1 = summit_R - edge_bp + 1
    inR_en1 = summit_R

    outR_st1 = summit_R + 1
    outR_en1 = summit_R + flank_bp

    return (
        summit_L,
        summit_R,
        outL_st1,
        outL_en1,
        inL_st1,
        inL_en1,
        inR_st1,
        inR_en1,
        outR_st1,
        outR_en1,
    )


def make_motif_window_1based(summit_1based, side, outside_bp, inside_bp):
    if side == "L":
        st1 = summit_1based - outside_bp
        en1 = summit_1based + inside_bp
    else:
        st1 = summit_1based - inside_bp
        en1 = summit_1based + outside_bp
    return st1, en1


def parse_single_rel_range(s: str):
    s = (s or "").strip()
    if s == "":
        return None
    for sep in [":", ","]:
        if sep in s:
            a, b = s.split(sep, 1)
            try:
                return int(a.strip()), int(b.strip())
            except Exception:
                raise RuntimeError(f"Invalid relative range '{s}'. Use like '+275:+475' or '-200:+300'.")
    raise RuntimeError(f"Invalid relative range '{s}'. Use like '+275:+475' or '-200:+300'.")


def parse_rel_range_list(x: str):
    s = str(x or "").strip()
    if s == "":
        return []
    parts = [p.strip() for p in s.split(";") if p.strip()]
    out = []
    for i, p in enumerate(parts, start=1):
        rr = parse_single_rel_range(p)
        if rr is None:
            continue
        out.append((i, rr[0], rr[1], p))
    return out


def shift_rel_coord_away_from_zero(x: int, offset_bp: int, partner: int | None = None):
    x = int(x)
    off = abs(int(offset_bp))
    if off == 0:
        return x
    if x < 0:
        return x - off
    if x > 0:
        return x + off
    if partner is None:
        return 0
    partner = int(partner)
    if partner < 0:
        return -off
    if partner > 0:
        return off
    return 0


def shift_rel_range(rel_a: int, rel_b: int, offset_bp: int):
    a = shift_rel_coord_away_from_zero(rel_a, offset_bp, partner=rel_b)
    b = shift_rel_coord_away_from_zero(rel_b, offset_bp, partner=rel_a)
    return int(a), int(b)


def make_motif_window_1based_relative(summit_1based: int, side: str, rel_a: int, rel_b: int):
    u1, u2 = rel_a, rel_b
    if u1 > u2:
        u1, u2 = u2, u1

    if side == "L":
        p1 = summit_1based - u1
        p2 = summit_1based - u2
    else:
        p1 = summit_1based + u1
        p2 = summit_1based + u2

    st1 = min(p1, p2)
    en1 = max(p1, p2)
    return st1, en1


def build_equal_distance_bins(start_bp: int, end_bp: int, num_bins: int):
    start_bp = int(start_bp)
    end_bp = int(end_bp)
    num_bins = int(num_bins)

    if num_bins < 1:
        raise RuntimeError("--outside_num_bins must be >= 1")
    if start_bp <= 0:
        raise RuntimeError("--outside_start_bp must be > 0")
    if end_bp <= start_bp:
        raise RuntimeError("--outside_end_bp must be > --outside_start_bp")

    cuts = np.linspace(start_bp, end_bp, num_bins + 1)
    cuts = np.rint(cuts).astype(int)
    cuts[0] = start_bp
    cuts[-1] = end_bp

    for i in range(1, len(cuts)):
        if cuts[i] <= cuts[i - 1]:
            cuts[i] = cuts[i - 1] + 1

    bins = []
    for i in range(num_bins):
        a = int(cuts[i])
        b = int(cuts[i + 1])
        if b <= a:
            raise RuntimeError(f"Invalid outside-bin construction: bin{i+1} has non-positive width ({a}, {b})")
        label = f"bin{i+1}_{a}to{b}"
        bins.append((i + 1, a, b, label))
    return bins


def build_equal_distance_bins_rel(start_rel: int, end_rel: int, num_bins: int, err_prefix: str = "range"):
    start_rel = int(start_rel)
    end_rel = int(end_rel)
    num_bins = int(num_bins)

    if num_bins < 1:
        raise RuntimeError(f"--{err_prefix}_num_bins must be >= 1")
    if end_rel <= start_rel:
        raise RuntimeError(f"--{err_prefix}_end_rel must be > --{err_prefix}_start_rel")

    cuts = np.linspace(start_rel, end_rel, num_bins + 1)
    cuts = np.rint(cuts).astype(int)
    cuts[0] = start_rel
    cuts[-1] = end_rel

    for i in range(1, len(cuts)):
        if cuts[i] <= cuts[i - 1]:
            cuts[i] = cuts[i - 1] + 1

    bins = []
    for i in range(num_bins):
        a = int(cuts[i])
        b = int(cuts[i + 1])
        if b <= a:
            raise RuntimeError(f"Invalid {err_prefix}-bin construction: bin{i+1} has non-positive width ({a}, {b})")
        label = f"bin{i+1}_{a}to{b}"
        bins.append((i + 1, a, b, label))
    return bins


def build_windows_from_best_side(meta_df: pd.DataFrame, rel_a: int, rel_b: int, chrlen: dict, clip_to_array: bool = False):
    if meta_df.shape[0] == 0:
        return meta_df.iloc[0:0].copy(), pd.DataFrame(columns=["chr", "start0", "end1", "id"])

    arr_chr = meta_df["chr"].astype(str).to_numpy()
    arr_side = meta_df["best_side"].astype(str).to_numpy()
    arr_summit = meta_df["best_summit_1based"].to_numpy(dtype=np.int64)

    st1 = np.zeros(arr_summit.shape[0], dtype=np.int64)
    en1 = np.zeros(arr_summit.shape[0], dtype=np.int64)
    for i in range(arr_summit.size):
        s1, e1 = make_motif_window_1based_relative(int(arr_summit[i]), arr_side[i], int(rel_a), int(rel_b))
        st1[i] = s1
        en1[i] = e1

    ok = np.ones(arr_summit.shape[0], dtype=bool)
    st1c = np.array(st1, dtype=np.int64, copy=True)
    en1c = np.array(en1, dtype=np.int64, copy=True)

    if clip_to_array:
        if ("start" not in meta_df.columns) or ("end" not in meta_df.columns):
            raise RuntimeError("meta_df must contain start/end columns when clip_to_array=True")
        arr_start0 = meta_df["start"].to_numpy(dtype=np.int64)
        arr_end1 = meta_df["end"].to_numpy(dtype=np.int64)
    else:
        arr_start0 = None
        arr_end1 = None

    for i, c in enumerate(arr_chr):
        L = chrlen.get(c)
        if L is None:
            ok[i] = False
            continue

        s = int(st1c[i])
        e = int(en1c[i])

        if s < 1:
            s = 1
        if e > L:
            e = L

        if clip_to_array:
            arr_lo = int(arr_start0[i]) + 1
            arr_hi = int(arr_end1[i])
            if s < arr_lo:
                s = arr_lo
            if e > arr_hi:
                e = arr_hi

        if s > e:
            ok[i] = False
            continue

        st1c[i] = s
        en1c[i] = e

    keep = meta_df.iloc[np.where(ok)[0]].copy()
    chr2 = arr_chr[ok]
    st0 = (st1c[ok] - 1).astype(np.int64)
    en1o = en1c[ok].astype(np.int64)
    rid = (keep["region_id"].astype(str) + "|" + keep["best_side"].astype(str)).to_list()

    bed_df = pd.DataFrame({"chr": chr2, "start0": st0, "end1": en1o, "id": rid})
    return keep, bed_df


class BWMeanEngine:
    def __init__(self, bw_path: str):
        self.bw_path = bw_path
        self.backend = None
        self.bw = None

        try:
            import pyBigWig  # type: ignore

            self.backend = "pybigwig"
            self.bw = pyBigWig.open(bw_path)
            return
        except Exception:
            pass

        if shutil.which("bigWigAverageOverBed") is not None:
            self.backend = "ucsc"
            return

        raise RuntimeError("No backend to read bigWig. Install pyBigWig or provide bigWigAverageOverBed in PATH.")

    def close(self):
        if self.backend == "pybigwig" and self.bw is not None:
            try:
                self.bw.close()
            except Exception:
                pass

    def mean_intervals(self, chr_arr, st1_arr, en1_arr, chrlen):
        n = len(chr_arr)
        out = np.full(n, np.nan, dtype=np.float64)

        ok, st1c, en1c = clip_1based_inclusive(chr_arr, st1_arr, en1_arr, chrlen)
        idxs = np.where(ok)[0]
        if idxs.size == 0:
            return out

        if self.backend == "pybigwig":
            for i in idxs:
                c = chr_arr[i]
                st0 = int(st1c[i] - 1)
                en0 = int(en1c[i])
                if en0 <= st0:
                    continue
                try:
                    v = self.bw.stats(c, st0, en0, type="mean")[0]
                    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
                        out[i] = 0.0
                    else:
                        out[i] = float(v)
                except Exception:
                    out[i] = np.nan
            return out

        with tempfile.TemporaryDirectory() as td:
            bed = os.path.join(td, "q.bed")
            outtab = os.path.join(td, "o.tab")
            with open(bed, "w") as f:
                for i in idxs:
                    c = chr_arr[i]
                    st0 = int(st1c[i] - 1)
                    en0 = int(en1c[i])
                    if en0 <= st0:
                        continue
                    f.write(f"{c}\t{st0}\t{en0}\t{i}\n")

            cmd = ["bigWigAverageOverBed", self.bw_path, bed, outtab]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            tab = pd.read_csv(outtab, sep="\t", header=None)
            for _, row in tab.iterrows():
                idx = int(row.iloc[0])
                mean0 = row.iloc[4] if len(row) > 4 else np.nan
                mean = row.iloc[5] if len(row) > 5 else np.nan
                v = mean0
                try:
                    v = float(v)
                    if not np.isfinite(v):
                        v = float(mean)
                except Exception:
                    try:
                        v = float(mean)
                    except Exception:
                        v = np.nan
                out[idx] = v
        return out


def score_edges(dir_mode, m_outL, m_inL, m_inR, m_outR):
    if dir_mode == "inside_low":
        sL = m_outL - m_inL
        sR = m_outR - m_inR
    elif dir_mode == "inside_high":
        sL = m_inL - m_outL
        sR = m_inR - m_outR
    else:
        sL = np.abs(m_inL - m_outL)
        sR = np.abs(m_inR - m_outR)
    return sL, sR


def preferred_sample_order(sample_values):
    pref = ["HPNE", "PANC-1", "Capan-1"]
    sset = set([str(x) for x in sample_values])
    return list(dict.fromkeys([s for s in pref if s in sset] + sorted(sset - set(pref))))


def parse_args():
    import argparse

    ap = argparse.ArgumentParser(
        description=(
            "Run motif enrichment on phased-array subtypes by directly consuming "
            "cluster_array_subtype.R outputs. The script can automatically read "
            "best_k_selection.tsv, locate the best-k cluster assignment table, build "
            "target/background BEDs, and optionally run HOMER."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io = ap.add_argument_group("Input / output")
    io.add_argument("--cluster_run_dir", type=str, default=None,
                    help="Output directory from cluster_array_subtype.R. If provided, the script reads best_k_selection.tsv and auto-locates the best-k cluster assignment TSV.")
    io.add_argument("--cluster_tsv", type=str, default=None,
                    help="Cluster assignment TSV. If set, it overrides --cluster_run_dir auto-detection.")
    io.add_argument("-O", "--out_dir", required=True,
                    help="Output directory for BEDs, metadata, and HOMER results.")
    io.add_argument("--annotation_file", default=None,
                    help="Annotation file for optional TSS filtering. Supports BED/BED.gz/GTF/GTF.gz. Required only when --tss_near_bp >= 0.")
    io.add_argument("--chrom_sizes", required=True,
                    help="Chromosome sizes file with two columns: chrom and length.")

    cl = ap.add_argument_group("Cluster selection")
    cl.add_argument("--target_clusters", default="ALL",
                    help="Comma-separated cluster labels to analyze, e.g. C1,C2. Use ALL to analyze all detected clusters.")
    cl.add_argument("--cluster_col_candidates", default="k2_cluster,cluster,pattern_cluster",
                    help="Comma-separated candidate cluster-column names. The first existing column will be used.")
    cl.add_argument("--fit_space_preference", default="within_sample_scaled",
                    help="Preferred fit-space label when multiple best-k cluster assignment TSVs exist under cluster_run_dir.")

    mode = ap.add_argument_group("Motif-region mode")
    mode.add_argument("--mode", choices=["both", "internal_linker", "edge", "edge_outside", "edge_inside", "edge_both", "all"], default="edge",
                      help=(
                          "How target motif regions are built from phased arrays. "
                          "edge reproduces the older bash workflow: one explicit relative range (or one semicolon-separated list of ranges) corresponds to one independent HOMER job per sample-cluster."
                      ))
    mode.add_argument("--tss_near_bp", type=int, default=-1,
                      help="Keep arrays overlapping a TSS-centered window of +/- this size. Use -1 to disable TSS filtering entirely. Use 0 to keep only arrays that overlap the exact TSS base.")
    mode.add_argument("--min_n_target_arrays", type=int, default=30,
                      help="General minimum target-array threshold. In edge mode, the actual skip check follows the older bash behavior and is applied after score/q_keep filtering.")
    mode.add_argument("--bg_multiple", type=float, default=2.0,
                      help="Number of background arrays relative to the number of target arrays.")
    mode.add_argument("--bg_match_location", action="store_true", default=True,
                      help="Sample background arrays matched by the 'location' column when available.")
    mode.add_argument("--no_bg_match_location", dest="bg_match_location", action="store_false",
                      help="Do not match background arrays by the 'location' column.")
    mode.add_argument("--homer_auto_bg", dest="homer_auto_bg", action="store_true", default=True,
                      help="Use HOMER internal background selection instead of writing explicit background BEDs.")
    mode.add_argument("--no_homer_auto_bg", dest="homer_auto_bg", action="store_false",
                      help="Write explicit background BEDs and pass -bg to HOMER.")
    mode.add_argument("--bg_sampling_strategy", choices=["per_bin", "shared"], default="per_bin",
                      help="When explicit backgrounds are used in edge modes, sample backgrounds independently for each bin or once per sample-cluster and share across bins.")

    il = ap.add_argument_group("internal_linker mode")
    il.add_argument("--positions_map", type=str, default=None,
                    help="Sample-to-position-file map in the form 'sample1=/path/a.positions.xls;sample2=/path/b.positions.xls'. Preferred over --dpos_map when both are provided.")
    il.add_argument("--dpos_map", type=str, default=None,
                    help="Backward-compatible alias for position-file map in the form 'sample1=/path/a.xls;sample2=/path/b.xls'.")
    il.add_argument("--dpos_coord", choices=["1based_closed", "0based_halfopen"], default="1based_closed",
                    help="Coordinate convention for DANPOS smt_pos.")
    il.add_argument("--nuc_core_half_bp", type=int, default=70,
                    help="Half-width of the nucleosome core around each dyad.")
    il.add_argument("--min_linker_len", type=int, default=20,
                    help="Minimum linker length to keep.")
    il.add_argument("--internal_pad_bp", type=int, default=0,
                    help="Pad each internal linker segment on both sides by this many bp.")
    il.add_argument("--internal_region_mode", choices=["linker", "nucleosome", "both"], default="linker",
                    help="Internal-region composition used in internal_linker mode: linker only, nucleosome only, or both.")
    il.add_argument("--internal_nuc_pad_bp", type=int, default=0,
                    help="Pad each internal nucleosome segment on both sides by this many bp.")
    il.add_argument("--require_two_dyads_for_linker", action="store_true", default=True,
                    help="Require at least two dyads within an array to define internal linkers. This setting only affects linker extraction.")
    il.add_argument("--allow_single_dyad_linker", dest="require_two_dyads_for_linker", action="store_false",
                    help="Allow linker extraction from arrays with a single dyad. This setting only affects linker extraction.")

    eo = ap.add_argument_group("edge-based modes")
    eo.add_argument("--ins_bw_map", type=str, default=None,
                    help="Sample-to-insertion-bigWig map in the form 'sample1=/path/a.bw;sample2=/path/b.bw'. Required for all edge modes.")
    eo.add_argument("--cluster_dir_default", type=str, default="neutral",
                    choices=["neutral", "inside_low", "inside_high"],
                    help="Default edge-scoring direction for clusters not listed in --cluster_dir_map.")
    eo.add_argument("--cluster_dir_map", type=str, default="",
                    help="Optional cluster-to-direction map, e.g. 'C1=inside_low;C2=neutral'.")
    eo.add_argument("--require_positive_score", action="store_true", default=True,
                    help="Require positive best edge score for inside_low/inside_high modes.")
    eo.add_argument("--allow_nonpositive_score", dest="require_positive_score", action="store_false",
                    help="Keep windows even if the edge score is <= 0.")
    eo.add_argument("--score_flank_bp", type=int, default=800,
                    help="Flanking distance used to score left/right array edges.")
    eo.add_argument("--score_edge_bp", type=int, default=200,
                    help="Edge window size used to compare inside/outside ATAC signal.")

    eo.add_argument("--outside_start_bp", type=int, default=75,
                    help="Start distance of the outside boundary interval.")
    eo.add_argument("--outside_end_bp", type=int, default=675,
                    help="End distance of the outside boundary interval.")
    eo.add_argument("--outside_num_bins", type=int, default=3,
                    help="Number of equal-width bins to split the outside interval into. Each bin will be analyzed separately by HOMER.")

    eo.add_argument("--strict_outside_only", action="store_true", default=True,
                    help="Require outside intervals to start at positive distance from the boundary summit.")
    eo.add_argument("--allow_inside_overlap", dest="strict_outside_only", action="store_false",
                    help="Allow intervals that overlap the boundary or interior.")

    eo.add_argument("--rel_range", type=str, default=None,
                    help="When --mode edge, one or more semicolon-separated ranges, e.g. '+400:+600' or '+400:+600;+600:+800'. Each range is run as an independent HOMER job per sample-cluster. Ignored by other edge modes.")
    eo.add_argument("--range_offset_bp", type=int, default=70,
                    help="Default offset added to user-specified relative coordinates in edge mode. Set 0 to disable.")
    eo.add_argument("--motif_outside_bp", type=int, default=375,
                    help="Outside window size used in edge mode when --rel_range is empty.")
    eo.add_argument("--motif_inside_bp", type=int, default=275,
                    help="Inside window size used in edge mode when --rel_range is empty.")
    eo.add_argument("--q_keep", type=float, default=0.95,
                    help="Keep the top fraction of target arrays by best edge score.")

    ei = ap.add_argument_group("edge_inside mode")
    ei.add_argument("--inside_start_rel", type=int, default=-475,
                    help="Start relative coordinate of the inside-boundary interval.")
    ei.add_argument("--inside_end_rel", type=int, default=-75,
                    help="End relative coordinate of the inside-boundary interval. Must be greater than --inside_start_rel.")
    ei.add_argument("--inside_num_bins", type=int, default=2,
                    help="Number of equal-width bins to split the inside interval into. Each bin will be analyzed separately by HOMER.")
    ei.add_argument("--inside_rel_range", type=str, default=None,
                    help="Deprecated compatibility parameter. Ignored when inside bins are used.")
    ei.add_argument("--inside_clip_to_array", action="store_true", default=False,
                    help="Clip edge_inside windows to the array body.")

    hm = ap.add_argument_group("HOMER")
    hm.add_argument("--genome", required=True,
                    help="Genome argument passed to findMotifsGenome.pl, e.g. hg19, hg38, mm10.")
    hm.add_argument("--homer_exec", default="findMotifsGenome.pl",
                    help="Path to findMotifsGenome.pl.")
    hm.add_argument("--threads", type=int, default=8,
                    help="Threads passed to HOMER with -p.")
    hm.add_argument("--homer_len", default="8,10,12,15,20",
                    help="Motif lengths passed to HOMER via -len.")
    hm.add_argument("--run_known_only", action="store_true", default=True,
                    help="Pass -nomotif to HOMER and only run known-motif enrichment.")
    hm.add_argument("--run_denovo_too", dest="run_known_only", action="store_false",
                    help="Do not pass -nomotif; run de novo motif discovery too.")
    hm.add_argument("--build_only", action="store_true",
                    help="Only build target/background BEDs and do not run HOMER.")
    hm.add_argument("--extra_homer_args", default="",
                    help="Additional raw arguments appended to each HOMER command line.")

    return ap.parse_args()


def split_shell_like(s: str):
    import shlex
    if s is None:
        return []
    s = str(s).strip()
    return [] if s == "" else shlex.split(s)


def infer_cluster_tsv(cluster_run_dir: str, fit_space_preference: str = "within_sample_scaled"):
    best_k_path = os.path.join(cluster_run_dir, "best_k_selection.tsv")
    if not os.path.exists(best_k_path):
        raise RuntimeError(f"Missing best_k_selection.tsv under cluster_run_dir: {cluster_run_dir}")
    bk = pd.read_csv(best_k_path, sep="\t", low_memory=False)
    if "best_k" not in bk.columns or bk.shape[0] == 0:
        raise RuntimeError(f"best_k_selection.tsv does not contain a valid best_k column: {best_k_path}")
    best_k = int(pd.to_numeric(bk.loc[0, "best_k"], errors="coerce"))
    if not np.isfinite(best_k):
        raise RuntimeError(f"Invalid best_k value in {best_k_path}")

    cluster_dir = os.path.join(cluster_run_dir, "cluster_tables")
    if not os.path.isdir(cluster_dir):
        raise RuntimeError(f"Missing cluster_tables directory under cluster_run_dir: {cluster_run_dir}")

    patt = os.path.join(cluster_dir, f"k{best_k}_cluster_assignment.z1z2.*.tsv")
    cands = sorted(glob.glob(patt))
    if not cands:
        exact = os.path.join(cluster_dir, f"k{best_k}_cluster_assignment.z1z2.{fit_space_preference}.tsv")
        raise RuntimeError(
            f"Could not find a best-k cluster assignment TSV. Tried pattern: {patt}; preferred exact path: {exact}"
        )

    preferred = [p for p in cands if fit_space_preference in os.path.basename(p)]
    if preferred:
        cluster_tsv = preferred[0]
    else:
        cluster_tsv = cands[0]
    return cluster_tsv, best_k_path, best_k


def natural_key_cluster(x: str):
    m = re.match(r"^C(\d+)$", str(x))
    if m:
        return (0, int(m.group(1)), str(x))
    return (1, str(x))


def run_homer_from_summary(args, summary_df, outdir):
    need_cmd(args.homer_exec)
    extra_homer_args = split_shell_like(args.extra_homer_args)
    total = 0
    ok = 0

    for row in summary_df.itertuples(index=False):
        if str(getattr(row, "note", "")) not in ("ok", "ok_auto_bg"):
            continue
        tar = getattr(row, "target_bed", None)
        if tar is None or not isinstance(tar, str) or not os.path.exists(tar) or os.path.getsize(tar) == 0:
            continue

        total += 1
        tag = os.path.basename(tar)
        if tag.startswith("target."):
            tag = tag[len("target."):]
        if tag.endswith(".bed"):
            tag = tag[:-4]

        out = os.path.join(outdir, "homer_out", tag)
        os.makedirs(out, exist_ok=True)
        log1 = os.path.join(outdir, "homer_logs", f"{tag}.out.txt")
        log2 = os.path.join(outdir, "homer_logs", f"{tag}.err.txt")

        cmd = [
            args.homer_exec,
            tar,
            args.genome,
            out,
            "-size",
            "given",
            "-p",
            str(args.threads),
            "-len",
            str(args.homer_len),
        ]
        if args.run_known_only:
            cmd.append("-nomotif")
        if not args.homer_auto_bg:
            bg = getattr(row, "bg_bed", None)
            if bg is None or not isinstance(bg, str) or (not os.path.exists(bg)) or os.path.getsize(bg) == 0:
                with open(log2, "a") as fh:
                    fh.write("[SKIP] Missing or empty background BED for explicit-background HOMER run.\n")
                continue
            cmd.extend(["-bg", bg])
        cmd.extend(extra_homer_args)

        eprint(f">>> HOMER ({ok+1}/{total}) tag={tag} threads={args.threads} auto_bg={int(args.homer_auto_bg)}")
        with open(log1, "w") as fo, open(log2, "w") as fe:
            fe.write("COMMAND:\n" + " ".join(cmd) + "\n")
            subprocess.run(cmd, stdout=fo, stderr=fe, check=True)
        ok += 1

    return {"n_homer_jobs_total": int(total), "n_homer_jobs_completed": int(ok)}


def build_edge_signal_meta(df_source: pd.DataFrame, samples, ins_map, chrlen, score_flank_bp, score_edge_bp):
    meta_by_sample = {}
    for sm in samples:
        if sm not in ins_map:
            raise RuntimeError(f"Missing insertion bigWig for sample={sm} in --ins_bw_map")
        bw_path = ins_map[sm]
        if not os.path.exists(bw_path):
            raise RuntimeError(f"Insertion bigWig not found: {bw_path}")

        eprint(f"\n>>> Sample={sm} | mode=edge_based")
        df_sm = df_source[df_source["sample"] == sm].copy()
        if df_sm.shape[0] == 0:
            continue

        chr_arr = df_sm["chr"].to_numpy(dtype=str)
        start0 = df_sm["start"].to_numpy(dtype=np.int64)
        end1 = df_sm["end"].to_numpy(dtype=np.int64)
        (
            summit_L,
            summit_R,
            outL_st1,
            outL_en1,
            inL_st1,
            inL_en1,
            inR_st1,
            inR_en1,
            outR_st1,
            outR_en1,
        ) = make_edge_windows(start0, end1, score_flank_bp, score_edge_bp)

        engine = BWMeanEngine(bw_path)
        try:
            eprint(f"  BW backend: {engine.backend} ({bw_path})")
            m_outL = engine.mean_intervals(chr_arr, outL_st1, outL_en1, chrlen)
            m_inL = engine.mean_intervals(chr_arr, inL_st1, inL_en1, chrlen)
            m_inR = engine.mean_intervals(chr_arr, inR_st1, inR_en1, chrlen)
            m_outR = engine.mean_intervals(chr_arr, outR_st1, outR_en1, chrlen)
        finally:
            engine.close()

        df_sm = df_sm.reset_index(drop=True)
        df_sm["summit_L_1based"] = summit_L
        df_sm["summit_R_1based"] = summit_R
        df_sm["m_outL"] = m_outL
        df_sm["m_inL"] = m_inL
        df_sm["m_inR"] = m_inR
        df_sm["m_outR"] = m_outR
        meta_by_sample[sm] = df_sm
    return meta_by_sample


def main():
    args = parse_args()

    if args.cluster_tsv is None and args.cluster_run_dir is None:
        raise RuntimeError("One of --cluster_tsv or --cluster_run_dir must be provided")

    if args.cluster_tsv is not None:
        cluster_tsv = args.cluster_tsv
        best_k_path = None
        best_k = None
    else:
        cluster_tsv, best_k_path, best_k = infer_cluster_tsv(args.cluster_run_dir, args.fit_space_preference)

    mode = args.mode
    do_internal = mode in ("both", "internal_linker", "all")
    do_edge_outside = mode in ("both", "edge_outside", "edge_both", "all")
    do_edge_inside = mode in ("edge_inside", "edge_both", "all")
    do_edge = mode == "edge"
    do_edge_any = do_edge_outside or do_edge_inside or do_edge

    outdir = args.out_dir
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "bed_target"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "bed_bg"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "meta"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "homer_out"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "homer_logs"), exist_ok=True)

    chrlen = load_chrom_sizes(args.chrom_sizes)

    if not os.path.exists(cluster_tsv):
        raise RuntimeError(f"Missing cluster TSV: {cluster_tsv}")
    if args.annotation_file is not None and not os.path.exists(args.annotation_file):
        raise RuntimeError(f"Missing annotation file: {args.annotation_file}")

    if args.tss_near_bp >= 0 and not args.annotation_file:
        raise RuntimeError("--annotation_file is required when --tss_near_bp >= 0")

    tss_filter_enabled = args.tss_near_bp >= 0
    tss_tag = f"nearTSS{args.tss_near_bp}" if tss_filter_enabled else "noTSSFilter"

    eprint(f">>> Reading cluster TSV: {cluster_tsv}")
    df = pd.read_csv(cluster_tsv, sep="\t", low_memory=False)

    cluster_candidates = parse_csv(args.cluster_col_candidates)
    cluster_col = pick_first_existing(df.columns, cluster_candidates)
    if cluster_col is None:
        raise RuntimeError(f"No cluster column found. Tried: {cluster_candidates}")

    need_cols = {"sample", "region_id", "chr", "start", "end", cluster_col}
    miss = [c for c in need_cols if c not in df.columns]
    if miss:
        raise RuntimeError(f"Missing columns in cluster TSV: {miss}")

    if "location" not in df.columns:
        df["location"] = "NA"

    df["sample"] = df["sample"].astype(str)
    df["region_id"] = df["region_id"].astype(str)
    df["cluster"] = df[cluster_col].astype(str).str.strip()
    df["chr"] = df["chr"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df["location"] = df["location"].astype(str)

    df = df[df["cluster"].notna()].copy()
    df = df[df["cluster"].str.len() > 0].copy()
    df = df[~df["cluster"].isin(["NA", "na", "Na", "nA", "nan", "None"])].copy()
    df = df[df["start"].notna() & df["end"].notna()].copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["chr"].isin(chrlen.keys())].copy()

    if tss_filter_enabled:
        eprint(f">>> Building TSS-near windows from: {args.annotation_file}")
        _, tss_starts_idx, tss_ends_idx = read_tss_as_windows(args.annotation_file, args.tss_near_bp, chrlen)

        near_tss = []
        for r in df.itertuples(index=False):
            near_tss.append(overlaps_any_merged(r.chr, int(r.start), int(r.end), tss_starts_idx, tss_ends_idx))
        df["near_tss"] = np.array(near_tss, dtype=bool)
        df_tss = df[df["near_tss"]].copy()
        df_tss.to_csv(os.path.join(outdir, "meta", "cluster_table.filtered_nearTSS.tsv"), sep="\t", index=False)

        eprint(f">>> Arrays total: {df.shape[0]}")
        eprint(f">>> Arrays near TSS (+/-{args.tss_near_bp}): {df_tss.shape[0]}")

        if df_tss.shape[0] == 0:
            raise RuntimeError("No arrays remain after TSS filtering")
    else:
        eprint(">>> TSS filtering: disabled (using all arrays)")
        df["near_tss"] = True
        df_tss = df.copy()
        df_tss.to_csv(os.path.join(outdir, "meta", "cluster_table.no_TSS_filter.tsv"), sep="\t", index=False)

        eprint(f">>> Arrays total: {df.shape[0]}")
        eprint(f">>> Arrays used (no TSS filter): {df_tss.shape[0]}")

    samples = preferred_sample_order(df_tss["sample"].astype(str).tolist())
    all_clusters = sorted(df_tss["cluster"].dropna().astype(str).unique().tolist(), key=natural_key_cluster)
    if str(args.target_clusters).strip().upper() == "ALL":
        target_clusters = all_clusters
    else:
        target_clusters = parse_csv(args.target_clusters)
    if len(target_clusters) == 0:
        raise RuntimeError("No target clusters selected")

    summary_rows = []
    positions_map_used_raw = None

    if do_internal:
        positions_map_used_raw = args.positions_map if str(args.positions_map or "").strip() else args.dpos_map
        if not positions_map_used_raw:
            raise RuntimeError("--positions_map or --dpos_map is required for internal_linker mode")
        dpos_map = parse_map_semicolon(positions_map_used_raw)
        if not dpos_map:
            raise RuntimeError("Parsed --positions_map/--dpos_map is empty")

        for sm in samples:
            if sm not in dpos_map:
                raise RuntimeError(f"Missing DANPOS positions file for sample={sm} in --positions_map/--dpos_map")
            dpos_path = dpos_map[sm]
            if not os.path.exists(dpos_path):
                raise RuntimeError(f"DANPOS positions file not found: {dpos_path}")

            eprint(f"\n>>> Sample={sm} | mode=internal_linker | internal_region_mode={args.internal_region_mode}")
            df_sm = df_tss[df_tss["sample"] == sm].copy()
            if df_sm.shape[0] == 0:
                continue

            dfp = read_danpos_positions(dpos_path)
            dfp = dfp[dfp["chr"].isin(chrlen.keys())].copy()
            if dfp.shape[0] == 0:
                eprint(f"[WARN] {sm}: DANPOS empty after chrom filter, skip")
                continue

            mask_intervals_by_chr, _, _ = build_mask_intervals(dfp, args.nuc_core_half_bp, chrlen, args.dpos_coord)
            chr_to_dyads = build_chr_to_dyads(dfp, chrlen, args.dpos_coord)

            for cl in target_clusters:
                df_tar = df_sm[df_sm["cluster"] == cl].copy()
                n_tar_arrays = df_tar.shape[0]
                tag = (
                    f"{sm}.{cl}.{tss_tag}.mode_internalRegion.{args.internal_region_mode}"
                    f".core{args.nuc_core_half_bp}.minLink{args.min_linker_len}"
                    f".linkPad{args.internal_pad_bp}.nucPad{args.internal_nuc_pad_bp}"
                    f".bgx{args.bg_multiple}.loc{int(args.bg_match_location)}"
                    f".minTar{args.min_n_target_arrays}.autoBG{int(args.homer_auto_bg)}"
                )
                meta_path = os.path.join(outdir, "meta", f"meta.{tag}.tsv")
                df_tar.to_csv(meta_path, sep="\t", index=False)

                if n_tar_arrays < args.min_n_target_arrays:
                    eprint(f"[SKIP] {sm} {cl}: n_target_arrays={n_tar_arrays} < {args.min_n_target_arrays}")
                    summary_rows.append({
                        "sample": sm,
                        "cluster": cl,
                        "mode": "internal_linker",
                        "internal_region_mode": args.internal_region_mode,
                        "note": "skip_min_target_arrays",
                        "n_target_arrays": int(n_tar_arrays),
                    })
                    continue

                tar_segments = []
                n_tar_arrays_with_segments = 0
                for r in df_tar.itertuples(index=False):
                    segs = internal_segments_for_region(
                        chrom=r.chr,
                        start0=int(r.start),
                        end1=int(r.end),
                        region_id=str(r.region_id),
                        chr_to_dyads=chr_to_dyads,
                        mask_intervals_by_chr=mask_intervals_by_chr,
                        internal_region_mode=args.internal_region_mode,
                        require_two_dyads=args.require_two_dyads_for_linker,
                        min_linker_len=args.min_linker_len,
                        linker_pad_bp=args.internal_pad_bp,
                        nuc_core_half_bp=args.nuc_core_half_bp,
                        nuc_pad_bp=args.internal_nuc_pad_bp,
                        chrlen=chrlen,
                    )
                    if segs:
                        n_tar_arrays_with_segments += 1
                        tar_segments.extend(segs)

                tar_bed = os.path.join(outdir, "bed_target", f"target.{tag}.bed")
                tar_bed_df = pd.DataFrame(tar_segments, columns=["chr", "start0", "end1", "id"])
                write_bed(tar_bed_df, tar_bed)
                if tar_bed_df.shape[0] == 0:
                    eprint(f"[SKIP] {sm} {cl}: no internal segments")
                    summary_rows.append({
                        "sample": sm,
                        "cluster": cl,
                        "mode": "internal_linker",
                        "internal_region_mode": args.internal_region_mode,
                        "note": "skip_no_target_segments",
                        "n_target_arrays": int(n_tar_arrays),
                        "n_target_arrays_with_segments": int(n_tar_arrays_with_segments),
                        "n_target_segments": 0,
                    })
                    continue

                if args.homer_auto_bg:
                    eprint(f"[OK] {sm} {cl}: target_arrays={n_tar_arrays}, target_segments={tar_bed_df.shape[0]} (AUTO_BG=1)")
                    summary_rows.append({
                        "sample": sm,
                        "cluster": cl,
                        "mode": "internal_linker",
                        "internal_region_mode": args.internal_region_mode,
                        "note": "ok_auto_bg",
                        "n_target_arrays": int(n_tar_arrays),
                        "n_target_arrays_with_segments": int(n_tar_arrays_with_segments),
                        "n_target_segments": int(tar_bed_df.shape[0]),
                        "target_bed": tar_bed,
                    })
                    continue

                df_bg_pool = df_sm[df_sm["cluster"] != cl].copy()
                if df_bg_pool.shape[0] == 0:
                    eprint(f"[WARN] {sm} {cl}: bg pool empty")
                    summary_rows.append({
                        "sample": sm,
                        "cluster": cl,
                        "mode": "internal_linker",
                        "internal_region_mode": args.internal_region_mode,
                        "note": "skip_bg_pool_empty",
                        "n_target_arrays": int(n_tar_arrays),
                        "n_target_segments": int(tar_bed_df.shape[0]),
                    })
                    continue

                n_bg_need = int(math.ceil(args.bg_multiple * float(n_tar_arrays)))
                df_bg_pick = sample_background_regions(df_bg_pool, df_tar, n_bg_need, args.bg_match_location, seed=1)
                bg_segments = []
                n_bg_arrays_with_segments = 0
                for r in df_bg_pick.itertuples(index=False):
                    segs = internal_segments_for_region(
                        chrom=r.chr,
                        start0=int(r.start),
                        end1=int(r.end),
                        region_id=str(r.region_id),
                        chr_to_dyads=chr_to_dyads,
                        mask_intervals_by_chr=mask_intervals_by_chr,
                        internal_region_mode=args.internal_region_mode,
                        require_two_dyads=args.require_two_dyads_for_linker,
                        min_linker_len=args.min_linker_len,
                        linker_pad_bp=args.internal_pad_bp,
                        nuc_core_half_bp=args.nuc_core_half_bp,
                        nuc_pad_bp=args.internal_nuc_pad_bp,
                        chrlen=chrlen,
                    )
                    if segs:
                        n_bg_arrays_with_segments += 1
                        bg_segments.extend(segs)

                bg_bed = os.path.join(outdir, "bed_bg", f"bg.{tag}.bed")
                bg_bed_df = pd.DataFrame(bg_segments, columns=["chr", "start0", "end1", "id"])
                write_bed(bg_bed_df, bg_bed)
                if bg_bed_df.shape[0] == 0:
                    eprint(f"[WARN] {sm} {cl}: bg segments empty")
                    summary_rows.append({
                        "sample": sm,
                        "cluster": cl,
                        "mode": "internal_linker",
                        "internal_region_mode": args.internal_region_mode,
                        "note": "skip_bg_no_segments",
                        "n_target_arrays": int(n_tar_arrays),
                        "n_target_segments": int(tar_bed_df.shape[0]),
                        "n_bg_arrays": int(df_bg_pick.shape[0]),
                        "n_bg_arrays_with_segments": int(n_bg_arrays_with_segments),
                        "n_bg_segments": 0,
                    })
                    continue

                eprint(f"[OK] {sm} {cl}: target_arrays={n_tar_arrays}, target_segments={tar_bed_df.shape[0]}, bg_arrays={df_bg_pick.shape[0]}, bg_segments={bg_bed_df.shape[0]}")
                summary_rows.append({
                    "sample": sm,
                    "cluster": cl,
                    "mode": "internal_linker",
                    "internal_region_mode": args.internal_region_mode,
                    "note": "ok",
                    "n_target_arrays": int(n_tar_arrays),
                    "n_target_arrays_with_segments": int(n_tar_arrays_with_segments),
                    "n_target_segments": int(tar_bed_df.shape[0]),
                    "n_bg_arrays": int(df_bg_pick.shape[0]),
                    "n_bg_arrays_with_segments": int(n_bg_arrays_with_segments),
                    "n_bg_segments": int(bg_bed_df.shape[0]),
                    "target_bed": tar_bed,
                    "bg_bed": bg_bed,
                })

    if do_edge_any:
        if not args.ins_bw_map:
            raise RuntimeError("--ins_bw_map is required for edge modes")
        ins_map = parse_map_semicolon(args.ins_bw_map)
        if not ins_map:
            raise RuntimeError("Parsed --ins_bw_map is empty")
        dir_map = parse_map_semicolon(args.cluster_dir_map)

        meta_by_sample = build_edge_signal_meta(
            df_source=df_tss,
            samples=samples,
            ins_map=ins_map,
            chrlen=chrlen,
            score_flank_bp=args.score_flank_bp,
            score_edge_bp=args.score_edge_bp,
        )

        if do_edge:
            rel_specs = parse_rel_range_list(args.rel_range)
            if len(rel_specs) == 0:
                rel_specs = [(1, None, None, "legacyOutsideInside")]
            eprint(f">>> edge ranges: {len(rel_specs)} | range_offset_bp={args.range_offset_bp} (away-from-zero)")

            for sm in samples:
                df_sm = meta_by_sample.get(sm)
                if df_sm is None or df_sm.shape[0] == 0:
                    continue

                for cl in target_clusters:
                    df_cl = df_sm[df_sm["cluster"] == cl].copy()
                    if df_cl.shape[0] == 0:
                        continue

                    dir_cl = dir_map.get(str(cl), args.cluster_dir_default)
                    if dir_cl not in ("inside_low", "inside_high", "neutral"):
                        raise RuntimeError(f"Invalid dir mode for cluster={cl}: {dir_cl}")

                    sL, sR = score_edges(
                        dir_cl,
                        df_cl["m_outL"].to_numpy(),
                        df_cl["m_inL"].to_numpy(),
                        df_cl["m_inR"].to_numpy(),
                        df_cl["m_outR"].to_numpy(),
                    )
                    df_cl["score_left"] = sL
                    df_cl["score_right"] = sR
                    df_cl["best_side"] = np.where(df_cl["score_left"] >= df_cl["score_right"], "L", "R")
                    df_cl["best_score"] = np.maximum(df_cl["score_left"], df_cl["score_right"])
                    df_cl["best_summit_1based"] = np.where(
                        df_cl["best_side"] == "L",
                        df_cl["summit_L_1based"],
                        df_cl["summit_R_1based"],
                    ).astype(np.int64)
                    df_cl["dir_cluster"] = dir_cl

                    meta = df_cl[[
                        "sample", "cluster", "region_id", "chr", "start", "end", "location",
                        "dir_cluster", "score_left", "score_right", "best_side", "best_score", "best_summit_1based",
                        "m_outL", "m_inL", "m_inR", "m_outR"
                    ]].copy()

                    meta_f = meta[np.isfinite(meta["best_score"])].copy()
                    if args.require_positive_score and dir_cl in ("inside_low", "inside_high"):
                        meta_f = meta_f[meta_f["best_score"] > 0].copy()

                    if meta_f.shape[0] == 0:
                        eprint(f"[SKIP] {sm} {cl} no valid regions after score filtering (dir={dir_cl})")
                        for _, rel_a_raw, rel_b_raw, _raw_text in rel_specs:
                            if rel_a_raw is None:
                                rel_str = "NA"
                            else:
                                rel_a_use, rel_b_use = shift_rel_range(rel_a_raw, rel_b_raw, args.range_offset_bp)
                                rel_str = f"{rel_a_use:+d}to{rel_b_use:+d}"
                            tag0 = (
                                f"{sm}.{cl}.dir_{dir_cl}.scoreF{args.score_flank_bp}.scoreE{args.score_edge_bp}"
                                f".motifOut{args.motif_outside_bp}.motifIn{args.motif_inside_bp}.rel{rel_str}"
                                f".q{args.q_keep}.pos{int(args.require_positive_score)}.autoBG{int(args.homer_auto_bg)}"
                            )
                            meta_path0 = os.path.join(outdir, "meta", f"meta.{tag0}.tsv")
                            meta.to_csv(meta_path0, sep="\t", index=False)
                        continue

                    thr = meta_f["best_score"].quantile(1.0 - float(args.q_keep))
                    meta_keep = meta_f[(meta_f["best_score"] >= thr) & np.isfinite(meta_f["best_score"])].copy()

                    if meta_keep.shape[0] < int(args.min_n_target_arrays):
                        eprint(f"[SKIP motif] {sm} {cl} kept={meta_keep.shape[0]} < {args.min_n_target_arrays} (dir={dir_cl})")
                        for _, rel_a_raw, rel_b_raw, _raw_text in rel_specs:
                            if rel_a_raw is None:
                                rel_str = "NA"
                            else:
                                rel_a_use, rel_b_use = shift_rel_range(rel_a_raw, rel_b_raw, args.range_offset_bp)
                                rel_str = f"{rel_a_use:+d}to{rel_b_use:+d}"
                            tag_skip = (
                                f"{sm}.{cl}.dir_{dir_cl}.scoreF{args.score_flank_bp}.scoreE{args.score_edge_bp}"
                                f".motifOut{args.motif_outside_bp}.motifIn{args.motif_inside_bp}.rel{rel_str}"
                                f".q{args.q_keep}.bgx{args.bg_multiple}.minTar{args.min_n_target_arrays}.loc{int(args.bg_match_location)}"
                                f".pos{int(args.require_positive_score)}.autoBG{int(args.homer_auto_bg)}"
                            )
                            meta_path_skip = os.path.join(outdir, "meta", f"meta.{tag_skip}.tsv")
                            meta.to_csv(meta_path_skip, sep="\t", index=False)
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge",
                                "rel_range_used": rel_str,
                                "rel_range_used_input": "" if rel_a_raw is None else f"{rel_a_raw:+d}:{rel_b_raw:+d}",
                                "note": "skip_min_target_after_qkeep",
                                "n_target_arrays": int(df_cl.shape[0]),
                                "n_target_arrays_kept": int(meta_keep.shape[0]),
                            })
                        continue

                    for _, rel_a_raw, rel_b_raw, raw_text in rel_specs:
                        if rel_a_raw is None:
                            rel_a_use = None
                            rel_b_use = None
                            rel_str = "NA"
                        else:
                            rel_a_use, rel_b_use = shift_rel_range(rel_a_raw, rel_b_raw, args.range_offset_bp)
                            rel_str = f"{rel_a_use:+d}to{rel_b_use:+d}"

                        tag = (
                            f"{sm}.{cl}.dir_{dir_cl}.scoreF{args.score_flank_bp}.scoreE{args.score_edge_bp}"
                            f".motifOut{args.motif_outside_bp}.motifIn{args.motif_inside_bp}.rel{rel_str}"
                            f".q{args.q_keep}.bgx{args.bg_multiple}.minTar{args.min_n_target_arrays}.loc{int(args.bg_match_location)}"
                            f".pos{int(args.require_positive_score)}.autoBG{int(args.homer_auto_bg)}"
                        )
                        meta_path = os.path.join(outdir, "meta", f"meta.{tag}.tsv")
                        meta.to_csv(meta_path, sep="\t", index=False)

                        tar_chr = meta_keep["chr"].astype(str).to_numpy()
                        tar_side = meta_keep["best_side"].astype(str).to_numpy()
                        tar_summit = meta_keep["best_summit_1based"].to_numpy(dtype=np.int64)

                        tar_st1 = np.zeros_like(tar_summit)
                        tar_en1 = np.zeros_like(tar_summit)
                        for i in range(tar_summit.size):
                            if rel_a_use is not None:
                                s1, e1 = make_motif_window_1based_relative(int(tar_summit[i]), tar_side[i], rel_a_use, rel_b_use)
                            else:
                                s1, e1 = make_motif_window_1based(int(tar_summit[i]), tar_side[i], args.motif_outside_bp, args.motif_inside_bp)
                            tar_st1[i] = s1
                            tar_en1[i] = e1

                        ok, tar_st1c, tar_en1c = clip_1based_inclusive(tar_chr, tar_st1, tar_en1, chrlen)
                        meta_keep2 = meta_keep.iloc[np.where(ok)[0]].copy()
                        tar_chr2 = tar_chr[ok]
                        tar_st0 = (tar_st1c[ok] - 1).astype(np.int64)
                        tar_en1o = tar_en1c[ok].astype(np.int64)
                        tar_id = (meta_keep2["region_id"].astype(str) + "|" + meta_keep2["best_side"].astype(str)).to_list()

                        tar_df = pd.DataFrame({"chr": tar_chr2, "start0": tar_st0, "end1": tar_en1o, "id": tar_id})
                        tar_df = tar_df.sort_values(["chr", "start0"])
                        tar_bed = os.path.join(outdir, "bed_target", f"target.{tag}.bed")
                        tar_df.to_csv(tar_bed, sep="\t", index=False, header=False)

                        if args.homer_auto_bg:
                            eprint(f"[OK] {sm} {cl} dir={dir_cl} rel={rel_str} target={meta_keep2.shape[0]} tag={tag} (AUTO_BG=1)")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge",
                                "rel_range_used": rel_str,
                                "rel_range_used_input": "" if rel_a_raw is None else raw_text,
                                "note": "ok_auto_bg",
                                "n_target_arrays": int(df_cl.shape[0]),
                                "n_target_arrays_kept": int(meta_keep2.shape[0]),
                                "n_target_windows": int(tar_df.shape[0]),
                                "target_bed": tar_bed,
                            })
                            continue

                        df_bg_pool = df_sm[df_sm["cluster"] != cl].copy()
                        if df_bg_pool.shape[0] == 0:
                            eprint(f"[SKIP bg] {sm} {cl} rel={rel_str} no bg pool")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge",
                                "rel_range_used": rel_str,
                                "note": "skip_bg_pool_empty",
                                "n_target_arrays": int(df_cl.shape[0]),
                                "n_target_arrays_kept": int(meta_keep2.shape[0]),
                                "n_target_windows": int(tar_df.shape[0]),
                            })
                            continue

                        b_sL, b_sR = score_edges(
                            dir_cl,
                            df_bg_pool["m_outL"].to_numpy(),
                            df_bg_pool["m_inL"].to_numpy(),
                            df_bg_pool["m_inR"].to_numpy(),
                            df_bg_pool["m_outR"].to_numpy()
                        )
                        bg_meta = df_bg_pool[["region_id", "chr", "start", "end", "location", "summit_L_1based", "summit_R_1based"]].copy()
                        bg_meta["score_left"] = b_sL
                        bg_meta["score_right"] = b_sR
                        bg_meta["best_side"] = np.where(bg_meta["score_left"] >= bg_meta["score_right"], "L", "R")
                        bg_meta["best_score"] = np.maximum(bg_meta["score_left"], bg_meta["score_right"])
                        bg_meta["best_summit_1based"] = np.where(
                            bg_meta["best_side"] == "L",
                            bg_meta["summit_L_1based"],
                            bg_meta["summit_R_1based"]
                        ).astype(np.int64)

                        bg_meta = bg_meta[np.isfinite(bg_meta["best_score"])].copy()
                        if args.require_positive_score and dir_cl in ("inside_low", "inside_high"):
                            bg_meta = bg_meta[bg_meta["best_score"] > 0].copy()

                        if bg_meta.shape[0] == 0:
                            eprint(f"[SKIP bg] {sm} {cl} rel={rel_str} bg_meta=0 after filtering")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge",
                                "rel_range_used": rel_str,
                                "note": "skip_bg_empty_after_filter",
                                "n_target_arrays": int(df_cl.shape[0]),
                                "n_target_arrays_kept": int(meta_keep2.shape[0]),
                                "n_target_windows": int(tar_df.shape[0]),
                            })
                            continue

                        n_tar = meta_keep2.shape[0]
                        n_bg_need = int(math.ceil(float(args.bg_multiple) * float(n_tar)))
                        bg_pick = sample_background_regions(bg_meta, meta_keep2, n_bg_need, args.bg_match_location, seed=1)
                        if bg_pick.shape[0] == 0:
                            eprint(f"[SKIP bg] {sm} {cl} rel={rel_str} bg_pick=0")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge",
                                "rel_range_used": rel_str,
                                "note": "skip_bg_pick_empty",
                                "n_target_arrays": int(df_cl.shape[0]),
                                "n_target_arrays_kept": int(meta_keep2.shape[0]),
                                "n_target_windows": int(tar_df.shape[0]),
                            })
                            continue

                        bg_chr = bg_pick["chr"].astype(str).to_numpy()
                        bg_side = bg_pick["best_side"].astype(str).to_numpy()
                        bg_summit = bg_pick["best_summit_1based"].to_numpy(dtype=np.int64)

                        bg_st1 = np.zeros_like(bg_summit)
                        bg_en1 = np.zeros_like(bg_summit)
                        for i in range(bg_summit.size):
                            if rel_a_use is not None:
                                s1, e1 = make_motif_window_1based_relative(int(bg_summit[i]), bg_side[i], rel_a_use, rel_b_use)
                            else:
                                s1, e1 = make_motif_window_1based(int(bg_summit[i]), bg_side[i], args.motif_outside_bp, args.motif_inside_bp)
                            bg_st1[i] = s1
                            bg_en1[i] = e1

                        okb, bg_st1c, bg_en1c = clip_1based_inclusive(bg_chr, bg_st1, bg_en1, chrlen)
                        bg_pick2 = bg_pick.iloc[np.where(okb)[0]].copy()
                        bg_chr2 = bg_chr[okb]
                        bg_st0 = (bg_st1c[okb] - 1).astype(np.int64)
                        bg_en1o = bg_en1c[okb].astype(np.int64)
                        bg_id = (bg_pick2["region_id"].astype(str) + "|" + bg_pick2["best_side"].astype(str)).to_list()

                        bg_df = pd.DataFrame({"chr": bg_chr2, "start0": bg_st0, "end1": bg_en1o, "id": bg_id})
                        bg_df = bg_df.sort_values(["chr", "start0"])
                        bg_bed = os.path.join(outdir, "bed_bg", f"bg.{tag}.bed")
                        bg_df.to_csv(bg_bed, sep="\t", index=False, header=False)

                        eprint(f"[OK] {sm} {cl} dir={dir_cl} rel={rel_str} target={meta_keep2.shape[0]} bg={bg_pick2.shape[0]} tag={tag}")
                        summary_rows.append({
                            "sample": sm,
                            "cluster": cl,
                            "mode": "edge",
                            "rel_range_used": rel_str,
                            "rel_range_used_input": "" if rel_a_raw is None else raw_text,
                            "note": "ok",
                            "n_target_arrays": int(df_cl.shape[0]),
                            "n_target_arrays_kept": int(meta_keep2.shape[0]),
                            "n_target_windows": int(tar_df.shape[0]),
                            "n_bg_arrays": int(bg_pick2.shape[0]),
                            "n_bg_windows": int(bg_df.shape[0]),
                            "target_bed": tar_bed,
                            "bg_bed": bg_bed,
                        })

        if do_edge_outside or do_edge_inside:
            outside_bins = build_equal_distance_bins(args.outside_start_bp, args.outside_end_bp, args.outside_num_bins)
            inside_bins = build_equal_distance_bins_rel(args.inside_start_rel, args.inside_end_rel, args.inside_num_bins, err_prefix="inside")

            if args.strict_outside_only and args.outside_start_bp <= 0:
                raise RuntimeError("--strict_outside_only requires --outside_start_bp > 0")
            if args.inside_end_rel > 0:
                raise RuntimeError("--inside_end_rel should be <= 0 for inward windows")

            for sm in samples:
                df_sm = meta_by_sample.get(sm)
                if df_sm is None or df_sm.shape[0] == 0:
                    continue

                for cl in target_clusters:
                    df_tar = df_sm[df_sm["cluster"] == cl].copy()
                    n_tar_arrays = df_tar.shape[0]
                    dir_cl = dir_map.get(cl, args.cluster_dir_default)

                    base_tag = (
                        f"{sm}.{cl}.{tss_tag}.mode_edgeBased"
                        f".scoreF{args.score_flank_bp}.scoreE{args.score_edge_bp}"
                        f".q{args.q_keep}.bgx{args.bg_multiple}.loc{int(args.bg_match_location)}"
                        f".minTar{args.min_n_target_arrays}.pos{int(args.require_positive_score)}"
                        f".autoBG{int(args.homer_auto_bg)}"
                    )

                    if n_tar_arrays < args.min_n_target_arrays:
                        if do_edge_outside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_outside", "note": "skip_min_target_arrays", "n_target_arrays": int(n_tar_arrays)})
                        if do_edge_inside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_inside", "note": "skip_min_target_arrays", "n_target_arrays": int(n_tar_arrays)})
                        eprint(f"[SKIP] {sm} {cl}: n_target_arrays={n_tar_arrays} < {args.min_n_target_arrays}")
                        continue

                    sL, sR = score_edges(
                        dir_cl,
                        df_tar["m_outL"].to_numpy(),
                        df_tar["m_inL"].to_numpy(),
                        df_tar["m_inR"].to_numpy(),
                        df_tar["m_outR"].to_numpy(),
                    )

                    meta_tar = df_tar[[
                        "sample",
                        "cluster",
                        "region_id",
                        "chr",
                        "start",
                        "end",
                        "location",
                        "near_tss",
                        "summit_L_1based",
                        "summit_R_1based",
                    ]].copy()
                    meta_tar["score_left"] = sL
                    meta_tar["score_right"] = sR
                    meta_tar["best_side"] = np.where(meta_tar["score_left"] >= meta_tar["score_right"], "L", "R")
                    meta_tar["best_score"] = np.maximum(meta_tar["score_left"], meta_tar["score_right"])
                    meta_tar["best_summit_1based"] = np.where(
                        meta_tar["best_side"] == "L",
                        meta_tar["summit_L_1based"],
                        meta_tar["summit_R_1based"],
                    ).astype(np.int64)

                    meta_tar = meta_tar[np.isfinite(meta_tar["best_score"])].copy()
                    if args.require_positive_score and dir_cl in ("inside_low", "inside_high"):
                        meta_tar = meta_tar[meta_tar["best_score"] > 0].copy()
                    if meta_tar.shape[0] == 0:
                        eprint(f"[SKIP] {sm} {cl}: no target arrays left after score filtering")
                        if do_edge_outside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_outside", "note": "skip_no_arrays_after_score_filter", "n_target_arrays": int(n_tar_arrays)})
                        if do_edge_inside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_inside", "note": "skip_no_arrays_after_score_filter", "n_target_arrays": int(n_tar_arrays)})
                        continue

                    if float(args.q_keep) < 1:
                        cutoff = float(meta_tar["best_score"].quantile(1 - float(args.q_keep)))
                        meta_keep = meta_tar[meta_tar["best_score"] >= cutoff].copy()
                    else:
                        meta_keep = meta_tar.copy()

                    meta_keep = meta_keep.sort_values(["best_score", "region_id"], ascending=[False, True]).reset_index(drop=True)
                    meta_keep.to_csv(os.path.join(outdir, "meta", f"meta.{base_tag}.tsv"), sep="\t", index=False)

                    if meta_keep.shape[0] == 0:
                        eprint(f"[SKIP] {sm} {cl}: no target arrays left after q_keep filter")
                        if do_edge_outside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_outside", "note": "skip_no_arrays_after_qkeep", "n_target_arrays": int(n_tar_arrays)})
                        if do_edge_inside:
                            summary_rows.append({"sample": sm, "cluster": cl, "mode": "edge_inside", "note": "skip_no_arrays_after_qkeep", "n_target_arrays": int(n_tar_arrays)})
                        continue

                    bg_meta_pool = None
                    bg_pick_shared = None

                    if not args.homer_auto_bg:
                        df_bg_pool = df_sm[df_sm["cluster"] != cl].copy()
                        if df_bg_pool.shape[0] == 0:
                            bg_meta_pool = pd.DataFrame()
                        else:
                            b_sL, b_sR = score_edges(
                                dir_cl,
                                df_bg_pool["m_outL"].to_numpy(),
                                df_bg_pool["m_inL"].to_numpy(),
                                df_bg_pool["m_inR"].to_numpy(),
                                df_bg_pool["m_outR"].to_numpy(),
                            )
                            bg_meta_pool = df_bg_pool[[
                                "region_id",
                                "chr",
                                "start",
                                "end",
                                "location",
                                "summit_L_1based",
                                "summit_R_1based",
                            ]].copy()
                            bg_meta_pool["score_left"] = b_sL
                            bg_meta_pool["score_right"] = b_sR
                            bg_meta_pool["best_side"] = np.where(bg_meta_pool["score_left"] >= bg_meta_pool["score_right"], "L", "R")
                            bg_meta_pool["best_score"] = np.maximum(bg_meta_pool["score_left"], bg_meta_pool["score_right"])
                            bg_meta_pool["best_summit_1based"] = np.where(
                                bg_meta_pool["best_side"] == "L",
                                bg_meta_pool["summit_L_1based"],
                                bg_meta_pool["summit_R_1based"],
                            ).astype(np.int64)

                            bg_meta_pool = bg_meta_pool[np.isfinite(bg_meta_pool["best_score"])].copy()
                            if args.require_positive_score and dir_cl in ("inside_low", "inside_high"):
                                bg_meta_pool = bg_meta_pool[bg_meta_pool["best_score"] > 0].copy()

                            if args.bg_sampling_strategy == "shared" and bg_meta_pool.shape[0] > 0:
                                n_bg_need_shared = int(math.ceil(args.bg_multiple * float(meta_keep.shape[0])))
                                bg_pick_shared = sample_background_regions(
                                    bg_meta_pool,
                                    meta_keep,
                                    n_bg_need_shared,
                                    args.bg_match_location,
                                    seed=1,
                                )
                                if bg_pick_shared.shape[0] > 0:
                                    bg_pick_shared.to_csv(
                                        os.path.join(outdir, "meta", f"meta.{base_tag}.background.shared.tsv"),
                                        sep="\t",
                                        index=False,
                                    )

                    if do_edge_outside:
                        for bin_idx, rel_a, rel_b, bin_label in outside_bins:
                            tag = f"{base_tag}.mode_edgeOutside.{bin_label}"

                            tar_keep2, tar_bed_df = build_windows_from_best_side(meta_keep, rel_a, rel_b, chrlen, clip_to_array=False)
                            tar_bed = os.path.join(outdir, "bed_target", f"target.{tag}.bed")
                            write_bed(tar_bed_df, tar_bed)

                            if tar_bed_df.shape[0] == 0:
                                eprint(f"[SKIP] {sm} {cl} outside {bin_label}: no target windows after clipping")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_outside",
                                    "outside_bin_index": int(bin_idx),
                                    "outside_bin_label": bin_label,
                                    "outside_start_bp_bin": int(rel_a),
                                    "outside_end_bp_bin": int(rel_b),
                                    "note": "skip_no_target_windows",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(meta_keep.shape[0]),
                                    "n_target_windows": 0,
                                })
                                continue

                            if args.homer_auto_bg:
                                eprint(f"[OK] {sm} {cl} outside {bin_label}: target_arrays={n_tar_arrays}, kept={tar_keep2.shape[0]}, target_windows={tar_bed_df.shape[0]} (AUTO_BG=1)")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_outside",
                                    "outside_bin_index": int(bin_idx),
                                    "outside_bin_label": bin_label,
                                    "outside_start_bp_bin": int(rel_a),
                                    "outside_end_bp_bin": int(rel_b),
                                    "note": "ok_auto_bg",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                    "target_bed": tar_bed,
                                })
                                continue

                            if bg_meta_pool is None or bg_meta_pool.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} outside {bin_label}: bg empty after score filtering")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_outside",
                                    "outside_bin_index": int(bin_idx),
                                    "outside_bin_label": bin_label,
                                    "outside_start_bp_bin": int(rel_a),
                                    "outside_end_bp_bin": int(rel_b),
                                    "note": "skip_bg_empty_after_filter",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                })
                                continue

                            if args.bg_sampling_strategy == "shared":
                                bg_pick = bg_pick_shared.copy() if bg_pick_shared is not None else bg_meta_pool.iloc[0:0].copy()
                            else:
                                n_bg_need = int(math.ceil(args.bg_multiple * float(tar_keep2.shape[0])))
                                bg_pick = sample_background_regions(bg_meta_pool, tar_keep2, n_bg_need, args.bg_match_location, seed=1)

                            if bg_pick.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} outside {bin_label}: bg_pick empty")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_outside",
                                    "outside_bin_index": int(bin_idx),
                                    "outside_bin_label": bin_label,
                                    "outside_start_bp_bin": int(rel_a),
                                    "outside_end_bp_bin": int(rel_b),
                                    "note": "skip_bg_pick_empty",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                })
                                continue

                            bg_pick2, bg_bed_df = build_windows_from_best_side(bg_pick, rel_a, rel_b, chrlen, clip_to_array=False)
                            bg_bed = os.path.join(outdir, "bed_bg", f"bg.{tag}.bed")
                            write_bed(bg_bed_df, bg_bed)

                            if bg_bed_df.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} outside {bin_label}: bg windows empty after clipping")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_outside",
                                    "outside_bin_index": int(bin_idx),
                                    "outside_bin_label": bin_label,
                                    "outside_start_bp_bin": int(rel_a),
                                    "outside_end_bp_bin": int(rel_b),
                                    "note": "skip_bg_windows_empty",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                    "n_bg_arrays": int(bg_pick.shape[0]),
                                    "n_bg_windows": 0,
                                })
                                continue

                            eprint(f"[OK] {sm} {cl} outside {bin_label}: target_arrays={n_tar_arrays}, kept={tar_keep2.shape[0]}, bg_arrays={bg_pick2.shape[0]}, target_windows={tar_bed_df.shape[0]}, bg_windows={bg_bed_df.shape[0]}")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge_outside",
                                "outside_bin_index": int(bin_idx),
                                "outside_bin_label": bin_label,
                                "outside_start_bp_bin": int(rel_a),
                                "outside_end_bp_bin": int(rel_b),
                                "note": "ok",
                                "n_target_arrays": int(n_tar_arrays),
                                "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                "n_target_windows": int(tar_bed_df.shape[0]),
                                "n_bg_arrays": int(bg_pick2.shape[0]),
                                "n_bg_windows": int(bg_bed_df.shape[0]),
                                "target_bed": tar_bed,
                                "bg_bed": bg_bed,
                            })

                    if do_edge_inside:
                        for bin_idx, rel_a, rel_b, bin_label in inside_bins:
                            tag = f"{base_tag}.mode_edgeInside.{bin_label}"

                            tar_keep2, tar_bed_df = build_windows_from_best_side(
                                meta_keep,
                                rel_a,
                                rel_b,
                                chrlen,
                                clip_to_array=args.inside_clip_to_array,
                            )
                            tar_bed = os.path.join(outdir, "bed_target", f"target.{tag}.bed")
                            write_bed(tar_bed_df, tar_bed)

                            if tar_bed_df.shape[0] == 0:
                                eprint(f"[SKIP] {sm} {cl} inside {bin_label}: no target windows after clipping")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_inside",
                                    "inside_bin_index": int(bin_idx),
                                    "inside_bin_label": bin_label,
                                    "inside_start_rel_bin": int(rel_a),
                                    "inside_end_rel_bin": int(rel_b),
                                    "note": "skip_no_target_windows",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(meta_keep.shape[0]),
                                    "n_target_windows": 0,
                                })
                                continue

                            if args.homer_auto_bg:
                                eprint(f"[OK] {sm} {cl} inside {bin_label}: target_arrays={n_tar_arrays}, kept={tar_keep2.shape[0]}, target_windows={tar_bed_df.shape[0]} (AUTO_BG=1)")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_inside",
                                    "inside_bin_index": int(bin_idx),
                                    "inside_bin_label": bin_label,
                                    "inside_start_rel_bin": int(rel_a),
                                    "inside_end_rel_bin": int(rel_b),
                                    "note": "ok_auto_bg",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                    "target_bed": tar_bed,
                                })
                                continue

                            if bg_meta_pool is None or bg_meta_pool.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} inside {bin_label}: bg empty after score filtering")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_inside",
                                    "inside_bin_index": int(bin_idx),
                                    "inside_bin_label": bin_label,
                                    "inside_start_rel_bin": int(rel_a),
                                    "inside_end_rel_bin": int(rel_b),
                                    "note": "skip_bg_empty_after_filter",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                })
                                continue

                            if args.bg_sampling_strategy == "shared":
                                bg_pick = bg_pick_shared.copy() if bg_pick_shared is not None else bg_meta_pool.iloc[0:0].copy()
                            else:
                                n_bg_need = int(math.ceil(args.bg_multiple * float(tar_keep2.shape[0])))
                                bg_pick = sample_background_regions(bg_meta_pool, tar_keep2, n_bg_need, args.bg_match_location, seed=1)

                            if bg_pick.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} inside {bin_label}: bg_pick empty")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_inside",
                                    "inside_bin_index": int(bin_idx),
                                    "inside_bin_label": bin_label,
                                    "inside_start_rel_bin": int(rel_a),
                                    "inside_end_rel_bin": int(rel_b),
                                    "note": "skip_bg_pick_empty",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                })
                                continue

                            bg_pick2, bg_bed_df = build_windows_from_best_side(
                                bg_pick,
                                rel_a,
                                rel_b,
                                chrlen,
                                clip_to_array=args.inside_clip_to_array,
                            )
                            bg_bed = os.path.join(outdir, "bed_bg", f"bg.{tag}.bed")
                            write_bed(bg_bed_df, bg_bed)

                            if bg_bed_df.shape[0] == 0:
                                eprint(f"[WARN] {sm} {cl} inside {bin_label}: bg windows empty after clipping")
                                summary_rows.append({
                                    "sample": sm,
                                    "cluster": cl,
                                    "mode": "edge_inside",
                                    "inside_bin_index": int(bin_idx),
                                    "inside_bin_label": bin_label,
                                    "inside_start_rel_bin": int(rel_a),
                                    "inside_end_rel_bin": int(rel_b),
                                    "note": "skip_bg_windows_empty",
                                    "n_target_arrays": int(n_tar_arrays),
                                    "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                    "n_target_windows": int(tar_bed_df.shape[0]),
                                    "n_bg_arrays": int(bg_pick.shape[0]),
                                    "n_bg_windows": 0,
                                })
                                continue

                            eprint(f"[OK] {sm} {cl} inside {bin_label}: target_arrays={n_tar_arrays}, kept={tar_keep2.shape[0]}, bg_arrays={bg_pick2.shape[0]}, target_windows={tar_bed_df.shape[0]}, bg_windows={bg_bed_df.shape[0]}")
                            summary_rows.append({
                                "sample": sm,
                                "cluster": cl,
                                "mode": "edge_inside",
                                "inside_bin_index": int(bin_idx),
                                "inside_bin_label": bin_label,
                                "inside_start_rel_bin": int(rel_a),
                                "inside_end_rel_bin": int(rel_b),
                                "note": "ok",
                                "n_target_arrays": int(n_tar_arrays),
                                "n_target_arrays_kept": int(tar_keep2.shape[0]),
                                "n_target_windows": int(tar_bed_df.shape[0]),
                                "n_bg_arrays": int(bg_pick2.shape[0]),
                                "n_bg_windows": int(bg_bed_df.shape[0]),
                                "target_bed": tar_bed,
                                "bg_bed": bg_bed,
                            })

    summary_path = os.path.join(outdir, "meta", "build_summary.tsv")
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_path, sep="\t", index=False)

    meta_info = {
        "cluster_tsv": cluster_tsv,
        "cluster_run_dir": args.cluster_run_dir,
        "best_k_selection_tsv": best_k_path,
        "best_k": best_k,
        "cluster_column_used": cluster_col,
        "mode": args.mode,
        "target_clusters": ",".join(target_clusters),
        "annotation_file": args.annotation_file,
        "tss_near_bp": args.tss_near_bp,
        "tss_filter_enabled": int(tss_filter_enabled),
        "genome": args.genome,
        "chrom_sizes": args.chrom_sizes,
        "positions_map_used": positions_map_used_raw,
        "internal_region_mode": args.internal_region_mode,
        "internal_nuc_pad_bp": args.internal_nuc_pad_bp,
        "outside_start_bp": args.outside_start_bp,
        "outside_end_bp": args.outside_end_bp,
        "outside_num_bins": args.outside_num_bins,
        "inside_start_rel": args.inside_start_rel,
        "inside_end_rel": args.inside_end_rel,
        "inside_num_bins": args.inside_num_bins,
        "inside_clip_to_array": int(args.inside_clip_to_array),
        "homer_auto_bg": int(args.homer_auto_bg),
        "bg_sampling_strategy": args.bg_sampling_strategy,
        "run_known_only": int(args.run_known_only),
        "rel_range": args.rel_range,
        "range_offset_bp": args.range_offset_bp,
        "motif_outside_bp": args.motif_outside_bp,
        "motif_inside_bp": args.motif_inside_bp,
    }
    pd.DataFrame([meta_info]).to_csv(os.path.join(outdir, "meta", "run_metadata.tsv"), sep="\t", index=False)

    homer_stat = {"n_homer_jobs_total": 0, "n_homer_jobs_completed": 0}
    if not args.build_only:
        homer_stat = run_homer_from_summary(args, summary_df, outdir)

    eprint("\n=== DONE ===")
    eprint(f"Summary: {summary_path}")
    eprint(f"Target BED dir: {os.path.join(outdir, 'bed_target')}")
    eprint(f"Background BED dir: {os.path.join(outdir, 'bed_bg')}")
    eprint(f"HOMER output dir: {os.path.join(outdir, 'homer_out')}")
    eprint(f"HOMER logs dir: {os.path.join(outdir, 'homer_logs')}")
    eprint(f"HOMER jobs: {homer_stat['n_homer_jobs_completed']}/{homer_stat['n_homer_jobs_total']} completed")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        eprint("FATAL:", str(e))
        sys.exit(1)
