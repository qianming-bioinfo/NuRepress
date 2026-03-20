#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(mclust)
  library(cluster)
  library(ggplot2)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

print_help <- function() {
  cat(
"Usage:
  Rscript cluster_array_subtype_autoK.R --sample_sheet SAMPLE_SHEET --out_dir OUT_DIR [options]

Required:
  --sample_sheet FILE        Tabular sample sheet with columns:
                             sample, regions_bed, atac_bw
  --out_dir DIR             Output directory

Optional:
  --sample_order STR        Comma-separated sample order. Default: sample_sheet order
  --ks STR                  Comma-separated k values to evaluate. Default: 2,3,4,5,6
  --best_k_source STR       Which silhouette curve determines best k:
                             within_sample | rank | global
                             Default: within_sample
  --seed INT                Random seed. Default: 1
  --eps FLOAT               Epsilon for log-ratio. Default: 1e-9
  --outside_flank_bp INT    Outside flank length. Default: 1000
  --sr_flank_bp INT         Metaplot flank length. Default: 1000
  --sr_flank_step INT       Metaplot flank step. Default: 10
  --sr_body_bins INT        Metaplot body bins. Default: 100
  --sr_min_regions_per_cluster_meta INT
                            Minimum regions per cluster for metaplot. Default: 50
  --sr_max_regions_per_cluster_meta INT
                            Maximum sampled regions per cluster for metaplot. Default: 10000
  --sil_subsample_n INT     Subsample size for silhouette. Default: 8000
  --mclust_max_fit_n INT    Max rows for direct mclust fitting. Default: 60000
  --max_points_per_sample INT
                            Max points per sample for genome scatter. Default: 200000
  --max_points_per_sample_cs INT
                            Max points per sample for cluster-space scatter. Default: 200000
  --sample_sheet_delim STR  auto | tab | comma. Default: auto
  -h, --help                Show help

Example:
  Rscript cluster_array_subtype_autoK.R \
    --sample_sheet /path/to/sample_sheet.tsv \
    --out_dir /path/to/output_dir \
    --sample_order sample1,sample2,sample3 \
    --ks 2,3,4,5,6

Sample sheet example:
  sample\tregions_bed\tatac_bw
  sample1\t/path/to/sample1_regions.bed\t/path/to/sample1_atac.bw
  sample2\t/path/to/sample2_regions.bed\t/path/to/sample2_atac.bw
",
  sep = "")
}

parse_args <- function(args) {
  opt <- list(
    sample_sheet = NULL,
    out_dir = NULL,
    sample_order = NULL,
    ks = "2,3,4,5,6",
    best_k_source = "within_sample",
    seed = 1L,
    eps = 1e-9,
    outside_flank_bp = 1000L,
    sr_flank_bp = 1000L,
    sr_flank_step = 10L,
    sr_body_bins = 100L,
    sr_min_regions_per_cluster_meta = 50L,
    sr_max_regions_per_cluster_meta = 10000L,
    sil_subsample_n = 8000L,
    mclust_max_fit_n = 60000L,
    max_points_per_sample = 200000L,
    max_points_per_sample_cs = 200000L,
    sample_sheet_delim = "auto"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      print_help()
      quit(save = "no", status = 0)
    }
    if (!startsWith(key, "-")) stop("Unknown positional argument: ", key)
    if (i == length(args)) stop("Missing value for option: ", key)
    val <- args[[i + 1L]]
    switch(key,
      "--sample_sheet" = { opt$sample_sheet <- val; i <- i + 2L },
      "--out_dir" = { opt$out_dir <- val; i <- i + 2L },
      "--sample_order" = { opt$sample_order <- val; i <- i + 2L },
      "--ks" = { opt$ks <- val; i <- i + 2L },
      "--best_k_source" = { opt$best_k_source <- val; i <- i + 2L },
      "--seed" = { opt$seed <- as.integer(val); i <- i + 2L },
      "--eps" = { opt$eps <- as.numeric(val); i <- i + 2L },
      "--outside_flank_bp" = { opt$outside_flank_bp <- as.integer(val); i <- i + 2L },
      "--sr_flank_bp" = { opt$sr_flank_bp <- as.integer(val); i <- i + 2L },
      "--sr_flank_step" = { opt$sr_flank_step <- as.integer(val); i <- i + 2L },
      "--sr_body_bins" = { opt$sr_body_bins <- as.integer(val); i <- i + 2L },
      "--sr_min_regions_per_cluster_meta" = { opt$sr_min_regions_per_cluster_meta <- as.integer(val); i <- i + 2L },
      "--sr_max_regions_per_cluster_meta" = { opt$sr_max_regions_per_cluster_meta <- as.integer(val); i <- i + 2L },
      "--sil_subsample_n" = { opt$sil_subsample_n <- as.integer(val); i <- i + 2L },
      "--mclust_max_fit_n" = { opt$mclust_max_fit_n <- as.integer(val); i <- i + 2L },
      "--max_points_per_sample" = { opt$max_points_per_sample <- as.integer(val); i <- i + 2L },
      "--max_points_per_sample_cs" = { opt$max_points_per_sample_cs <- as.integer(val); i <- i + 2L },
      "--sample_sheet_delim" = { opt$sample_sheet_delim <- val; i <- i + 2L },
      stop("Unknown option: ", key)
    )
  }

  if (is.null(opt$sample_sheet) || is.null(opt$out_dir)) {
    print_help()
    stop("--sample_sheet and --out_dir are required")
  }

  opt$best_k_source <- match.arg(opt$best_k_source, c("within_sample", "rank", "global"))
  opt$sample_sheet_delim <- match.arg(opt$sample_sheet_delim, c("auto", "tab", "comma"))

  ks <- as.integer(strsplit(opt$ks, ",", fixed = TRUE)[[1]])
  ks <- ks[is.finite(ks)]
  ks <- sort(unique(ks))
  if (length(ks) == 0) stop("No valid k values found in --ks")
  if (any(ks < 2L)) stop("All k values in --ks must be >= 2")
  opt$ks <- ks

  if (!is.null(opt$sample_order)) {
    so <- strsplit(opt$sample_order, ",", fixed = TRUE)[[1]]
    so <- trimws(so)
    so <- so[nzchar(so)]
    if (length(so) == 0) stop("--sample_order is empty after parsing")
    opt$sample_order <- so
  }

  opt
}

trim_ws_df <- function(df) {
  for (nm in names(df)) {
    if (is.character(df[[nm]])) df[[nm]] <- trimws(df[[nm]])
  }
  df
}

read_sample_sheet <- function(path, delim_mode = "auto") {
  if (!file.exists(path)) stop("File not found: ", path)
  delim <- delim_mode
  if (delim == "auto") {
    if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) delim <- "tab" else delim <- "comma"
  }
  df <- switch(delim,
    tab = readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    comma = readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- trim_ws_df(df)
  req <- c("sample", "regions_bed", "atac_bw")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("sample_sheet missing columns: ", paste(miss, collapse = ", "))
  df
}

stop_if_missing <- function(p) if (!file.exists(p)) stop("File not found: ", p)

bw_import_as_rlelist <- function(path) {
  x <- rtracklayer::import(path, as = "RleList")
  if (inherits(x, "RleList")) return(x)
  if (inherits(x, "SimpleRleList")) return(as(x, "RleList"))
  x2 <- rtracklayer::import(path)
  if (inherits(x2, "RleList")) return(x2)
  if (inherits(x2, "SimpleRleList")) return(as(x2, "RleList"))
  if (inherits(x2, "GRanges")) {
    sc <- S4Vectors::mcols(x2)$score
    if (is.null(sc)) sc <- rep(1, length(x2))
    seqlens <- GenomeInfoDb::seqlengths(x2)
    if (all(is.na(seqlens))) {
      seqlens <- tapply(GenomicRanges::end(x2), as.character(GenomicRanges::seqnames(x2)), max)
    }
    seqlens <- seqlens[!is.na(seqlens)]
    cov <- GenomicRanges::coverage(x2, weight = sc, width = seqlens)
    return(as(cov, "RleList"))
  }
  stop("bw_import_as_rlelist: unsupported import type: ", paste(class(x2), collapse = ", "))
}

read_bed_any_to_df <- function(bed_path, sample) {
  df <- readr::read_tsv(
    bed_path,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  )
  if (ncol(df) < 3) stop("BED must have >=3 columns: ", bed_path)
  chr   <- as.character(df[[1]])
  start <- as.integer(df[[2]])
  end   <- as.integer(df[[3]])
  id <- if (ncol(df) >= 4) as.character(df[[4]]) else NA_character_
  if (all(is.na(id)) || all(id == "")) {
    id <- paste0(sample, "_region_", seq_along(chr))
  } else {
    id[id == "" | is.na(id)] <- paste0(sample, "_region_", which(id == "" | is.na(id)))
  }
  out <- data.frame(chr = chr, start = start, end = end, region_id = id, stringsAsFactors = FALSE)
  out <- out[is.finite(out$start) & is.finite(out$end) & out$start <= out$end, , drop = FALSE]
  out
}

meanSignalFast <- function(gr, cov_rlelist) {
  gr <- as(gr, "GRanges")
  res <- rep(NA_real_, length(gr))
  chrs <- as.character(seqnames(gr))
  idx_by_chr <- split(seq_along(gr), chrs)

  for (chr in names(idx_by_chr)) {
    ii <- idx_by_chr[[chr]]
    if (!chr %in% names(cov_rlelist)) next
    covchr <- cov_rlelist[[chr]]
    chr_len <- length(covchr)

    r <- ranges(gr)[ii]
    s <- start(r); e <- end(r)
    s[s < 1] <- 1L
    e[e > chr_len] <- chr_len

    ok <- s <= e
    tmp <- rep(NA_real_, length(ii))
    if (any(ok)) {
      vv <- IRanges::Views(covchr, start = s[ok], end = e[ok])
      tmp[ok] <- as.numeric(IRanges::viewMeans(vv, na.rm = TRUE))
    }
    res[ii] <- tmp
  }
  res
}

fill_na_inf_median <- function(df_num) {
  df_num <- as.data.frame(df_num)
  for (cc in colnames(df_num)) {
    x <- df_num[[cc]]
    x[is.infinite(x)] <- NA_real_
    if (all(is.na(x))) {
      x <- rep(0, length(x))
    } else {
      med <- median(x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      x[is.na(x)] <- med
    }
    df_num[[cc]] <- x
  }
  df_num
}

safe_scale <- function(mat) {
  mat <- as.matrix(mat)
  out <- mat
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    mu <- mean(x)
    sdv <- stats::sd(x)
    if (!is.finite(sdv) || sdv == 0) out[, j] <- 0 else out[, j] <- (x - mu) / sdv
  }
  out
}

percent_rank_vec <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  n_ok <- sum(ok)
  out <- rep(NA_real_, length(x))
  if (n_ok <= 1) {
    out[ok] <- 0
    return(out)
  }
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- (r - 1) / (n_ok - 1)
  out
}

get_fixed_sub_idx <- function(n, max_n, seed, save_path = NULL) {
  if (n <= max_n) {
    idx <- seq_len(n)
  } else {
    set.seed(seed)
    idx <- sort(sample.int(n, max_n))
  }
  if (!is.null(save_path)) writeLines(as.character(idx), con = save_path)
  idx
}

sil_curve_mclust_subsample <- function(mat, ks, seed = 1, max_n = 8000L, idx_save = NULL) {
  n <- nrow(mat)
  idx_sub <- get_fixed_sub_idx(n, max_n = max_n, seed = seed, save_path = idx_save)
  mat_s <- mat[idx_sub, , drop = FALSE]
  d <- stats::dist(mat_s)
  sapply(ks, function(k) {
    mc <- mclust::Mclust(mat_s, G = k)
    ss <- cluster::silhouette(mc$classification, d)
    mean(ss[, "sil_width"])
  })
}

save_sil_plot <- function(ks, sil, main, ylab, out_pdf, out_png = NULL) {
  grDevices::pdf(out_pdf, width = 6.8, height = 4.8)
  plot(ks, sil, type = "b", xlab = "k", ylab = ylab, main = main)
  grDevices::dev.off()
  if (!is.null(out_png)) {
    grDevices::png(out_png, width = 1200, height = 850, res = 150)
    plot(ks, sil, type = "b", xlab = "k", ylab = ylab, main = main)
    grDevices::dev.off()
  }
}

fit_mclust_predict <- function(mat, G, seed = 1, max_fit_n = 60000L) {
  n <- nrow(mat)
  if (n <= max_fit_n) {
    set.seed(seed)
    return(mclust::Mclust(mat, G = G))
  }
  set.seed(seed)
  idx <- sample.int(n, max_fit_n)
  mc_sub <- mclust::Mclust(mat[idx, , drop = FALSE], G = G)
  pr <- predict(mc_sub, newdata = mat)
  mc_sub$classification <- pr$classification
  mc_sub$z <- pr$z
  mc_sub
}

sr_make_labels_col <- function(flank_bp, flank_step, body_bins) {
  n_flank_bins <- as.integer(ceiling(flank_bp / flank_step))
  n_total <- n_flank_bins + body_bins + n_flank_bins
  lab <- rep("", n_total)

  lab[1] <- paste0("-", n_flank_bins * flank_step)
  mid_up <- max(1, floor(n_flank_bins/2))
  lab[mid_up] <- paste0("-", mid_up * flank_step)
  lab[n_flank_bins] <- paste0("-", flank_step)

  b0 <- n_flank_bins + 1
  b1 <- n_flank_bins + body_bins
  lab[b0] <- "Body:0%"
  lab[b0 + floor((body_bins-1)/2)] <- "50%"
  lab[b1] <- "100%"

  d0 <- n_flank_bins + body_bins + 1
  d1 <- n_total
  lab[d0] <- paste0("+", flank_step)
  mid_dn <- d0 + floor(n_flank_bins/2)
  if (mid_dn <= d1) lab[mid_dn] <- paste0("+", (mid_dn - d0 + 1) * flank_step)
  lab[d1] <- paste0("+", n_flank_bins * flank_step)

  list(
    labels_col = lab,
    n_flank_bins = n_flank_bins,
    body_start = b0,
    body_end = b1,
    down_start = d0,
    total_cols = n_total
  )
}

sr_signal_matrix_scale_regions <- function(gr, cov_rlelist,
                                           flank_bp = 1000,
                                           flank_step = 10,
                                           body_bins = 100) {
  if (!inherits(cov_rlelist, "RleList")) {
    if (inherits(cov_rlelist, "SimpleRleList")) cov_rlelist <- as(cov_rlelist, "RleList")
    else stop("cov_rlelist must be RleList")
  }
  gr <- as(gr, "GRanges")
  n <- length(gr)
  if (n == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  n_flank_bins <- as.integer(ceiling(flank_bp / flank_step))
  up_offsets <- seq(-n_flank_bins * flank_step, -flank_step, by = flank_step)
  down_offsets <- seq(flank_step, n_flank_bins * flank_step, by = flank_step)

  coln <- c(as.character(up_offsets), paste0("B", seq_len(body_bins)), as.character(down_offsets))
  mat <- matrix(NA_real_, nrow = n, ncol = length(coln))
  colnames(mat) <- coln

  chrs <- as.character(seqnames(gr))
  idx_by_chr <- split(seq_len(n), chrs)

  for (chr in names(idx_by_chr)) {
    ii <- idx_by_chr[[chr]]
    if (!chr %in% names(cov_rlelist)) next
    covchr <- cov_rlelist[[chr]]
    chr_len <- length(covchr)

    for (k in ii) {
      s0 <- start(gr)[k]
      e0 <- end(gr)[k]
      if (!is.finite(s0) || !is.finite(e0) || s0 > e0) next

      up_pos <- s0 + up_offsets
      L <- e0 - s0 + 1L
      centers <- s0 + floor(((seq_len(body_bins) - 0.5) * L) / body_bins)
      centers[centers < s0] <- s0
      centers[centers > e0] <- e0
      down_pos <- e0 + down_offsets

      pos <- c(up_pos, centers, down_pos)
      ok <- (pos >= 1) & (pos <= chr_len)
      vals <- rep(NA_real_, length(pos))
      if (any(ok)) vals[ok] <- as.numeric(covchr[pos[ok]])
      mat[k, ] <- vals
    }
  }
  mat
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

  samples <- read_sample_sheet(opt$sample_sheet, opt$sample_sheet_delim)
  if (is.null(opt$sample_order)) {
    sample_order <- samples$sample
  } else {
    sample_order <- opt$sample_order
    missing_order <- setdiff(sample_order, samples$sample)
    extra_sheet <- setdiff(samples$sample, sample_order)
    if (length(missing_order) > 0) stop("sample_order contains values absent from sample_sheet: ", paste(missing_order, collapse = ", "))
    if (length(extra_sheet) > 0) stop("sample_sheet contains values absent from sample_order: ", paste(extra_sheet, collapse = ", "))
    samples <- samples[match(sample_order, samples$sample), , drop = FALSE]
  }

  for (i in seq_len(nrow(samples))) {
    stop_if_missing(samples$regions_bed[i])
    stop_if_missing(samples$atac_bw[i])
  }
  stopifnot(all(samples$sample == sample_order))

  metadata_file <- file.path(opt$out_dir, "run_metadata.txt")
  writeLines(c(
    paste0("sample_sheet=", normalizePath(opt$sample_sheet, winslash = "/", mustWork = FALSE)),
    paste0("out_dir=", normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE)),
    paste0("sample_order=", paste(sample_order, collapse = ",")),
    paste0("ks=", paste(opt$ks, collapse = ",")),
    paste0("best_k_source=", opt$best_k_source),
    paste0("seed=", opt$seed),
    paste0("eps=", opt$eps),
    paste0("outside_flank_bp=", opt$outside_flank_bp),
    paste0("sr_flank_bp=", opt$sr_flank_bp),
    paste0("sr_flank_step=", opt$sr_flank_step),
    paste0("sr_body_bins=", opt$sr_body_bins),
    paste0("sr_min_regions_per_cluster_meta=", opt$sr_min_regions_per_cluster_meta),
    paste0("sr_max_regions_per_cluster_meta=", opt$sr_max_regions_per_cluster_meta),
    paste0("sil_subsample_n=", opt$sil_subsample_n),
    paste0("mclust_max_fit_n=", opt$mclust_max_fit_n)
  ), con = metadata_file)

  message(">>> [1] Read BEDs and compute ATAC_inside/left/right from bigWig")

  atac_cov_map <- vector("list", length(sample_order))
  names(atac_cov_map) <- sample_order
  sample_long_list <- vector("list", length(sample_order))

  for (i in seq_along(sample_order)) {
    ct <- sample_order[i]
    message("  - sample: ", ct)

    bed_df <- read_bed_any_to_df(samples$regions_bed[i], ct)
    if (nrow(bed_df) == 0) stop("Empty BED for sample: ", ct)

    message("    importing ATAC bw as RleList: ", samples$atac_bw[i])
    atac_cov <- bw_import_as_rlelist(samples$atac_bw[i])
    atac_cov_map[[ct]] <- atac_cov

    gr_inside <- GRanges(
      seqnames = bed_df$chr,
      ranges = IRanges(start = as.integer(bed_df$start) + 1L, end = as.integer(bed_df$end)),
      region_id = bed_df$region_id
    )

    s1 <- start(gr_inside)
    e1 <- end(gr_inside)
    left_start <- s1 - opt$outside_flank_bp
    left_end   <- s1 - 1L
    right_start <- e1 + 1L
    right_end   <- e1 + opt$outside_flank_bp

    gr_left <- GRanges(seqnames = seqnames(gr_inside), ranges = IRanges(start = left_start, end = left_end), region_id = bed_df$region_id)
    gr_right <- GRanges(seqnames = seqnames(gr_inside), ranges = IRanges(start = right_start, end = right_end), region_id = bed_df$region_id)

    inside_atac <- meanSignalFast(gr_inside, atac_cov)
    left_atac   <- meanSignalFast(gr_left, atac_cov)
    right_atac  <- meanSignalFast(gr_right, atac_cov)

    df <- data.frame(
      sample    = ct,
      region_id = bed_df$region_id,
      chr       = bed_df$chr,
      start     = bed_df$start,
      end       = bed_df$end,
      width_bp  = as.integer(bed_df$end) - as.integer(bed_df$start),
      ATAC_inside  = inside_atac,
      ATAC_left1k  = left_atac,
      ATAC_right1k = right_atac,
      stringsAsFactors = FALSE
    )

    df$outside_left  <- df$ATAC_left1k
    df$outside_right <- df$ATAC_right1k
    df$outside_mean  <- rowMeans(cbind(df$outside_left, df$outside_right), na.rm = TRUE)
    df$z1 <- log2((df$ATAC_inside + opt$eps) / (df$outside_mean + opt$eps))
    df$z2 <- abs(log2((df$outside_right + opt$eps) / (df$outside_left + opt$eps)))
    df$edge_logratio <- log2((df$outside_mean + opt$eps) / (df$ATAC_inside + opt$eps))

    sample_long_list[[i]] <- df
  }

  sample_long <- dplyr::bind_rows(sample_long_list)
  sample_long$sample <- factor(sample_long$sample, levels = sample_order)

  na_stat <- sample_long %>%
    summarise(
      n = n(),
      na_inside = sum(!is.finite(ATAC_inside)),
      na_left = sum(!is.finite(outside_left)),
      na_right = sum(!is.finite(outside_right)),
      na_outmean = sum(!is.finite(outside_mean))
    )
  write_tsv(na_stat, file.path(opt$out_dir, "NA_diagnostic.ATAC_inside_left_right_outmean.tsv"))

  message(">>> [2] Build feature matrices")

  z1_rank <- rep(NA_real_, nrow(sample_long))
  z2_rank <- rep(NA_real_, nrow(sample_long))
  for (ct in sample_order) {
    idx <- which(as.character(sample_long$sample) == ct)
    z1_rank[idx] <- percent_rank_vec(sample_long$z1[idx])
    z2_rank[idx] <- percent_rank_vec(sample_long$z2[idx])
  }
  feat_rank_df <- data.frame(z1_rank = z1_rank, z2_rank = z2_rank, stringsAsFactors = FALSE)
  feat_rank_df <- fill_na_inf_median(feat_rank_df)
  feat_rank <- safe_scale(feat_rank_df)

  feat_raw_df <- sample_long[, c("z1", "z2"), drop = FALSE]
  feat_raw_df <- fill_na_inf_median(feat_raw_df)
  feat_norank_global <- safe_scale(feat_raw_df)

  z1_ws <- rep(0, nrow(sample_long))
  z2_ws <- rep(0, nrow(sample_long))
  for (ct in sample_order) {
    idx <- which(as.character(sample_long$sample) == ct)
    x1 <- sample_long$z1[idx]; x2 <- sample_long$z2[idx]
    x1[!is.finite(x1)] <- NA_real_
    x2[!is.finite(x2)] <- NA_real_
    med1 <- median(x1, na.rm = TRUE); if (!is.finite(med1)) med1 <- 0
    med2 <- median(x2, na.rm = TRUE); if (!is.finite(med2)) med2 <- 0
    z1_ws[idx] <- as.numeric(scale(ifelse(is.na(x1), med1, x1)))
    z2_ws[idx] <- as.numeric(scale(ifelse(is.na(x2), med2, x2)))
  }
  z1_ws[!is.finite(z1_ws)] <- 0
  z2_ws[!is.finite(z2_ws)] <- 0
  feat_ws <- cbind(z1_ws = z1_ws, z2_ws = z2_ws)

  message(">>> [3] Silhouette vs k")
  set.seed(opt$seed)

  sil_rank <- sil_curve_mclust_subsample(
    feat_rank, opt$ks, seed = opt$seed,
    max_n = opt$sil_subsample_n,
    idx_save = file.path(opt$out_dir, "silhouette_subsample_idx.rank_z1z2.txt")
  )
  df_sil_rank <- data.frame(k = opt$ks, mean_silhouette_rank = as.numeric(sil_rank))
  save_sil_plot(
    opt$ks, sil_rank,
    main = "Mean silhouette vs k (mclust, rank z1/z2) [subsampled]",
    ylab = "mean silhouette (rank z1/z2)",
    out_pdf = file.path(opt$out_dir, "silhouette.rank_z1z2.mclust.pdf"),
    out_png = file.path(opt$out_dir, "silhouette.rank_z1z2.mclust.png")
  )
  write_tsv(df_sil_rank, file.path(opt$out_dir, "silhouette.rank_z1z2.mclust.tsv"))

  sil_norank_global <- sil_curve_mclust_subsample(
    feat_norank_global, opt$ks, seed = opt$seed,
    max_n = opt$sil_subsample_n,
    idx_save = file.path(opt$out_dir, "silhouette_subsample_idx.norank_global_z1z2.txt")
  )
  df_sil_norank_global <- data.frame(k = opt$ks, mean_silhouette_norank_global = as.numeric(sil_norank_global))
  save_sil_plot(
    opt$ks, sil_norank_global,
    main = "Mean silhouette vs k (mclust, global z-score z1/z2) [subsampled]",
    ylab = "mean silhouette (global z-score z1/z2)",
    out_pdf = file.path(opt$out_dir, "silhouette.norank_global_z1z2.mclust.pdf"),
    out_png = file.path(opt$out_dir, "silhouette.norank_global_z1z2.mclust.png")
  )
  write_tsv(df_sil_norank_global, file.path(opt$out_dir, "silhouette.norank_global_z1z2.mclust.tsv"))

  sil_ws <- sil_curve_mclust_subsample(
    feat_ws, opt$ks, seed = opt$seed,
    max_n = opt$sil_subsample_n,
    idx_save = file.path(opt$out_dir, "silhouette_subsample_idx.within_sample_z1z2.txt")
  )
  df_sil_ws <- data.frame(k = opt$ks, mean_silhouette_within_sample = as.numeric(sil_ws))
  save_sil_plot(
    opt$ks, sil_ws,
    main = "Mean silhouette vs k (mclust, within-sample scaled z1/z2) [subsampled]",
    ylab = "mean silhouette (within-sample scaled z1/z2)",
    out_pdf = file.path(opt$out_dir, "silhouette.norank_within_sample_z1z2.mclust.pdf"),
    out_png = file.path(opt$out_dir, "silhouette.norank_within_sample_z1z2.mclust.png")
  )
  write_tsv(df_sil_ws, file.path(opt$out_dir, "silhouette.norank_within_sample_z1z2.mclust.tsv"))

  df_summary <- df_sil_rank %>%
    left_join(df_sil_norank_global, by = "k") %>%
    left_join(df_sil_ws, by = "k")
  write_tsv(df_summary, file.path(opt$out_dir, "silhouette_evidence_summary.tsv"))

  best_k <- switch(
    opt$best_k_source,
    rank = df_sil_rank$k[which.max(df_sil_rank$mean_silhouette_rank)],
    global = df_sil_norank_global$k[which.max(df_sil_norank_global$mean_silhouette_norank_global)],
    within_sample = df_sil_ws$k[which.max(df_sil_ws$mean_silhouette_within_sample)]
  )
  best_k <- as.integer(best_k)

  fit_mat <- switch(
    opt$best_k_source,
    rank = feat_rank,
    global = feat_norank_global,
    within_sample = feat_ws
  )
  fit_space_label <- switch(
    opt$best_k_source,
    rank = "rank_scaled_z1z2",
    global = "global_scaled_z1z2",
    within_sample = "within_sample_scaled_z1z2"
  )
  fit_space_xlab <- switch(
    opt$best_k_source,
    rank = "rank-scaled z1",
    global = "global-scaled z1",
    within_sample = "within-sample scaled z1"
  )
  fit_space_ylab <- switch(
    opt$best_k_source,
    rank = "rank-scaled z2",
    global = "global-scaled z2",
    within_sample = "within-sample scaled z2"
  )
  fit_space_coord1 <- fit_mat[, 1]
  fit_space_coord2 <- fit_mat[, 2]

  best_k_tbl <- data.frame(
    best_k_source = opt$best_k_source,
    best_k = best_k,
    best_silhouette = switch(
      opt$best_k_source,
      rank = max(df_sil_rank$mean_silhouette_rank, na.rm = TRUE),
      global = max(df_sil_norank_global$mean_silhouette_norank_global, na.rm = TRUE),
      within_sample = max(df_sil_ws$mean_silhouette_within_sample, na.rm = TRUE)
    ),
    fit_space = fit_space_label
  )
  write_tsv(best_k_tbl, file.path(opt$out_dir, "best_k_selection.tsv"))
  write(paste0("best_k=", best_k), file = metadata_file, append = TRUE)
  write(paste0("best_k_source=", opt$best_k_source), file = metadata_file, append = TRUE)

  message(">>> [4] Fit best k = ", best_k, " on ", fit_space_label)
  mc_best <- fit_mclust_predict(fit_mat, G = best_k, seed = opt$seed, max_fit_n = opt$mclust_max_fit_n)
  cl_best <- paste0("C", mc_best$classification)

  tab <- table(as.character(sample_long$sample), cl_best)
  prop <- prop.table(tab, 1)
  txt_out <- file.path(opt$out_dir, paste0("k", best_k, "_sample_association.txt"))
  cat(paste0("=== table(sample, k=", best_k, " cluster) ===\n"), file = txt_out)
  capture.output(tab, file = txt_out, append = TRUE)
  cat("\n=== prop.table(tab, 1) ===\n", file = txt_out, append = TRUE)
  capture.output(prop, file = txt_out, append = TRUE)

  message(">>> [5] Plot z1/z2 scatter colored by best k")
  df_scatter <- sample_long
  df_scatter$cluster <- cl_best
  df_scatter$sample <- factor(as.character(df_scatter$sample), levels = sample_order)

  p_sc <- ggplot(df_scatter, aes(x = z1, y = z2, color = cluster)) +
    geom_point(size = 0.4, alpha = 0.25) +
    facet_wrap(~ sample, nrow = 1, drop = FALSE) +
    theme_bw() +
    labs(
      title = paste0("z1/z2 log-ratio embedding colored by k=", best_k, " cluster"),
      x = "z1 = log2((inside+eps)/(outside_mean+eps))",
      y = "z2 = abs(log2((right+eps)/(left+eps)))",
      color = paste0("k=", best_k, " cluster")
    )
  pdf(file.path(opt$out_dir, paste0("z1z2.scatter.bySample.k", best_k, ".pdf")), width = 12.5, height = 4.8)
  print(p_sc)
  dev.off()

  message(">>> [6] ATAC metaplot by sample x cluster")
  meta_dir <- file.path(opt$out_dir, "ATAC_metaplot_by_sample_cluster")
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

  labinfo <- sr_make_labels_col(opt$sr_flank_bp, opt$sr_flank_step, opt$sr_body_bins)
  tick_pos <- which(labinfo$labels_col != "")
  tick_lab <- labinfo$labels_col[tick_pos]

  meta_list_all <- list()
  for (ct in sample_order) {
    message("  - metaplot sample: ", ct)
    df_s <- df_scatter[df_scatter$sample == ct & !is.na(df_scatter$cluster), , drop = FALSE]
    if (nrow(df_s) == 0) next

    atac_cov <- atac_cov_map[[ct]]
    if (is.null(atac_cov) || !inherits(atac_cov, "RleList")) stop("ATAC cov missing or invalid for sample: ", ct)

    gr_s <- GRanges(seqnames = df_s$chr, ranges = IRanges(start = as.integer(df_s$start) + 1L, end = as.integer(df_s$end)))

    for (cl in sort(unique(df_s$cluster))) {
      idx <- which(df_s$cluster == cl)
      if (length(idx) < opt$sr_min_regions_per_cluster_meta) next
      sub_idx <- idx
      if (length(sub_idx) > opt$sr_max_regions_per_cluster_meta) {
        set.seed(opt$seed)
        sub_idx <- sample(sub_idx, opt$sr_max_regions_per_cluster_meta)
      }

      mat <- sr_signal_matrix_scale_regions(
        gr = gr_s[sub_idx],
        cov_rlelist = atac_cov,
        flank_bp = opt$sr_flank_bp,
        flank_step = opt$sr_flank_step,
        body_bins = opt$sr_body_bins
      )
      prof <- colMeans(mat, na.rm = TRUE)
      meta_list_all[[length(meta_list_all) + 1L]] <- data.frame(
        sample = ct,
        cluster = cl,
        x = seq_along(prof),
        mean_signal = as.numeric(prof),
        n_regions = length(sub_idx),
        stringsAsFactors = FALSE
      )
    }
  }

  meta_all <- bind_rows(meta_list_all)
  if (nrow(meta_all) > 0) {
    write_tsv(meta_all, file.path(meta_dir, paste0("ATAC_metaplot.k", best_k, ".scaleRegions.tsv")))
    for (ct in unique(meta_all$sample)) {
      df_s <- meta_all[meta_all$sample == ct, , drop = FALSE]
      p <- ggplot(df_s, aes(x = x, y = mean_signal, color = cluster)) +
        geom_line(linewidth = 0.8) +
        geom_vline(xintercept = c(labinfo$body_start - 0.5, labinfo$body_end + 0.5), linetype = "dashed", linewidth = 0.35) +
        theme_bw() +
        labs(
          title = paste0(ct, " ATAC metaplot by k=", best_k, " cluster"),
          x = "Upstream bp / Body scaled bins / Downstream bp",
          y = "Mean ATAC insertion (CPM bw)",
          color = paste0("k=", best_k, " cluster")
        ) +
        scale_x_continuous(breaks = tick_pos, labels = tick_lab)
      pdf(file.path(meta_dir, paste0("ATAC_metaplot.", ct, ".k", best_k, ".scaleRegions.pdf")), width = 9.0, height = 4.8)
      print(p)
      dev.off()
    }

    meta_all$sample <- factor(meta_all$sample, levels = sample_order)
    y_vec <- meta_all$mean_signal
    y_vec <- y_vec[is.finite(y_vec)]
    y_lim <- if (length(y_vec) == 0) c(0, 1) else as.numeric(quantile(y_vec, probs = c(0.001, 0.999), na.rm = TRUE))
    if (!all(is.finite(y_lim)) || y_lim[1] == y_lim[2]) y_lim <- range(y_vec)

    p_all <- ggplot(meta_all, aes(x = x, y = mean_signal, color = cluster)) +
      geom_line(linewidth = 0.7) +
      geom_vline(xintercept = c(labinfo$body_start - 0.5, labinfo$body_end + 0.5), linetype = "dashed", linewidth = 0.35) +
      facet_wrap(~ sample, nrow = 1, drop = FALSE) +
      theme_bw() +
      labs(
        title = paste0("ATAC metaplot by k=", best_k, " cluster"),
        x = "Upstream bp / Body scaled bins / Downstream bp",
        y = "Mean ATAC insertion (CPM bw)",
        color = paste0("k=", best_k, " cluster")
      ) +
      scale_x_continuous(breaks = tick_pos, labels = tick_lab) +
      coord_cartesian(ylim = y_lim)
    pdf(file.path(meta_dir, paste0("ATAC_metaplot.ALLsamples.k", best_k, ".scaleRegions.sharedY.pdf")), width = 12.5, height = 5.0)
    print(p_all)
    dev.off()
  }

  message(">>> [7] Export cluster tables")
  cluster_out_dir <- file.path(opt$out_dir, "cluster_tables")
  dir.create(cluster_out_dir, recursive = TRUE, showWarnings = FALSE)

  df_cluster_assign <- df_scatter %>%
    dplyr::select(any_of(c("sample","region_id","chr","start","end","width_bp")),
                  ATAC_inside, outside_mean, outside_left, outside_right,
                  z1, z2, edge_logratio, cluster)

  write_tsv(df_cluster_assign, file.path(cluster_out_dir, paste0("k", best_k, "_cluster_assignment.z1z2.", fit_space_label, ".tsv")))

  df_cluster_summary <- df_cluster_assign %>%
    dplyr::count(sample, cluster, name = "n_regions") %>%
    group_by(sample) %>%
    mutate(prop = n_regions / sum(n_regions)) %>%
    ungroup() %>%
    arrange(sample, cluster)
  write_tsv(df_cluster_summary, file.path(cluster_out_dir, paste0("k", best_k, "_cluster_summary.by_sample.tsv")))

  df_cluster_proto <- df_cluster_assign %>%
    group_by(sample, cluster) %>%
    summarise(
      n_regions = n(),
      z1_median = median(z1, na.rm = TRUE),
      z2_median = median(z2, na.rm = TRUE),
      ATAC_inside_median = median(ATAC_inside, na.rm = TRUE),
      outside_mean_median = median(outside_mean, na.rm = TRUE),
      edge_logratio_median = median(edge_logratio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(sample, cluster)
  write_tsv(df_cluster_proto, file.path(cluster_out_dir, paste0("k", best_k, "_cluster_feature_summary.tsv")))

  message(">>> [8] Genome-wide scatter colored by cluster")
  plot_dir <- file.path(opt$out_dir, "genome_scatter")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  set.seed(opt$seed)
  df_genome <- df_cluster_assign %>%
    mutate(
      sample = factor(as.character(sample), levels = sample_order),
      chr = as.character(chr),
      mid_bp = floor((as.numeric(start) + as.numeric(end)) / 2),
      mid_mb = mid_bp / 1e6,
      cluster = factor(cluster)
    ) %>%
    filter(chr %in% chr_levels) %>%
    mutate(chr = factor(chr, levels = chr_levels))

  df_genome_plot <- df_genome %>%
    group_by(sample) %>%
    group_modify(~{
      n <- nrow(.x)
      if (n > opt$max_points_per_sample) .x[sample.int(n, opt$max_points_per_sample), , drop = FALSE] else .x
    }) %>%
    ungroup()

  if (nrow(df_genome_plot) > 0) {
    p_gw <- ggplot(df_genome_plot, aes(x = mid_mb, y = chr, color = cluster)) +
      geom_point(size = 0.25, alpha = 0.25) +
      facet_wrap(~ sample, nrow = 1, drop = FALSE) +
      theme_bw() +
      labs(
        title = paste0("Genome-wide region distribution colored by k=", best_k, " cluster"),
        x = "Region midpoint (Mb)",
        y = "Chromosome",
        color = paste0("k=", best_k, " cluster")
      )
    pdf(file.path(plot_dir, paste0("genome_scatter.ALLsamples.chr_vs_midMb.k", best_k, ".pdf")), width = 13.5, height = 7.2)
    print(p_gw)
    dev.off()

    for (ct in levels(df_genome_plot$sample)) {
      df_s <- df_genome_plot[df_genome_plot$sample == ct, , drop = FALSE]
      if (nrow(df_s) == 0) next
      p_s <- ggplot(df_s, aes(x = mid_mb, y = chr, color = cluster)) +
        geom_point(size = 0.28, alpha = 0.28) +
        theme_bw() +
        labs(
          title = paste0(ct, " genome-wide region distribution colored by k=", best_k, " cluster"),
          x = "Region midpoint (Mb)",
          y = "Chromosome",
          color = paste0("k=", best_k, " cluster")
        )
      pdf(file.path(plot_dir, paste0("genome_scatter.", ct, ".chr_vs_midMb.k", best_k, ".pdf")), width = 10.5, height = 7.2)
      print(p_s)
      dev.off()
    }
  }

  message("=== ALL DONE ===")
  message("Output dir: ", opt$out_dir)
  message("Best k source: ", opt$best_k_source)
  message("Best k: ", best_k)
  message("Best-k table: ", file.path(opt$out_dir, "best_k_selection.tsv"))
  message("Sample association: ", txt_out)
}

main()
