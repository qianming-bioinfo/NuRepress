#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

print_help <- function() {
  cat(
"Usage:
  Rscript describe_array_subtype.R --cluster_tsv FILE --out_dir DIR [options]

Required:
  --cluster_tsv FILE          Output table from cluster_array_subtype.R, typically:
                             cluster_tables/k{best_k}_cluster_assignment.z1z2.*.tsv
  --out_dir DIR              Output directory

Optional core:
  --cluster_col STR          Cluster column in cluster_tsv. Default: auto-detect
  --sample_order STR         Comma-separated sample order. Default: order in cluster_tsv
  --eps FLOAT                Small constant for ratio-derived features. Default: 1e-9

Optional modules:
  --region_stats_sheet FILE  TSV with columns: sample, stats_tsv
                             Each stats_tsv is the *.stats.tsv from identify_well_phased_array.py
  --gtf FILE                 GTF/GFF file for TSS annotation module
  --tss_windows STR          Comma-separated windows in bp. Default: 0,250,500,1000,2000,3000,4000,5000
  --use_transcript_tss STR   true|false. Default: true
  --ref_fa FILE              Indexed reference FASTA for GC module (requires .fai)
  --chunk_size INT           GC extraction chunk size. Default: 50000
  --peak_sheet FILE          TSV with columns: sample, peak_bed
                             peak_bed may be bed/narrowPeak/broadPeak/gappedPeak/scoreisland
                             Blank/NA/null means skipping that sample
  --skip_peak_module STR     true|false. Default: false
  --peak_start_is_0based STR true|false. Default: true
  --min_overlap_frac_region FLOAT
                             Peak module threshold. Default: 0.30
  --boundary_frac FLOAT      Peak boundary-zone fraction. Default: 0.10
  --region_start_is_0based STR
                             true|false. Default: true
  --union_long_tsv FILE      Union-region long table for remodeling module
  --union_cluster_col STR    Cluster column in union_long_tsv. Default: auto-detect
  --top_frac FLOAT           Top absolute-change fraction for union subsets. Default: 0.15
  --max_hist_dist_to_tss INT Max distance for TSS histogram. Default: 5000

Optional plotting:
  --plot_width FLOAT         Default: 9
  --plot_height FLOAT        Default: 4.8

Other:
  -h, --help                Show help

Examples:
  Rscript describe_array_subtype.R \\
    --cluster_tsv /path/to/cluster_tables/k3_cluster_assignment.z1z2.within_sample_scaled.tsv \\
    --out_dir /path/to/array_subtype_description

  Rscript describe_array_subtype.R \\
    --cluster_tsv /path/to/cluster_tables/k3_cluster_assignment.z1z2.within_sample_scaled.tsv \\
    --out_dir /path/to/array_subtype_description \\
    --gtf /path/to/annotation.gtf \\
    --ref_fa /path/to/genome.fa \\
    --peak_sheet /path/to/peak_sheet.tsv \\
    --region_stats_sheet /path/to/region_stats_sheet.tsv \\
    --union_long_tsv /path/to/union_long.mapped.tsv

  Rscript describe_array_subtype.R \\
    --cluster_tsv /path/to/cluster_tables/k3_cluster_assignment.z1z2.within_sample_scaled.tsv \\
    --out_dir /path/to/array_subtype_description \\
    --skip_peak_module true

region_stats_sheet format:
  sample\\tstats_tsv
  sample1\\t/path/to/sample1.phased_regions.stats.tsv
  sample2\\t/path/to/sample2.phased_regions.stats.tsv

peak_sheet format:
  sample\\tpeak_bed
  sample1\\t/path/to/sample1_peaks.bed
  sample2\\t/path/to/sample2.scoreisland
  sample3\\tNA
",
    sep = ""
  )
}

parse_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

split_csv <- function(x, mode = c("character", "integer", "numeric")) {
  mode <- match.arg(mode)
  y <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  y <- trimws(y)
  y <- y[nzchar(y)]
  switch(
    mode,
    character = y,
    integer = as.integer(y),
    numeric = as.numeric(y)
  )
}

parse_args <- function(args) {
  normalize_nullable <- function(x) {
    if (is.null(x)) return(NULL)
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(NULL)
    if (tolower(x) %in% c("na", "null", "none")) return(NULL)
    x
  }

  opt <- list(
    cluster_tsv = NULL,
    out_dir = NULL,
    cluster_col = NULL,
    sample_order = NULL,
    eps = 1e-9,
    region_stats_sheet = NULL,
    gtf = NULL,
    tss_windows = c(0L, 250L, 500L, 1000L, 2000L, 3000L, 4000L, 5000L),
    use_transcript_tss = TRUE,
    ref_fa = NULL,
    chunk_size = 50000L,
    peak_sheet = NULL,
    skip_peak_module = FALSE,
    peak_start_is_0based = TRUE,
    min_overlap_frac_region = 0.30,
    boundary_frac = 0.10,
    region_start_is_0based = TRUE,
    union_long_tsv = NULL,
    union_cluster_col = NULL,
    top_frac = 0.15,
    max_hist_dist_to_tss = 5000L,
    plot_width = 9,
    plot_height = 4.8
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
    switch(
      key,
      "--cluster_tsv" = { opt$cluster_tsv <- val; i <- i + 2L },
      "--out_dir" = { opt$out_dir <- val; i <- i + 2L },
      "--cluster_col" = { opt$cluster_col <- val; i <- i + 2L },
      "--sample_order" = { opt$sample_order <- split_csv(val, "character"); i <- i + 2L },
      "--eps" = { opt$eps <- as.numeric(val); i <- i + 2L },
      "--region_stats_sheet" = { opt$region_stats_sheet <- val; i <- i + 2L },
      "--gtf" = { opt$gtf <- val; i <- i + 2L },
      "--tss_windows" = { opt$tss_windows <- split_csv(val, "integer"); i <- i + 2L },
      "--use_transcript_tss" = { opt$use_transcript_tss <- parse_bool(val); i <- i + 2L },
      "--ref_fa" = { opt$ref_fa <- val; i <- i + 2L },
      "--chunk_size" = { opt$chunk_size <- as.integer(val); i <- i + 2L },
      "--peak_sheet" = { opt$peak_sheet <- val; i <- i + 2L },
      "--skip_peak_module" = { opt$skip_peak_module <- parse_bool(val); i <- i + 2L },
      "--peak_start_is_0based" = { opt$peak_start_is_0based <- parse_bool(val); i <- i + 2L },
      "--min_overlap_frac_region" = { opt$min_overlap_frac_region <- as.numeric(val); i <- i + 2L },
      "--boundary_frac" = { opt$boundary_frac <- as.numeric(val); i <- i + 2L },
      "--region_start_is_0based" = { opt$region_start_is_0based <- parse_bool(val); i <- i + 2L },
      "--union_long_tsv" = { opt$union_long_tsv <- val; i <- i + 2L },
      "--union_cluster_col" = { opt$union_cluster_col <- val; i <- i + 2L },
      "--top_frac" = { opt$top_frac <- as.numeric(val); i <- i + 2L },
      "--max_hist_dist_to_tss" = { opt$max_hist_dist_to_tss <- as.integer(val); i <- i + 2L },
      "--plot_width" = { opt$plot_width <- as.numeric(val); i <- i + 2L },
      "--plot_height" = { opt$plot_height <- as.numeric(val); i <- i + 2L },
      stop("Unknown option: ", key)
    )
  }

  if (is.null(opt$cluster_tsv) || is.null(opt$out_dir)) {
    print_help()
    stop("--cluster_tsv and --out_dir are required")
  }

  for (nm in c("region_stats_sheet", "gtf", "ref_fa", "peak_sheet", "union_long_tsv")) {
    opt[[nm]] <- normalize_nullable(opt[[nm]])
  }

  opt
}

stop_if_missing <- function(x) if (!file.exists(x)) stop("File not found: ", x)
ensure_dir <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)

trim_ws_df <- function(df) {
  for (nm in names(df)) {
    if (is.character(df[[nm]])) df[[nm]] <- trimws(df[[nm]])
  }
  df
}

read_sheet_auto <- function(path) {
  stop_if_missing(path)
  if (grepl("\\.(csv)$", path, ignore.case = TRUE)) {
    df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
  trim_ws_df(as.data.frame(df, stringsAsFactors = FALSE))
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "NA")
}

infer_cluster_levels <- function(x) {
  u <- unique(as.character(x))
  u <- u[!is.na(u) & nzchar(u)]
  nums <- suppressWarnings(as.integer(stringr::str_extract(u, "\\d+")))
  if (length(u) == 0) return(character(0))
  if (all(!is.na(nums))) u[order(nums, u)] else sort(u)
}

detect_cluster_col <- function(df, prefer = NULL) {
  if (!is.null(prefer)) {
    if (!prefer %in% colnames(df)) stop("Requested cluster column not found: ", prefer)
    return(prefer)
  }
  cand <- c("cluster", "pattern_cluster", grep("(^|_)cluster$", colnames(df), value = TRUE))
  cand <- unique(cand[cand %in% colnames(df)])
  if (length(cand) == 0) stop("Cannot auto-detect cluster column")
  non_na <- vapply(cand, function(nm) sum(!is.na(df[[nm]]) & nzchar(as.character(df[[nm]]))), numeric(1))
  cand[[which.max(non_na)]]
}

standardize_cluster_df <- function(df, cluster_col, sample_order = NULL, eps = 1e-9) {
  need <- c("sample", "region_id", "chr", "start", "end")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("cluster_tsv missing columns: ", paste(miss, collapse = ", "))

  df <- df %>%
    mutate(
      sample = as.character(sample),
      region_id = as.character(region_id),
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end),
      cluster = as.character(.data[[cluster_col]])
    )

  if (!"width_bp" %in% colnames(df)) {
    df$width_bp <- as.integer(df$end - df$start)
  }
  df$region_len_bp <- ifelse(is.na(df$width_bp), as.integer(df$end - df$start), as.integer(df$width_bp))

  if (!"outside_mean" %in% colnames(df) && all(c("outside_left", "outside_right") %in% colnames(df))) {
    df$outside_mean <- rowMeans(cbind(df$outside_left, df$outside_right), na.rm = TRUE)
  }
  if (!"boundary_mean" %in% colnames(df) && all(c("outside_left", "outside_right") %in% colnames(df))) {
    df$boundary_mean <- rowMeans(cbind(df$outside_left, df$outside_right), na.rm = TRUE)
  }
  if (all(c("outside_left", "outside_right") %in% colnames(df))) {
    df$boundary_balance <- (as.numeric(df$outside_right) - as.numeric(df$outside_left)) /
      (as.numeric(df$outside_right) + as.numeric(df$outside_left) + eps)
    df$boundary_logratio_signed <- log2((as.numeric(df$outside_right) + eps) / (as.numeric(df$outside_left) + eps))
    if (!"z2" %in% colnames(df)) df$z2 <- abs(df$boundary_logratio_signed)
  } else {
    df$boundary_balance <- NA_real_
    df$boundary_logratio_signed <- NA_real_
  }
  if (!"z1" %in% colnames(df) && all(c("ATAC_inside", "outside_mean") %in% colnames(df))) {
    df$z1 <- log2((as.numeric(df$ATAC_inside) + eps) / (as.numeric(df$outside_mean) + eps))
  }

  if (is.null(sample_order)) sample_order <- unique(df$sample)
  sample_order <- unique(sample_order)
  sample_order <- c(sample_order, setdiff(unique(df$sample), sample_order))
  df$sample <- factor(df$sample, levels = sample_order)
  clv <- infer_cluster_levels(df$cluster)
  df$cluster <- factor(df$cluster, levels = clv)
  df
}

write_table <- function(x, path) readr::write_tsv(x, path)

save_boxplot <- function(df, feature, out_dir, sample_levels, cluster_levels, width = 9, height = 4.8) {
  if (!feature %in% colnames(df)) return(invisible(NULL))
  d <- df %>% filter(!is.na(.data[[feature]]))
  if (nrow(d) == 0) return(invisible(NULL))
  d$sample <- factor(as.character(d$sample), levels = sample_levels)
  d$cluster <- factor(as.character(d$cluster), levels = cluster_levels)
  p <- ggplot(d, aes(x = cluster, y = .data[[feature]], fill = cluster)) +
    geom_boxplot(outlier.size = 0.35) +
    facet_wrap(~ sample, nrow = 1, drop = FALSE, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = "Cluster", y = feature, title = paste0(feature, " by cluster"))
  ggsave(file.path(out_dir, paste0(feature, ".boxplot.by_cluster.pdf")), p, width = width, height = height)
}

save_violin <- function(df, feature, out_dir, sample_levels, cluster_levels, width = 9, height = 4.8) {
  if (!feature %in% colnames(df)) return(invisible(NULL))
  d <- df %>% filter(!is.na(.data[[feature]]))
  if (nrow(d) == 0) return(invisible(NULL))
  d$sample <- factor(as.character(d$sample), levels = sample_levels)
  d$cluster <- factor(as.character(d$cluster), levels = cluster_levels)
  p <- ggplot(d, aes(x = cluster, y = .data[[feature]], fill = cluster)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.size = 0.3) +
    facet_wrap(~ sample, nrow = 1, drop = FALSE, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = "Cluster", y = feature, title = paste0(feature, " by cluster"))
  ggsave(file.path(out_dir, paste0(feature, ".violin.by_cluster.pdf")), p, width = width, height = height)
}

summarise_feature_matrix <- function(df, features) {
  features <- intersect(features, colnames(df))
  if (length(features) == 0) return(tibble())
  df %>%
    group_by(sample, cluster) %>%
    summarise(
      n_regions = n(),
      across(all_of(features), list(mean = ~mean(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
}

run_basic_module <- function(df, out_dir, opt) {
  mod_dir <- file.path(out_dir, "01_basic")
  ensure_dir(mod_dir)

  write_table(df, file.path(mod_dir, "cluster_assignment.standardized.tsv"))

  counts <- df %>%
    dplyr::count(sample, cluster, name = "n_regions") %>%
    group_by(sample) %>%
    mutate(prop = n_regions / sum(n_regions)) %>%
    ungroup()
  write_table(counts, file.path(mod_dir, "cluster_counts.by_sample.tsv"))

  feats <- c("region_len_bp", "ATAC_inside", "outside_left", "outside_right", "outside_mean",
             "boundary_mean", "boundary_balance", "boundary_logratio_signed", "edge_logratio", "z1", "z2")
  feat_sum <- summarise_feature_matrix(df, feats)
  write_table(feat_sum, file.path(mod_dir, "base_feature_summary.by_sample_cluster.tsv"))

  p_bar <- ggplot(counts, aes(x = cluster, y = prop, fill = cluster)) +
    geom_col() +
    facet_wrap(~ sample, nrow = 1, drop = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Cluster composition by sample", x = "Cluster", y = "Fraction")
  ggsave(file.path(mod_dir, "cluster_composition.by_sample.pdf"), p_bar, width = opt$plot_width, height = opt$plot_height)

  if (all(c("z1", "z2") %in% colnames(df))) {
    p_sc <- ggplot(df, aes(x = z1, y = z2, color = cluster)) +
      geom_point(size = 0.35, alpha = 0.25) +
      facet_wrap(~ sample, nrow = 1, drop = FALSE) +
      theme_bw() +
      labs(title = "z1/z2 embedding colored by cluster", x = "z1", y = "z2", color = "Cluster")
    ggsave(file.path(mod_dir, "z1z2.scatter.by_sample.pdf"), p_sc, width = 12.5, height = 4.8)
  }

  for (ft in feats) save_boxplot(df, ft, mod_dir, levels(df$sample), levels(df$cluster), opt$plot_width, opt$plot_height)

  invisible(list(counts = counts, feature_summary = feat_sum))
}

read_region_stats_sheet <- function(path) {
  df <- read_sheet_auto(path)
  req <- c("sample", "stats_tsv")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("region_stats_sheet missing columns: ", paste(miss, collapse = ", "))
  df
}

read_one_region_stats <- function(stats_tsv, sample_name) {
  stop_if_missing(stats_tsv)
  x <- readr::read_tsv(stats_tsv, show_col_types = FALSE)
  nm <- names(x)
  if ("region_start" %in% nm && !"start" %in% nm) x$start <- x$region_start
  if ("region_end" %in% nm && !"end" %in% nm) x$end <- x$region_end
  if (!"sample" %in% names(x)) x$sample <- sample_name
  x$sample <- sample_name
  x
}

merge_region_stats <- function(df_cluster, sheet_path, out_dir, opt) {
  mod_dir <- file.path(out_dir, "02_region_phasing_stats")
  ensure_dir(mod_dir)

  sh <- read_region_stats_sheet(sheet_path)
  stat_list <- lapply(seq_len(nrow(sh)), function(i) read_one_region_stats(sh$stats_tsv[[i]], sh$sample[[i]]))
  st <- bind_rows(stat_list)
  st <- trim_ws_df(st)

  stats_keep <- intersect(c("sample", "region_id", "chr", "start", "end", "length_bp", "n_nuc", "mean_spacing",
                            "cv", "good_ratio", "mean_signal", "p20_signal", "min_signal", "median_fuzz",
                            "p80_fuzz", "score", "L", "R"), colnames(st))
  st2 <- st %>% select(all_of(stats_keep))

  if (all(c("sample", "region_id") %in% colnames(st2))) {
    merged <- df_cluster %>% left_join(st2, by = c("sample", "region_id"), suffix = c("", ".stats"))
  } else {
    join_cols <- intersect(c("sample", "chr", "start", "end"), colnames(st2))
    if (length(join_cols) < 4) stop("region stats cannot be joined: need sample+region_id or sample+chr+start+end")
    merged <- df_cluster %>% left_join(st2, by = join_cols, suffix = c("", ".stats"))
  }

  write_table(merged, file.path(mod_dir, "cluster_assignment.with_region_stats.tsv"))

  feats <- c("n_nuc", "mean_spacing", "cv", "good_ratio", "mean_signal", "p20_signal",
             "min_signal", "median_fuzz", "p80_fuzz", "score", "length_bp")
  sum_tbl <- summarise_feature_matrix(merged, feats)
  write_table(sum_tbl, file.path(mod_dir, "region_phasing_summary.by_sample_cluster.tsv"))

  for (ft in feats) save_boxplot(merged, ft, mod_dir, levels(df_cluster$sample), levels(df_cluster$cluster), opt$plot_width, opt$plot_height)
  invisible(merged)
}

require_bioc_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing required packages: ", paste(miss, collapse = ", "))
}

build_id2gene_map_from_gtf <- function(gtf_file, id_type = c("transcript", "gene")) {
  id_type <- match.arg(id_type)
  gtf_gr <- rtracklayer::import(gtf_file)
  mc <- as.data.frame(S4Vectors::mcols(gtf_gr))
  if ("type" %in% colnames(mc)) {
    keep_type <- mc$type %in% c("transcript", "mRNA")
    if (any(keep_type, na.rm = TRUE)) mc <- mc[keep_type, , drop = FALSE]
  }
  get_col <- function(df, nm) if (nm %in% colnames(df)) as.character(df[[nm]]) else NULL
  tx_id <- get_col(mc, "transcript_id")
  gene_id <- get_col(mc, "gene_id")
  gene_name <- get_col(mc, "gene_name")
  if (is.null(gene_name)) gene_name <- gene_id
  if (is.null(gene_name)) gene_name <- rep(NA_character_, nrow(mc))

  first_non_na <- function(x) {
    x <- as.character(x)
    y <- x[!is.na(x) & nzchar(x)]
    if (length(y) == 0) NA_character_ else y[[1]]
  }

  if (id_type == "transcript") {
    if (is.null(tx_id)) stop("GTF lacks transcript_id")
    dfm <- tibble(id = tx_id, gene_name = gene_name, gene_id = gene_id) %>%
      filter(!is.na(id), id != "") %>%
      mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)) %>%
      group_by(id) %>% summarise(gene_name = first_non_na(gene_name), .groups = "drop")
  } else {
    if (is.null(gene_id)) stop("GTF lacks gene_id")
    dfm <- tibble(id = gene_id, gene_name = gene_name) %>%
      filter(!is.na(id), id != "") %>%
      mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", id, gene_name)) %>%
      group_by(id) %>% summarise(gene_name = first_non_na(gene_name), .groups = "drop")
  }
  mp <- dfm$gene_name
  names(mp) <- dfm$id
  mp
}

make_tss_gr <- function(gtf_file, use_transcript = TRUE) {
  txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(
      gtf_file,
      format = ifelse(grepl("\\.gff", gtf_file, ignore.case = TRUE), "gff3", "gtf")
    )
  } else {
    GenomicFeatures::makeTxDbFromGFF(
      gtf_file,
      format = ifelse(grepl("\\.gff", gtf_file, ignore.case = TRUE), "gff3", "gtf")
    )
  }

  if (use_transcript) {
    tx <- GenomicFeatures::transcripts(txdb)
    tss <- GenomicRanges::promoters(tx, upstream = 0, downstream = 1)
    tx_name <- NULL
    if ("tx_name" %in% colnames(S4Vectors::mcols(tx))) tx_name <- as.character(S4Vectors::mcols(tx)$tx_name)
    if (is.null(tx_name) || all(is.na(tx_name))) tx_name <- names(tx)
    if (is.null(tx_name) || all(is.na(tx_name))) tx_name <- paste0("tx_", seq_along(tx))
    S4Vectors::mcols(tss)$tss_id <- tx_name
  } else {
    g <- GenomicFeatures::genes(txdb)
    tss <- GenomicRanges::promoters(g, upstream = 0, downstream = 1)
    gid <- names(g)
    if (is.null(gid) || all(is.na(gid))) {
      if ("gene_id" %in% colnames(S4Vectors::mcols(g))) gid <- as.character(S4Vectors::mcols(g)$gene_id)
    }
    if (is.null(gid) || all(is.na(gid))) gid <- paste0("gene_", seq_along(g))
    S4Vectors::mcols(tss)$tss_id <- gid
  }
  tss
}

calc_tss_assoc <- function(df_cluster, gtf_file, windows, use_transcript_tss, out_dir, opt) {
  require_bioc_pkgs(c("GenomicRanges", "IRanges", "GenomicFeatures", "rtracklayer", "S4Vectors"))
  mod_dir <- file.path(out_dir, "03_TSS_annotation")
  ensure_dir(mod_dir)

  tss_gr <- make_tss_gr(gtf_file, use_transcript_tss)
  id2gene <- build_id2gene_map_from_gtf(gtf_file, if (use_transcript_tss) "transcript" else "gene")

  gr <- GenomicRanges::GRanges(
    seqnames = df_cluster$chr,
    ranges = IRanges::IRanges(start = as.integer(df_cluster$start) + 1L, end = as.integer(df_cluster$end)),
    sample = as.character(df_cluster$sample),
    region_id = as.character(df_cluster$region_id),
    cluster = as.character(df_cluster$cluster)
  )

  near <- GenomicRanges::distanceToNearest(gr, tss_gr, ignore.strand = TRUE)
  dist <- rep(NA_integer_, length(gr))
  nearest_tss_idx <- rep(NA_integer_, length(gr))
  if (length(near) > 0) {
    qh <- S4Vectors::queryHits(near)
    sh <- S4Vectors::subjectHits(near)
    dist[qh] <- S4Vectors::mcols(near)$distance
    nearest_tss_idx[qh] <- sh
  }

  tss_ids <- as.character(S4Vectors::mcols(tss_gr)$tss_id)
  nearest_tss_id <- rep(NA_character_, length(gr))
  ok <- !is.na(nearest_tss_idx) & nearest_tss_idx >= 1 & nearest_tss_idx <= length(tss_ids)
  nearest_tss_id[ok] <- tss_ids[nearest_tss_idx[ok]]
  nearest_gene_name <- unname(id2gene[nearest_tss_id])

  assoc <- tibble(
    sample = as.character(S4Vectors::mcols(gr)$sample),
    region_id = as.character(S4Vectors::mcols(gr)$region_id),
    cluster = as.character(S4Vectors::mcols(gr)$cluster),
    chr = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr) - 1L,
    end = GenomicRanges::end(gr),
    width_bp = GenomicRanges::width(gr),
    dist_to_nearest_TSS = dist,
    nearest_TSS_idx = nearest_tss_idx,
    nearest_TSS_id = nearest_tss_id,
    nearest_TSS_gene_name = nearest_gene_name,
    TSS_in_region = GenomicRanges::countOverlaps(gr, tss_gr, ignore.strand = TRUE) > 0
  )

  for (w in windows) {
    if (w == 0L) {
      assoc[[paste0("overlap_TSSwin_", w)]] <- assoc$TSS_in_region
    } else {
      tss_win <- GenomicRanges::promoters(tss_gr, upstream = w, downstream = w)
      assoc[[paste0("overlap_TSSwin_", w)]] <- GenomicRanges::countOverlaps(gr, tss_win, ignore.strand = TRUE) > 0
    }
  }

  write_table(assoc, file.path(mod_dir, paste0("region_TSS_association.detail.", if (use_transcript_tss) "transcriptTSS" else "geneTSS", ".tsv")))

  gene_list <- assoc %>% filter(TSS_in_region) %>% pull(nearest_TSS_gene_name) %>% unique() %>% stats::na.omit() %>% sort()
  readr::write_lines(gene_list, file.path(mod_dir, paste0("gene_list.TSSinRegion.", if (use_transcript_tss) "transcriptTSS" else "geneTSS", ".txt")))

  sum_sample <- assoc %>%
    group_by(sample) %>%
    summarise(
      n_regions = n(),
      n_TSS_in_region = sum(TSS_in_region, na.rm = TRUE),
      frac_TSS_in_region = n_TSS_in_region / n_regions,
      median_dist = median(dist_to_nearest_TSS, na.rm = TRUE),
      .groups = "drop"
    )
  write_table(sum_sample, file.path(mod_dir, "summary_by_sample.tsv"))

  win_cols <- paste0("overlap_TSSwin_", windows)
  sum_sc <- assoc %>%
    pivot_longer(all_of(win_cols), names_to = "tss_window", values_to = "hit") %>%
    group_by(sample, cluster, tss_window) %>%
    summarise(n = n(), n_hit = sum(hit, na.rm = TRUE), frac_hit = n_hit / n, .groups = "drop")
  write_table(sum_sc, file.path(mod_dir, "summary_by_sample_cluster_TSSwindow.tsv"))

  pA_df <- sum_sc %>%
    group_by(sample, tss_window) %>%
    summarise(frac_hit = sum(n_hit) / sum(n), .groups = "drop") %>%
    mutate(window_bp = as.integer(str_extract(tss_window, "\\d+")))
  pA <- ggplot(pA_df, aes(x = window_bp, y = frac_hit, color = sample)) +
    geom_line(aes(group = sample), linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = sort(unique(pA_df$window_bp))) +
    theme_bw() +
    labs(title = "Fraction of regions associated with TSS", x = "TSS window (bp)", y = "Fraction", color = "Sample")
  ggsave(file.path(mod_dir, "plot_fraction_by_TSSwindow.pdf"), pA, width = 7.5, height = 4.5)

  use_w <- if (2000L %in% windows) 2000L else windows[windows > 0][1]
  col_w <- paste0("overlap_TSSwin_", use_w)
  plot_df <- assoc %>%
    group_by(sample, cluster) %>%
    summarise(n = n(), n_hit = sum(.data[[col_w]], na.rm = TRUE), frac_hit = n_hit / n, .groups = "drop")
  pB <- ggplot(plot_df, aes(x = cluster, y = frac_hit, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = paste0(n_hit, "/", n)), vjust = -0.3, size = 3) +
    facet_wrap(~ sample, nrow = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste0("TSS association by cluster (window=", use_w, " bp)"), x = "Cluster", y = "Fraction")
  ggsave(file.path(mod_dir, paste0("plot_by_cluster_TSSwindow_", use_w, ".pdf")), pB, width = 10, height = 4)

  pC_df <- assoc %>% filter(is.finite(dist_to_nearest_TSS), dist_to_nearest_TSS >= 0, dist_to_nearest_TSS <= opt$max_hist_dist_to_tss)
  pC <- ggplot(pC_df, aes(x = dist_to_nearest_TSS)) +
    geom_histogram(binwidth = max(50, floor(opt$max_hist_dist_to_tss / 50)), boundary = 0, closed = "left") +
    facet_wrap(~ sample, scales = "free_y", nrow = 1) +
    theme_bw() +
    labs(title = paste0("Distance to nearest TSS [0-", opt$max_hist_dist_to_tss, " bp]"), x = "Distance to nearest TSS (bp)", y = "Count")
  ggsave(file.path(mod_dir, "plot_dist_to_nearest_TSS.hist.pdf"), pC, width = 10, height = 4)

  invisible(assoc)
}

fix_chr_onebyone <- function(x, fa_seqnames) {
  x <- gsub("\\s+", "", x)
  out <- x
  ok <- out %in% fa_seqnames
  need <- !ok
  out[need] <- ifelse(grepl("^chr", out[need]), out[need], paste0("chr", out[need]))
  ok <- out %in% fa_seqnames
  need <- !ok
  out[need] <- sub("^chr", "", out[need])
  ok <- out %in% fa_seqnames
  out[!ok] <- NA_character_
  out
}

calc_gc_from_fasta <- function(df_cluster, ref_fa, out_dir, opt) {
  require_bioc_pkgs(c("GenomicRanges", "IRanges", "Rsamtools", "Biostrings"))
  mod_dir <- file.path(out_dir, "04_GC_content")
  ensure_dir(mod_dir)
  stop_if_missing(ref_fa)
  fai_path <- paste0(ref_fa, ".fai")
  stop_if_missing(fai_path)

  fai <- readr::read_tsv(
    fai_path,
    col_names = c("seqname", "seqlength", "offset", "line_bases", "line_width"),
    col_types = "cicii",
    show_col_types = FALSE
  )
  fa_seqnames <- fai$seqname
  fa_seqlen <- setNames(as.integer(fai$seqlength), fai$seqname)

  df <- df_cluster %>%
    mutate(
      chr_h = fix_chr_onebyone(chr, fa_seqnames),
      start_1 = as.integer(start) + 1L,
      end_1 = as.integer(end)
    ) %>%
    filter(!is.na(chr_h)) %>%
    mutate(chr_len = fa_seqlen[chr_h]) %>%
    filter(!is.na(chr_len), start_1 >= 1L, end_1 <= chr_len, start_1 <= end_1)

  gr <- GenomicRanges::GRanges(seqnames = df$chr_h, ranges = IRanges::IRanges(start = df$start_1, end = df$end_1))
  fa <- Rsamtools::FaFile(ref_fa)
  Rsamtools::open.FaFile(fa)
  on.exit(Rsamtools::close.FaFile(fa), add = TRUE)

  n <- length(gr)
  idx_starts <- seq(1L, n, by = opt$chunk_size)
  gc_vec <- rep(NA_real_, n)

  letter_gc <- function(dna) {
    g <- Biostrings::letterFrequency(dna, letters = "G", as.prob = FALSE)
    c_ <- Biostrings::letterFrequency(dna, letters = "C", as.prob = FALSE)
    n_ <- Biostrings::letterFrequency(dna, letters = "N", as.prob = FALSE)
    len <- Biostrings::width(dna)
    eff <- pmax(len - n_, 0L)
    (g + c_) / ifelse(eff == 0, NA_real_, eff)
  }

  for (i in seq_along(idx_starts)) {
    s <- idx_starts[i]
    e <- min(s + opt$chunk_size - 1L, n)
    dna <- tryCatch(Biostrings::getSeq(fa, gr[s:e]), error = function(e) NULL)
    if (is.null(dna)) {
      for (j in s:e) {
        one <- tryCatch(Biostrings::getSeq(fa, gr[j]), error = function(e) NULL)
        gc_vec[j] <- if (is.null(one)) NA_real_ else as.numeric(letter_gc(one))
      }
    } else {
      gc_vec[s:e] <- as.numeric(letter_gc(dna))
    }
  }

  df_out <- df %>% mutate(GC_frac = gc_vec, GC_pct = 100 * GC_frac)
  write_table(df_out, file.path(mod_dir, "cluster_assignment.with_GC.tsv"))

  mean_gc <- df_out %>% filter(!is.na(GC_frac)) %>% group_by(sample, cluster) %>%
    summarise(n_regions = n(), mean_GC = mean(GC_frac), mean_GC_pct = mean(GC_pct), .groups = "drop")
  write_table(mean_gc, file.path(mod_dir, "mean_GC_by_sample_and_cluster.tsv"))

  save_boxplot(df_out, "GC_pct", mod_dir, levels(df_cluster$sample), levels(df_cluster$cluster), opt$plot_width, opt$plot_height)
  save_boxplot(df_out, "region_len_bp", mod_dir, levels(df_cluster$sample), levels(df_cluster$cluster), opt$plot_width, opt$plot_height)
  invisible(df_out)
}

read_peak_like_table <- function(path, peak_start_is_0based = TRUE) {
  x <- tryCatch(
    readr::read_tsv(
      path,
      col_names = FALSE,
      comment = "#",
      show_col_types = FALSE,
      progress = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(x)) {
    stop("Cannot read peak-like file as tabular text: ", path)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (ncol(x) < 3) {
    stop("Peak-like file has fewer than 3 columns: ", path)
  }

  colnames(x)[1:3] <- c("chr", "start", "end")

  chr <- trimws(as.character(x$chr))
  st <- suppressWarnings(as.integer(as.character(x$start)))
  ed <- suppressWarnings(as.integer(as.character(x$end)))

  keep <- !is.na(chr) & nzchar(chr) & !is.na(st) & !is.na(ed)
  x <- x[keep, , drop = FALSE]
  if (nrow(x) == 0L) return(GenomicRanges::GRanges())

  st <- suppressWarnings(as.integer(as.character(x$start)))
  ed <- suppressWarnings(as.integer(as.character(x$end)))
  if (peak_start_is_0based) st <- st + 1L

  keep2 <- !is.na(st) & !is.na(ed) & st <= ed
  if (!any(keep2)) return(GenomicRanges::GRanges())

  gr <- GenomicRanges::GRanges(
    seqnames = as.character(x$chr[keep2]),
    ranges = IRanges::IRanges(start = st[keep2], end = ed[keep2])
  )

  extra_cols <- x[keep2, setdiff(colnames(x), c("chr", "start", "end")), drop = FALSE]
  if (ncol(extra_cols) > 0) {
    extra_cols[] <- lapply(extra_cols, function(v) {
      if (is.character(v)) trimws(v) else v
    })
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(extra_cols)
  }

  gr
}

import_peak_any <- function(path, peak_start_is_0based = TRUE) {
  if (is.null(path)) return(GenomicRanges::GRanges())
  path <- trimws(as.character(path))
  if (!nzchar(path) || tolower(path) %in% c("na", "null", "none")) {
    return(GenomicRanges::GRanges())
  }
  stop_if_missing(path)

  path_low <- tolower(path)

  bed_like_ext <- grepl("\\.(bed|narrowpeak|broadpeak|gappedpeak)(\\.gz)?$", path_low)
  scoreisland_ext <- grepl("\\.scoreisland(\\.gz)?$", path_low)

  if (scoreisland_ext) {
    gr <- tryCatch(
      rtracklayer::import(path, format = "BED"),
      error = function(e) NULL
    )
    if (!is.null(gr)) return(gr)

    return(read_peak_like_table(path, peak_start_is_0based = peak_start_is_0based))
  }

  if (bed_like_ext) {
    gr <- tryCatch(
      rtracklayer::import(path, format = "BED"),
      error = function(e) NULL
    )
    if (!is.null(gr)) return(gr)

    return(read_peak_like_table(path, peak_start_is_0based = peak_start_is_0based))
  }

  gr <- tryCatch(
    rtracklayer::import(path),
    error = function(e) NULL
  )
  if (!is.null(gr)) return(gr)

  gr <- tryCatch(
    rtracklayer::import(path, format = "BED"),
    error = function(e) NULL
  )
  if (!is.null(gr)) return(gr)

  read_peak_like_table(path, peak_start_is_0based = peak_start_is_0based)
}

read_peak_sheet <- function(path) {
  df <- read_sheet_auto(path)

  if (!"sample" %in% colnames(df)) {
    stop("peak_sheet missing required column: sample")
  }

  peak_col <- NULL
  if ("peak_bed" %in% colnames(df)) {
    peak_col <- "peak_bed"
  } else if ("peak_file" %in% colnames(df)) {
    peak_col <- "peak_file"
  } else {
    stop("peak_sheet must contain column 'peak_bed' or 'peak_file'")
  }

  df <- df %>%
    mutate(
      sample = trimws(as.character(sample)),
      peak_bed = trimws(as.character(.data[[peak_col]]))
    ) %>%
    select(sample, peak_bed)

  df$peak_bed[!nzchar(df$peak_bed)] <- NA_character_
  df$peak_bed[tolower(df$peak_bed) %in% c("na", "null", "none")] <- NA_character_

  bad <- !is.na(df$peak_bed) & !file.exists(df$peak_bed)
  if (any(bad)) {
    stop(
      "Some peak files in peak_sheet do not exist:\n",
      paste0("  - ", df$sample[bad], ": ", df$peak_bed[bad], collapse = "\n")
    )
  }

  df
}

annotate_peak_positions_one_sample <- function(df_sample,
                                               peaks_bed_path,
                                               sample_name,
                                               region_start_is_0based = TRUE,
                                               boundary_frac = 0.10,
                                               peak_start_is_0based = TRUE) {
  gr_reg <- GenomicRanges::GRanges(
    seqnames = df_sample$chr,
    ranges = IRanges::IRanges(
      start = if (region_start_is_0based) as.integer(df_sample$start) + 1L else as.integer(df_sample$start),
      end = as.integer(df_sample$end)
    )
  )

  gr_peak <- import_peak_any(peaks_bed_path, peak_start_is_0based = peak_start_is_0based)

  if (length(gr_peak) == 0L) {
    return(df_sample %>% mutate(
      sample = sample_name,
      peak_overlap_bp = 0L,
      overlap_frac_region = NA_real_,
      best_peak_start = NA_integer_,
      best_peak_end = NA_integer_,
      peak_width_bp = NA_integer_,
      center_norm_raw = NA_real_,
      center_norm = NA_real_,
      center_in_peak = FALSE,
      pos_start_norm_raw = NA_real_,
      pos_end_norm_raw = NA_real_,
      pos_start_norm = NA_real_,
      pos_end_norm = NA_real_,
      edge_proximity = NA_real_,
      edge_touch_norm = NA_real_,
      in_boundary_zone = FALSE,
      has_peak_overlap = FALSE
    ))
  }

  hits <- GenomicRanges::findOverlaps(gr_reg, gr_peak, ignore.strand = TRUE)
  n_reg <- length(gr_reg)
  best_peak_idx <- rep(NA_integer_, n_reg)
  best_ov_bp <- rep(0L, n_reg)

  if (length(hits) > 0L) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    ov_ir <- IRanges::pintersect(GenomicRanges::ranges(gr_reg[qh]), GenomicRanges::ranges(gr_peak[sh]))
    ov_bp <- IRanges::width(ov_ir)

    best <- tibble(region_idx = qh, peak_idx = sh, ov_bp = ov_bp) %>%
      group_by(region_idx) %>%
      slice_max(order_by = ov_bp, n = 1, with_ties = FALSE) %>%
      ungroup()

    best_peak_idx[best$region_idx] <- best$peak_idx
    best_ov_bp[best$region_idx] <- best$ov_bp
  }

  reg_start <- GenomicRanges::start(gr_reg)
  reg_end <- GenomicRanges::end(gr_reg)
  reg_width <- GenomicRanges::width(gr_reg)
  reg_center <- as.integer(floor((reg_start + reg_end) / 2))
  ok <- !is.na(best_peak_idx)

  peak_start <- rep(NA_integer_, n_reg)
  peak_end <- rep(NA_integer_, n_reg)
  peak_width <- rep(NA_integer_, n_reg)
  peak_start[ok] <- GenomicRanges::start(gr_peak)[best_peak_idx[ok]]
  peak_end[ok] <- GenomicRanges::end(gr_peak)[best_peak_idx[ok]]
  peak_width[ok] <- GenomicRanges::width(gr_peak)[best_peak_idx[ok]]

  denom <- rep(NA_real_, n_reg)
  denom[ok] <- pmax((peak_end[ok] - peak_start[ok]), 1)

  clamp01 <- function(x) pmin(pmax(x, 0), 1)

  center_norm_raw <- rep(NA_real_, n_reg)
  pos_start_norm_raw <- rep(NA_real_, n_reg)
  pos_end_norm_raw <- rep(NA_real_, n_reg)

  center_norm_raw[ok] <- (reg_center[ok] - peak_start[ok]) / denom[ok]
  pos_start_norm_raw[ok] <- (reg_start[ok] - peak_start[ok]) / denom[ok]
  pos_end_norm_raw[ok] <- (reg_end[ok] - peak_start[ok]) / denom[ok]

  center_in_peak <- !is.na(center_norm_raw) & center_norm_raw >= 0 & center_norm_raw <= 1

  center_norm <- pos_start_norm <- pos_end_norm <- rep(NA_real_, n_reg)
  center_norm[ok] <- clamp01(center_norm_raw[ok])
  pos_start_norm[ok] <- clamp01(pos_start_norm_raw[ok])
  pos_end_norm[ok] <- clamp01(pos_end_norm_raw[ok])

  edge_proximity <- rep(NA_real_, n_reg)
  edge_proximity[ok] <- pmin(center_norm[ok], 1 - center_norm[ok])

  edge_touch_norm <- rep(NA_real_, n_reg)
  edge_touch_norm[ok] <- pmin(pos_start_norm[ok], 1 - pos_end_norm[ok])
  edge_touch_norm[ok] <- clamp01(edge_touch_norm[ok])

  overlap_frac_region <- rep(NA_real_, n_reg)
  overlap_frac_region[ok] <- best_ov_bp[ok] / reg_width[ok]

  in_boundary_zone <- center_in_peak & (
    center_norm_raw <= boundary_frac | center_norm_raw >= (1 - boundary_frac)
  )

  df_sample %>% mutate(
    sample = sample_name,
    peak_overlap_bp = best_ov_bp,
    overlap_frac_region = overlap_frac_region,
    best_peak_start = peak_start,
    best_peak_end = peak_end,
    peak_width_bp = peak_width,
    center_norm_raw = center_norm_raw,
    center_norm = center_norm,
    center_in_peak = center_in_peak,
    pos_start_norm_raw = pos_start_norm_raw,
    pos_end_norm_raw = pos_end_norm_raw,
    pos_start_norm = pos_start_norm,
    pos_end_norm = pos_end_norm,
    edge_proximity = edge_proximity,
    edge_touch_norm = edge_touch_norm,
    in_boundary_zone = in_boundary_zone,
    has_peak_overlap = ok & (best_ov_bp > 0)
  )
}

run_peak_module <- function(df_cluster, peak_sheet, out_dir, opt) {
  require_bioc_pkgs(c("GenomicRanges", "IRanges", "rtracklayer", "S4Vectors"))

  if (isTRUE(opt$skip_peak_module)) {
    message(">>> Peak module skipped by --skip_peak_module true")
    return(invisible(NULL))
  }
  if (is.null(peak_sheet) || !nzchar(trimws(peak_sheet))) {
    message(">>> Peak module skipped because --peak_sheet is empty")
    return(invisible(NULL))
  }

  mod_dir <- file.path(out_dir, "05_peak_relative_position")
  ensure_dir(mod_dir)

  sh <- read_peak_sheet(peak_sheet)
  if (nrow(sh) == 0L) {
    message(">>> Peak module skipped because peak_sheet has no rows")
    return(invisible(NULL))
  }

  sample_order <- levels(df_cluster$sample)
  res_list <- list()

  for (s in sample_order) {
    df_s <- df_cluster %>% filter(as.character(sample) == s)
    if (nrow(df_s) == 0) next

    row_s <- sh[sh$sample == s, , drop = FALSE]
    if (nrow(row_s) == 0) {
      message(">>> Peak module: sample skipped because not found in peak_sheet: ", s)
      next
    }

    peak_path <- row_s$peak_bed[[1]]
    if (is.na(peak_path) || !nzchar(trimws(peak_path))) {
      message(">>> Peak module: sample skipped because peak path is empty: ", s)
      next
    }

    res_list[[s]] <- annotate_peak_positions_one_sample(
      df_s,
      peak_path,
      s,
      region_start_is_0based = opt$region_start_is_0based,
      boundary_frac = opt$boundary_frac,
      peak_start_is_0based = opt$peak_start_is_0based
    )
  }

  df_pos <- bind_rows(res_list)
  if (nrow(df_pos) == 0) {
    message(">>> Peak module finished with no valid samples to process")
    return(invisible(NULL))
  }

  df_overlap <- df_pos %>%
    filter(has_peak_overlap, !is.na(overlap_frac_region), overlap_frac_region >= opt$min_overlap_frac_region)

  df_center <- df_overlap %>%
    filter(center_in_peak, !is.na(center_norm_raw))

  df_edge <- df_overlap %>%
    filter(!is.na(edge_touch_norm), !is.na(cluster))

  write_table(df_pos, file.path(mod_dir, "region_peak_position.all.tsv"))
  write_table(df_overlap, file.path(mod_dir, paste0("region_peak_position.overlapFracGe", opt$min_overlap_frac_region, ".tsv")))
  write_table(df_center, file.path(mod_dir, paste0("region_peak_position.centerInPeak.overlapFracGe", opt$min_overlap_frac_region, ".tsv")))

  df_summary <- df_pos %>%
    mutate(
      pass_overlap = !is.na(overlap_frac_region) & overlap_frac_region >= opt$min_overlap_frac_region,
      use_overlap = has_peak_overlap & pass_overlap,
      use_center = use_overlap & center_in_peak
    ) %>%
    group_by(sample, cluster) %>%
    summarise(
      n_regions_total = n(),
      n_overlap = sum(use_overlap, na.rm = TRUE),
      frac_overlap = n_overlap / n_regions_total,
      n_center_in_peak = sum(use_center, na.rm = TRUE),
      frac_center_in_peak_among_overlap = ifelse(n_overlap > 0, n_center_in_peak / n_overlap, NA_real_),
      median_center_norm = median(center_norm_raw[use_center], na.rm = TRUE),
      median_edge_touch_norm = median(edge_touch_norm[use_overlap], na.rm = TRUE),
      n_boundary_zone = sum(in_boundary_zone & use_center, na.rm = TRUE),
      frac_boundary_zone = ifelse(n_center_in_peak > 0, n_boundary_zone / n_center_in_peak, NA_real_),
      .groups = "drop"
    )
  write_table(df_summary, file.path(mod_dir, paste0("summary_by_sample_cluster.overlapFracGe", opt$min_overlap_frac_region, ".boundary", opt$boundary_frac, ".tsv")))

  if (length(levels(df_cluster$cluster)) >= 2) {
    c1 <- levels(df_cluster$cluster)[1]
    c2 <- levels(df_cluster$cluster)[2]

    df_wilcox <- df_edge %>%
      group_by(sample) %>%
      summarise(
        n_C1 = sum(cluster == c1),
        n_C2 = sum(cluster == c2),
        p_value = {
          x <- edge_touch_norm[cluster == c1]
          y <- edge_touch_norm[cluster == c2]
          if (length(x) >= 5 && length(y) >= 5) wilcox.test(x, y)$p.value else NA_real_
        },
        median_C1 = median(edge_touch_norm[cluster == c1], na.rm = TRUE),
        median_C2 = median(edge_touch_norm[cluster == c2], na.rm = TRUE),
        .groups = "drop"
      )
    write_table(df_wilcox, file.path(mod_dir, paste0("wilcox_cluster1_vs_cluster2_edgeTouch.overlapFracGe", opt$min_overlap_frac_region, ".tsv")))

    if (nrow(df_center) > 0) {
      df_fisher <- df_center %>%
        group_by(sample) %>%
        summarise(
          A_boundary = sum(cluster == c1 & in_boundary_zone),
          A_nonboundary = sum(cluster == c1 & !in_boundary_zone),
          B_boundary = sum(cluster == c2 & in_boundary_zone),
          B_nonboundary = sum(cluster == c2 & !in_boundary_zone),
          p_value = {
            m <- matrix(c(A_boundary, A_nonboundary, B_boundary, B_nonboundary), nrow = 2, byrow = TRUE)
            if (all(m >= 0) && sum(m) > 0) fisher.test(m)$p.value else NA_real_
          },
          .groups = "drop"
        )
      write_table(df_fisher, file.path(mod_dir, paste0("fisher_cluster1_vs_cluster2_boundaryZone.centerInPeak.overlapFracGe", opt$min_overlap_frac_region, ".boundary", opt$boundary_frac, ".tsv")))
    }
  }

  if (nrow(df_center) > 0) {
    p_center <- ggplot(df_center, aes(x = center_norm_raw, color = cluster)) +
      geom_density(linewidth = 0.8, adjust = 1) +
      facet_wrap(~ sample, nrow = 1) +
      theme_bw() +
      labs(title = "Center position of regions within peaks", x = "center_norm_raw", y = "Density")
    ggsave(file.path(mod_dir, paste0("density_center_norm_raw.centerInPeak.overlapFracGe", opt$min_overlap_frac_region, ".pdf")), p_center, width = 12, height = 4)
  }

  if (nrow(df_edge) > 0) {
    p_edge <- ggplot(df_edge, aes(x = cluster, y = edge_touch_norm, fill = cluster)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.2, outlier.size = 0.35) +
      facet_wrap(~ sample, nrow = 1, scales = "free_y") +
      theme_bw() +
      theme(legend.position = "none") +
      labs(title = "Edge-touch of regions within peaks", x = "Cluster", y = "edge_touch_norm")
    ggsave(file.path(mod_dir, paste0("violin_edge_touch_norm.overlapFracGe", opt$min_overlap_frac_region, ".pdf")), p_edge, width = 12, height = 4)
  }

  invisible(df_pos)
}

write_bed4 <- function(df, out_bed, id_col = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  nm <- if (is.null(id_col) || !id_col %in% colnames(df)) rep(".", nrow(df)) else as.character(df[[id_col]])
  readr::write_tsv(tibble(chr = df$chr, start = as.integer(df$start), end = as.integer(df$end), name = nm), out_bed, col_names = FALSE)
}

top_abs_subset <- function(df, delta_col, frac = 0.15) {
  if (nrow(df) == 0 || !delta_col %in% colnames(df)) return(df[0, , drop = FALSE])
  df2 <- df %>% filter(!is.na(.data[[delta_col]])) %>% mutate(abs_delta = abs(.data[[delta_col]])) %>% arrange(desc(abs_delta))
  if (nrow(df2) == 0) return(df2)
  k <- max(1L, ceiling(nrow(df2) * frac))
  dplyr::slice_head(df2, n = k)
}

run_union_module <- function(path, sample_order, cluster_col, out_dir, opt) {
  mod_dir <- file.path(out_dir, "06_union_region_remodeling")
  ensure_dir(mod_dir)
  df <- readr::read_tsv(path, show_col_types = FALSE)
  cluster_col <- detect_cluster_col(df, cluster_col)

  req <- c("union_id", "sample", "chr", "start", "end")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("union_long_tsv missing columns: ", paste(miss, collapse = ", "))

  if (!"width_bp" %in% colnames(df)) df$width_bp <- as.integer(df$end - df$start)
  if (!"is_phased" %in% colnames(df)) stop("union_long_tsv must contain is_phased")
  if (!"edge_logratio" %in% colnames(df)) stop("union_long_tsv must contain edge_logratio")
  if (!"boundary_mean" %in% colnames(df) && all(c("ATAC_left1k", "ATAC_right1k") %in% colnames(df))) {
    df$boundary_mean <- rowMeans(cbind(df$ATAC_left1k, df$ATAC_right1k), na.rm = TRUE)
  }
  df$pattern_cluster <- as.character(df[[cluster_col]])
  df <- df %>% filter(sample %in% sample_order) %>% mutate(sample = factor(sample, levels = sample_order))

  if (length(sample_order) != 3L) {
    stop("union_long_tsv module currently requires exactly 3 samples in --sample_order")
  }

  safe <- safe_name(sample_order)
  names(safe) <- sample_order
  s1 <- sample_order[[1]]; s2 <- sample_order[[2]]; s3 <- sample_order[[3]]
  a1 <- safe[[s1]]; a2 <- safe[[s2]]; a3 <- safe[[s3]]

  metric_cols <- intersect(c("is_phased", "pattern_cluster", "edge_logratio", "edge_index",
                             "ATAC_inside", "ATAC_left1k", "ATAC_right1k", "boundary_mean"), colnames(df))
  core_cols <- c("union_id", "chr", "start", "end", "width_bp")
  if ("location" %in% colnames(df)) core_cols <- c(core_cols, "location")

  w <- df %>%
    mutate(sample_safe = safe_name(as.character(sample))) %>%
    select(all_of(core_cols), sample_safe, all_of(metric_cols)) %>%
    pivot_wider(
      id_cols = all_of(core_cols),
      names_from = sample_safe,
      values_from = all_of(metric_cols),
      names_glue = "{.value}__{sample_safe}"
    )

  p1 <- w[[paste0("is_phased__", a1)]] %in% TRUE
  p2 <- w[[paste0("is_phased__", a2)]] %in% TRUE
  p3 <- w[[paste0("is_phased__", a3)]] %in% TRUE
  c1 <- paste0("pattern_cluster__", a1); c2 <- paste0("pattern_cluster__", a2); c3 <- paste0("pattern_cluster__", a3)
  e1 <- paste0("edge_logratio__", a1); e2 <- paste0("edge_logratio__", a2); e3 <- paste0("edge_logratio__", a3)

  w <- w %>% mutate(
    presence_binary = paste0(as.integer(p1), as.integer(p2), as.integer(p3)),
    presence_label = case_when(
      p1 & p2 & p3 ~ paste(a1, a2, a3, sep = "_"),
      p1 & !p2 & !p3 ~ a1,
      !p1 & p2 & !p3 ~ a2,
      !p1 & !p2 & p3 ~ a3,
      p1 & p2 & !p3 ~ paste(a1, a2, "only", sep = "_"),
      !p1 & p2 & p3 ~ paste(a2, a3, "only", sep = "_"),
      p1 & !p2 & p3 ~ paste(a1, a3, "only", sep = "_"),
      TRUE ~ "none"
    ),
    state_1 = ifelse(p1, .data[[c1]], NA_character_),
    state_2 = ifelse(p2, .data[[c2]], NA_character_),
    state_3 = ifelse(p3, .data[[c3]], NA_character_),
    d12 = ifelse(p1 & p2, .data[[e2]] - .data[[e1]], NA_real_),
    d23 = ifelse(p2 & p3, .data[[e3]] - .data[[e2]], NA_real_),
    d13 = ifelse(p1 & p3, .data[[e3]] - .data[[e1]], NA_real_),
    switch_12 = ifelse(p1 & p2 & !is.na(state_1) & !is.na(state_2) & state_1 != state_2, paste0(state_1, "_to_", state_2), NA_character_),
    switch_23 = ifelse(p2 & p3 & !is.na(state_2) & !is.na(state_3) & state_2 != state_3, paste0(state_2, "_to_", state_3), NA_character_),
    switch_13 = ifelse(p1 & p3 & !is.na(state_1) & !is.na(state_3) & state_1 != state_3, paste0(state_1, "_to_", state_3), NA_character_),
    mono_123 = ifelse(
      p1 & p2 & p3,
      ifelse(.data[[e1]] < .data[[e2]] & .data[[e2]] < .data[[e3]], "increasing",
             ifelse(.data[[e1]] > .data[[e2]] & .data[[e2]] > .data[[e3]], "decreasing", "non_mono")),
      NA_character_
    )
  )

  write_table(w, file.path(mod_dir, "union_long.wide_with_remodeling_labels.tsv"))

  subsets <- list(
    shared_all = w %>% filter(p1 & p2 & p3),
    only_1 = w %>% filter(p1 & !p2 & !p3),
    only_2 = w %>% filter(!p1 & p2 & !p3),
    only_3 = w %>% filter(!p1 & !p2 & p3),
    pair_12_only = w %>% filter(p1 & p2 & !p3),
    pair_23_only = w %>% filter(!p1 & p2 & p3),
    pair_13_only = w %>% filter(p1 & !p2 & p3),
    switch_12 = w %>% filter(!is.na(switch_12)),
    switch_23 = w %>% filter(!is.na(switch_23)),
    switch_13 = w %>% filter(!is.na(switch_13)),
    mono_increasing = w %>% filter(mono_123 == "increasing"),
    mono_decreasing = w %>% filter(mono_123 == "decreasing"),
    d12_up = w %>% filter(!is.na(d12) & d12 > 0),
    d12_down = w %>% filter(!is.na(d12) & d12 < 0),
    d23_up = w %>% filter(!is.na(d23) & d23 > 0),
    d23_down = w %>% filter(!is.na(d23) & d23 < 0),
    d13_up = w %>% filter(!is.na(d13) & d13 > 0),
    d13_down = w %>% filter(!is.na(d13) & d13 < 0),
    top_abs_d12 = top_abs_subset(w %>% filter(p1 & p2), "d12", opt$top_frac),
    top_abs_d23 = top_abs_subset(w %>% filter(p2 & p3), "d23", opt$top_frac),
    top_abs_d13 = top_abs_subset(w %>% filter(p1 & p3), "d13", opt$top_frac)
  )

  for (nm in names(subsets)) {
    write_bed4(subsets[[nm]], file.path(mod_dir, paste0("union.", nm, ".bed")), "union_id")
  }

  summary_tbl <- tibble(
    subset = names(subsets),
    n_union = vapply(subsets, nrow, integer(1))
  )
  write_table(summary_tbl, file.path(mod_dir, "union_subset_summary.tsv"))
  invisible(w)
}

write_metadata <- function(opt, out_dir, detected_cluster_col, sample_order) {
  lines <- c(
    paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("cluster_tsv=", normalizePath(opt$cluster_tsv, mustWork = FALSE)),
    paste0("out_dir=", normalizePath(out_dir, mustWork = FALSE)),
    paste0("cluster_col=", detected_cluster_col),
    paste0("sample_order=", paste(sample_order, collapse = ",")),
    paste0("region_stats_sheet=", ifelse(is.null(opt$region_stats_sheet), "NA", normalizePath(opt$region_stats_sheet, mustWork = FALSE))),
    paste0("gtf=", ifelse(is.null(opt$gtf), "NA", normalizePath(opt$gtf, mustWork = FALSE))),
    paste0("ref_fa=", ifelse(is.null(opt$ref_fa), "NA", normalizePath(opt$ref_fa, mustWork = FALSE))),
    paste0("peak_sheet=", ifelse(is.null(opt$peak_sheet), "NA", normalizePath(opt$peak_sheet, mustWork = FALSE))),
    paste0("skip_peak_module=", opt$skip_peak_module),
    paste0("peak_start_is_0based=", opt$peak_start_is_0based),
    paste0("union_long_tsv=", ifelse(is.null(opt$union_long_tsv), "NA", normalizePath(opt$union_long_tsv, mustWork = FALSE)))
  )
  readr::write_lines(lines, file.path(out_dir, "run_metadata.txt"))
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  ensure_dir(opt$out_dir)
  stop_if_missing(opt$cluster_tsv)

  df_raw <- readr::read_tsv(opt$cluster_tsv, show_col_types = FALSE)
  cluster_col <- detect_cluster_col(df_raw, opt$cluster_col)
  sample_order <- opt$sample_order
  if (is.null(sample_order)) sample_order <- unique(as.character(df_raw$sample))

  df_cluster <- standardize_cluster_df(df_raw, cluster_col, sample_order, opt$eps)
  write_metadata(opt, opt$out_dir, cluster_col, levels(df_cluster$sample))

  message(">>> [1] Basic subtype description")
  run_basic_module(df_cluster, opt$out_dir, opt)

  if (!is.null(opt$region_stats_sheet)) {
    message(">>> [2] Region phasing-stat description")
    merge_region_stats(df_cluster, opt$region_stats_sheet, opt$out_dir, opt)
  }

  if (!is.null(opt$gtf)) {
    message(">>> [3] TSS annotation")
    calc_tss_assoc(df_cluster, opt$gtf, opt$tss_windows, opt$use_transcript_tss, opt$out_dir, opt)
  }

  if (!is.null(opt$ref_fa)) {
    message(">>> [4] GC content")
    calc_gc_from_fasta(df_cluster, opt$ref_fa, opt$out_dir, opt)
  }

  if (isTRUE(opt$skip_peak_module)) {
    message(">>> [5] Peak-relative position [SKIPPED]")
  } else if (!is.null(opt$peak_sheet)) {
    message(">>> [5] Peak-relative position")
    run_peak_module(df_cluster, opt$peak_sheet, opt$out_dir, opt)
  }

  if (!is.null(opt$union_long_tsv)) {
    message(">>> [6] Union-region remodeling")
    stop_if_missing(opt$union_long_tsv)
    run_union_module(opt$union_long_tsv, levels(df_cluster$sample), opt$union_cluster_col, opt$out_dir, opt)
  }

  message("=== ALL DONE ===")
  message("Output dir: ", opt$out_dir)
}

main()