#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

usage_text <- function() {
  paste(
    "Usage:",
    "  Rscript score_repressor_candidates.R \\",
    "    --cluster_run_dir /path/to/cluster_output \\",
    "    --motif_run_dir /path/to/motif_output \\",
    "    --tss_anno_tsv /path/to/region_TSS_association.detail.transcriptTSS.tsv \\",
    "    --expr_tsv /path/to/expression_matrix.tsv \\",
    "    --out_dir /path/to/repressor_score_out",
    "",
    "Required:",
    "  --motif_run_dir      Output directory from run_array_motif_enrichment.py",
    "  --tss_anno_tsv       TSS annotation TSV that contains sample + region_id",
    "  --expr_tsv           Expression matrix TSV (gene x replicate/sample columns)",
    "  --out_dir            Output directory",
    "",
    "Input selection:",
    "  --cluster_run_dir    Output directory from cluster_array_subtype.R",
    "  --cluster_tsv        Cluster assignment TSV; overrides --cluster_run_dir when both are given",
    "  --fit_space_preference  auto|within_sample_scaled|global_scaled|raw  [default: auto]",
    "",
    "Optional:",
    "  --promoter_flag_col  [default: overlap_TSSwin_2000]",
    "  --sample_order       Comma-separated sample order",
    "  --restrict_clusters  Comma-separated clusters to keep, or ALL [default: ALL]",
    "  --expr_id_col        Column name of gene ID in expression matrix [default: first column]",
    "  --expr_sample_regex_map  sample1=^sample1_;sample2=^sample2_",
    "  --gene_id_candidates Comma-separated candidate columns in TSS TSV",
    "                       [default: nearest_TSS_gene_name,gene_symbol_mapped,nearest_TSS_id]",
    "  --motif_q_cutoff     [default: 0.05]",
    "  --min_group_n        [default: 5]",
    "  --top_n_plot         [default: 20]",
    "  --homer_known_filename  [default: knownResults.txt]",
    "  --no_plots           Disable PDF/PNG plot output",
    "  --help               Show help",
    sep = "\n"
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

stop2 <- function(...) {
  stop(paste0(...), call. = FALSE)
}

parse_args <- function(args) {
  opt <- list(
    cluster_run_dir = NULL,
    cluster_tsv = NULL,
    motif_run_dir = NULL,
    tss_anno_tsv = NULL,
    expr_tsv = NULL,
    out_dir = NULL,
    fit_space_preference = "auto",
    promoter_flag_col = "overlap_TSSwin_2000",
    sample_order = NULL,
    restrict_clusters = "ALL",
    expr_id_col = NULL,
    expr_sample_regex_map = NULL,
    gene_id_candidates = "nearest_TSS_gene_name,gene_symbol_mapped,nearest_TSS_id",
    motif_q_cutoff = 0.05,
    min_group_n = 5L,
    top_n_plot = 20L,
    homer_known_filename = "knownResults.txt",
    no_plots = FALSE,
    help = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      opt$help <- TRUE
      i <- i + 1L
    } else if (key == "--no_plots") {
      opt$no_plots <- TRUE
      i <- i + 1L
    } else if (startsWith(key, "--")) {
      if (i == length(args)) stop2("Missing value for ", key)
      val <- args[[i + 1L]]
      nm <- sub("^--", "", key)
      if (!nm %in% names(opt)) stop2("Unknown option: ", key)
      opt[[nm]] <- val
      i <- i + 2L
    } else {
      stop2("Unknown positional argument: ", key)
    }
  }

  opt$motif_q_cutoff <- as.numeric(opt$motif_q_cutoff)
  opt$min_group_n <- as.integer(opt$min_group_n)
  opt$top_n_plot <- as.integer(opt$top_n_plot)
  opt
}

split_csv <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (x == "") return(character(0))
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[out != ""]
}

escape_regex <- function(x) {
  stringr::str_replace_all(x, "([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1")
}

natural_key_cluster <- function(x) {
  idx <- suppressWarnings(as.integer(stringr::str_extract(x, "(?<=C)\\d+")))
  idx[is.na(idx)] <- Inf
  ord <- order(idx, x)
  x[ord]
}

parse_semicolon_map <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (x == "") return(list())
  parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  parts <- parts[parts != ""]
  out <- list()
  for (p in parts) {
    if (!grepl("=", p, fixed = TRUE)) next
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    k <- trimws(kv[1])
    v <- trimws(paste(kv[-1], collapse = "="))
    if (k != "" && v != "") out[[k]] <- v
  }
  out
}

parse_homer_tag <- function(tag) {
  sample_name <- sub("\\.(C[0-9]+)\\..*$", "", tag)
  cluster_name <- stringr::str_match(tag, "\\.(C[0-9]+)\\.")[, 2]
  mode <- dplyr::case_when(
    stringr::str_detect(tag, "mode_internalLinker") ~ "internal_linker",
    stringr::str_detect(tag, "mode_edgeOutside") ~ "edge_outside",
    TRUE ~ NA_character_
  )
  rel_m <- stringr::str_match(tag, "rel([+-]?[0-9]+)to([+-]?[0-9]+)")
  interval_start <- suppressWarnings(as.integer(rel_m[, 2]))
  interval_end <- suppressWarnings(as.integer(rel_m[, 3]))
  list(
    sample = ifelse(is.na(cluster_name), NA_character_, sample_name),
    cluster = cluster_name,
    mode = mode,
    interval_dir = tag,
    interval_start = interval_start,
    interval_end = interval_end
  )
}

locate_cluster_inputs <- function(cluster_run_dir, fit_space_preference = "auto") {
  best_k_path <- file.path(cluster_run_dir, "best_k_selection.tsv")
  if (!file.exists(best_k_path)) stop2("Missing best_k_selection.tsv under: ", cluster_run_dir)

  bk <- data.table::fread(best_k_path) %>% as.data.frame()
  if (!"best_k" %in% colnames(bk) || nrow(bk) == 0) stop2("best_k_selection.tsv does not contain best_k")
  best_k <- suppressWarnings(as.integer(bk$best_k[1]))
  if (!is.finite(best_k)) stop2("Invalid best_k in: ", best_k_path)

  fit_space <- if ("fit_space" %in% colnames(bk)) as.character(bk$fit_space[1]) else NA_character_
  cluster_dir <- file.path(cluster_run_dir, "cluster_tables")
  if (!dir.exists(cluster_dir)) stop2("Missing cluster_tables dir under: ", cluster_run_dir)

  cluster_tsv <- NULL
  if (!is.null(fit_space_preference) && fit_space_preference != "auto") {
    cand <- file.path(cluster_dir, paste0("k", best_k, "_cluster_assignment.z1z2.", fit_space_preference, ".tsv"))
    if (file.exists(cand)) cluster_tsv <- cand
  }
  if (is.null(cluster_tsv) && !is.na(fit_space) && nzchar(fit_space)) {
    cand <- file.path(cluster_dir, paste0("k", best_k, "_cluster_assignment.z1z2.", fit_space, ".tsv"))
    if (file.exists(cand)) cluster_tsv <- cand
  }
  if (is.null(cluster_tsv)) {
    patt <- file.path(cluster_dir, paste0("k", best_k, "_cluster_assignment.z1z2.*.tsv"))
    hits <- Sys.glob(patt)
    if (length(hits) == 0) stop2("Cannot locate best-k cluster assignment TSV under: ", cluster_dir)
    cluster_tsv <- hits[1]
  }

  cluster_feature_tsv <- file.path(cluster_dir, paste0("k", best_k, "_cluster_feature_summary.tsv"))
  if (!file.exists(cluster_feature_tsv)) cluster_feature_tsv <- NA_character_

  list(
    cluster_tsv = cluster_tsv,
    best_k_path = best_k_path,
    best_k = best_k,
    fit_space = fit_space,
    cluster_feature_tsv = cluster_feature_tsv
  )
}

pick_first_existing <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normalize_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- toupper(trimws(as.character(x)))
  x %in% c("TRUE", "T", "1", "YES", "Y")
}

safe_neglog10 <- function(x, floor = .Machine$double.xmin) {
  -log10(pmax(x, floor))
}

safe_wilcox_less <- function(x, y, min_n = 5L) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < min_n || length(y) < min_n) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x = x, y = y, alternative = "less", exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_kw <- function(df) {
  df <- df %>% dplyr::filter(!is.na(group), is.finite(expr_log))
  if (dplyr::n_distinct(df$group) < 2) return(NA_real_)
  tryCatch(
    stats::kruskal.test(expr_log ~ group, data = df)$p.value,
    error = function(e) NA_real_
  )
}

extract_tf_name <- function(x) {
  first_part <- stringr::str_split_fixed(as.character(x), "/", 2)[, 1]
  first_part <- stringr::str_remove(first_part, "\\(.*$")
  first_part <- stringr::str_trim(first_part)
  toupper(first_part)
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

choose_gene_id_col <- function(tss_df, expr_gene_ids, candidates) {
  candidates <- intersect(candidates, colnames(tss_df))
  if (length(candidates) == 0) stop2("No candidate gene ID columns found in TSS annotation table")
  overlap_counts <- vapply(
    candidates,
    function(cc) {
      vals <- as.character(tss_df[[cc]])
      vals <- vals[!is.na(vals) & vals != "" & vals != "FALSE"]
      sum(unique(vals) %in% expr_gene_ids)
    },
    FUN.VALUE = numeric(1)
  )
  candidates[which.max(overlap_counts)]
}

infer_sample_from_expr_col <- function(colname, sample_levels, regex_map = list()) {
  if (length(regex_map) > 0) {
    for (nm in names(regex_map)) {
      rgx <- regex_map[[nm]]
      if (stringr::str_detect(colname, rgx)) return(nm)
    }
  }

  sample_levels2 <- sample_levels[order(nchar(sample_levels), decreasing = TRUE)]
  for (s in sample_levels2) {
    pat <- paste0("^", escape_regex(s), "(?:_|$)")
    if (stringr::str_detect(colname, pat)) return(s)
  }

  prefix <- strsplit(colname, "_", fixed = TRUE)[[1]][1]
  if (prefix %in% sample_levels) return(prefix)
  NA_character_
}

read_expression_long <- function(expr_tsv, expr_id_col = NULL, sample_levels, regex_map = list()) {
  expr_dt <- data.table::fread(expr_tsv) %>% as.data.frame()
  if (ncol(expr_dt) < 2) stop2("Expression TSV must contain at least 2 columns")
  gene_col <- expr_id_col %||% colnames(expr_dt)[1]
  if (!gene_col %in% colnames(expr_dt)) stop2("expr_id_col not found in expression TSV: ", gene_col)
  colnames(expr_dt)[colnames(expr_dt) == gene_col] <- "gene_id"

  expr_long <- expr_dt %>%
    tidyr::pivot_longer(cols = -gene_id, names_to = "replicate", values_to = "expr_raw") %>%
    dplyr::mutate(
      sample = vapply(replicate, infer_sample_from_expr_col, character(1), sample_levels = sample_levels, regex_map = regex_map),
      expr_raw = suppressWarnings(as.numeric(expr_raw))
    ) %>%
    dplyr::filter(!is.na(sample))

  if (nrow(expr_long) == 0) stop2("No expression columns could be mapped to sample names")

  expr_by_sample <- expr_long %>%
    dplyr::group_by(gene_id, sample) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      expr_mean = mean(expr_raw, na.rm = TRUE),
      expr_log = log2(expr_mean + 1),
      .groups = "drop"
    )

  list(expr_long = expr_long, expr_by_sample = expr_by_sample)
}

parse_homer_known_results <- function(f, known_filename = "knownResults.txt") {
  tag <- basename(dirname(f))
  info <- parse_homer_tag(tag)

  if (is.na(info$sample) || is.na(info$cluster)) {
    sample_cluster_dir <- basename(dirname(f))
    interval_dir <- basename(dirname(dirname(dirname(f))))
    sample_name <- sub("\\.(C[0-9]+)\\..*$", "", sample_cluster_dir)
    cluster_name <- stringr::str_match(sample_cluster_dir, "\\.(C[0-9]+)\\.")[, 2]
    rel_m <- stringr::str_match(interval_dir, "rel([+-]?[0-9]+)to([+-]?[0-9]+)")
    info <- list(
      sample = ifelse(is.na(cluster_name), NA_character_, sample_name),
      cluster = cluster_name,
      mode = NA_character_,
      interval_dir = interval_dir,
      interval_start = suppressWarnings(as.integer(rel_m[, 2])),
      interval_end = suppressWarnings(as.integer(rel_m[, 3]))
    )
  }

  dt <- tryCatch(
    data.table::fread(f, sep = "\t", header = TRUE, fill = TRUE, quote = "") %>% as.data.frame(),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  colnames(dt) <- make.names(colnames(dt), unique = TRUE)
  motif_col <- pick_first_existing(colnames(dt), c("Motif.Name", "Motif.Name.", "MotifName"))
  q_col <- grep("^q\\.value", colnames(dt), ignore.case = TRUE, value = TRUE)[1]
  p_col <- grep("^P\\.value", colnames(dt), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(motif_col) || is.na(q_col)) return(NULL)

  data.frame(
    sample = info$sample,
    cluster = info$cluster,
    mode = info$mode,
    tag = tag,
    interval_dir = info$interval_dir,
    interval_start = info$interval_start,
    interval_end = info$interval_end,
    motif_name = as.character(dt[[motif_col]]),
    tf = extract_tf_name(as.character(dt[[motif_col]])),
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else NA_real_,
    q_value = suppressWarnings(as.numeric(dt[[q_col]])),
    stringsAsFactors = FALSE
  )
}

plot_signed_score_one_sample <- function(sample_name, prediction_df, cluster_levels, plot_dir, top_n) {
  df_plot <- prediction_df %>%
    dplyr::filter(sample == sample_name, prediction_score > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = prediction_score, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()
  if (nrow(df_plot) == 0) return(invisible(NULL))

  ncol_fac <- min(3, max(1, length(unique(df_plot$cluster))))
  df_plot <- df_plot %>%
    dplyr::mutate(
      tf_plot = reorder_within(tf, prediction_score_signed, cluster),
      cluster = factor(cluster, levels = cluster_levels)
    )

  p <- ggplot(df_plot, aes(x = prediction_score_signed, y = tf_plot, fill = S_score)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ cluster, scales = "free", ncol = ncol_fac) +
    scale_y_reordered() +
    scale_fill_gradient2(low = "#3B82F6", mid = "grey90", high = "#D62728", midpoint = 0, name = "S[s,c,t]") +
    labs(
      title = paste0(sample_name, ": signed repressor prediction score"),
      subtitle = "Bar length = P[s,c,t] × sign(S[s,c,t]); fill = S[s,c,t]",
      x = "Signed repressor prediction score",
      y = "TF"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(file.path(plot_dir, paste0(sample_name, ".signed_repressor_prediction_score.top", top_n, ".pdf")), p, width = max(10, 4 * ncol_fac), height = 4)
  ggsave(file.path(plot_dir, paste0(sample_name, ".signed_repressor_prediction_score.top", top_n, ".png")), p, width = max(10, 4 * ncol_fac), height = 4, dpi = 300)
}

plot_priority_score_one_sample <- function(sample_name, priority_prediction_df, cluster_levels, plot_dir, top_n) {
  df_plot <- priority_prediction_df %>%
    dplyr::filter(sample == sample_name) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = priority_score, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()
  if (nrow(df_plot) == 0) return(invisible(NULL))

  ncol_fac <- min(3, max(1, length(unique(df_plot$cluster))))
  df_plot <- df_plot %>%
    dplyr::mutate(
      tf_plot = reorder_within(tf, priority_score, cluster),
      cluster = factor(cluster, levels = cluster_levels)
    )

  p <- ggplot(df_plot, aes(x = priority_score, y = tf_plot, fill = S_score)) +
    geom_col(width = 0.72) +
    facet_wrap(~ cluster, scales = "free", ncol = ncol_fac) +
    scale_y_reordered() +
    scale_fill_gradient2(low = "#3B82F6", mid = "grey90", high = "#D62728", midpoint = 0, name = "S[s,c,t]") +
    labs(
      title = paste0(sample_name, ": priority candidate repressors"),
      subtitle = "Only TFs with S[s,c,t] > 0 are retained",
      x = "Priority repressor prediction score",
      y = "TF"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(file.path(plot_dir, paste0(sample_name, ".priority_repressor_prediction_score.top", top_n, ".pdf")), p, width = max(10, 4 * ncol_fac), height = 4)
  ggsave(file.path(plot_dir, paste0(sample_name, ".priority_repressor_prediction_score.top", top_n, ".png")), p, width = max(10, 4 * ncol_fac), height = 4, dpi = 300)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(opt$help)) {
    cat(usage_text(), "\n")
    quit(save = "no", status = 0)
  }

  if (is.null(opt$cluster_tsv) && is.null(opt$cluster_run_dir)) stop2("One of --cluster_tsv or --cluster_run_dir must be provided")
  if (is.null(opt$motif_run_dir)) stop2("--motif_run_dir is required")
  if (is.null(opt$tss_anno_tsv)) stop2("--tss_anno_tsv is required")
  if (is.null(opt$expr_tsv)) stop2("--expr_tsv is required")
  if (is.null(opt$out_dir)) stop2("--out_dir is required")

  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(opt$out_dir, "tables")
  plot_dir <- file.path(opt$out_dir, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(opt$no_plots)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_info <- list(best_k_path = NA_character_, best_k = NA_integer_, fit_space = NA_character_, cluster_feature_tsv = NA_character_)
  if (!is.null(opt$cluster_tsv)) {
    cluster_tsv <- opt$cluster_tsv
  } else {
    cluster_info <- locate_cluster_inputs(opt$cluster_run_dir, opt$fit_space_preference)
    cluster_tsv <- cluster_info$cluster_tsv
  }

  if (!file.exists(cluster_tsv)) stop2("Missing cluster TSV: ", cluster_tsv)
  if (!file.exists(opt$tss_anno_tsv)) stop2("Missing TSS annotation TSV: ", opt$tss_anno_tsv)
  if (!file.exists(opt$expr_tsv)) stop2("Missing expression TSV: ", opt$expr_tsv)
  if (!dir.exists(opt$motif_run_dir)) stop2("Missing motif_run_dir: ", opt$motif_run_dir)

  cluster_dt <- data.table::fread(cluster_tsv) %>% as.data.frame()
  cluster_col <- pick_first_existing(colnames(cluster_dt), c("cluster", grep("_cluster$", colnames(cluster_dt), value = TRUE)))
  if (is.na(cluster_col)) stop2("No cluster column found in cluster TSV")

  need_cols <- c("sample", "region_id", "chr", "start", "end")
  miss <- setdiff(need_cols, colnames(cluster_dt))
  if (length(miss) > 0) stop2("Missing columns in cluster TSV: ", paste(miss, collapse = ", "))

  cluster_assign <- cluster_dt %>%
    dplyr::transmute(
      sample = as.character(sample),
      region_id = as.character(region_id),
      chr = as.character(chr),
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end)),
      width_bp = if ("width_bp" %in% colnames(cluster_dt)) suppressWarnings(as.numeric(width_bp)) else suppressWarnings(as.numeric(end) - as.numeric(start)),
      cluster = as.character(.data[[cluster_col]])
    )

  sample_levels <- split_csv(opt$sample_order)
  if (length(sample_levels) == 0) sample_levels <- unique(cluster_assign$sample)
  sample_levels <- sample_levels[sample_levels %in% unique(cluster_assign$sample)]
  if (length(sample_levels) == 0) sample_levels <- unique(cluster_assign$sample)

  cluster_levels <- sort(unique(na.omit(cluster_assign$cluster)))
  cluster_levels <- natural_key_cluster(cluster_levels)
  restrict_clusters <- split_csv(opt$restrict_clusters)
  if (!(length(restrict_clusters) == 1 && toupper(restrict_clusters) == "ALL") && length(restrict_clusters) > 0) {
    cluster_levels <- cluster_levels[cluster_levels %in% restrict_clusters]
  }
  if (length(cluster_levels) == 0) stop2("No valid cluster levels remain after filtering")
  cluster_assign <- cluster_assign %>% dplyr::filter(cluster %in% cluster_levels, sample %in% sample_levels)

  regex_map <- parse_semicolon_map(opt$expr_sample_regex_map)
  expr_res <- read_expression_long(opt$expr_tsv, opt$expr_id_col, sample_levels, regex_map)
  expr_long <- expr_res$expr_long
  expr_by_sample <- expr_res$expr_by_sample %>%
    dplyr::mutate(sample = factor(sample, levels = sample_levels))
  all_expr_gene_ids <- unique(expr_by_sample$gene_id)

  tss_dt <- data.table::fread(
    opt$tss_anno_tsv,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    data.table = FALSE,
    check.names = FALSE
  )

  if (!opt$promoter_flag_col %in% colnames(tss_dt)) {
    stop2("promoter flag column not found in TSS annotation TSV: ", opt$promoter_flag_col)
  }
  if (!"sample" %in% colnames(tss_dt)) {
    stop2("TSS annotation TSV is missing required column: sample")
  }
  if (!"region_id" %in% colnames(tss_dt)) {
    stop2("TSS annotation TSV is missing required column: region_id")
  }

  tss_dt <- tss_dt %>%
    dplyr::mutate(
      sample = as.character(sample),
      region_id = as.character(region_id)
    )
  gene_id_candidates <- split_csv(opt$gene_id_candidates)
  gene_id_col <- choose_gene_id_col(tss_dt, all_expr_gene_ids, gene_id_candidates)

  tss_use <- tss_dt %>%
    dplyr::left_join(cluster_assign %>% dplyr::select(sample, region_id, cluster), by = c("sample", "region_id"), suffix = c("", ".from_cluster")) %>%
    dplyr::mutate(
      cluster = dplyr::coalesce(as.character(cluster.from_cluster), as.character(cluster)),
      gene_id = as.character(.data[[gene_id_col]]),
      promoter_hit = normalize_bool(.data[[opt$promoter_flag_col]]),
      sample = factor(as.character(sample), levels = sample_levels)
    )

  promoter_hits <- tss_use %>%
    dplyr::filter(
      sample %in% sample_levels,
      promoter_hit,
      !is.na(gene_id),
      gene_id != "",
      gene_id != "FALSE",
      cluster %in% cluster_levels
    ) %>%
    dplyr::distinct(sample, gene_id, cluster, region_id)

  build_gene_groups_one_sample <- function(sample_name) {
    expr_s <- expr_by_sample %>%
      dplyr::filter(sample == sample_name) %>%
      dplyr::mutate(sample = as.character(sample))

    hits_s <- promoter_hits %>%
      dplyr::filter(sample == sample_name, gene_id %in% expr_s$gene_id) %>%
      dplyr::mutate(value = TRUE) %>%
      dplyr::distinct(gene_id, cluster, value) %>%
      tidyr::pivot_wider(names_from = cluster, values_from = value, values_fill = FALSE)

    if (nrow(hits_s) == 0) hits_s <- data.frame(gene_id = character(0), stringsAsFactors = FALSE)
    for (cc in cluster_levels) {
      if (!cc %in% colnames(hits_s)) hits_s[[cc]] <- logical(nrow(hits_s))
    }

    out <- expr_s %>% dplyr::left_join(hits_s, by = "gene_id")
    for (cc in cluster_levels) out[[cc]] <- ifelse(is.na(out[[cc]]), FALSE, out[[cc]])

    flag_mat <- as.matrix(out[, cluster_levels, drop = FALSE])
    n_present <- rowSums(flag_mat)
    only_cluster <- apply(flag_mat, 1, function(z) if (sum(z) == 1) cluster_levels[which(z)] else NA_character_)

    out$group <- dplyr::case_when(
      n_present == 0 ~ "absent",
      n_present == 1 ~ paste0(only_cluster, "_only"),
      length(cluster_levels) == 2 & n_present == 2 ~ "both",
      n_present > 1 ~ "multi",
      TRUE ~ "other"
    )
    out$n_cluster_present <- n_present
    out
  }

  gene_group_df <- purrr::map_dfr(sample_levels, build_gene_groups_one_sample) %>%
    dplyr::mutate(sample = factor(sample, levels = sample_levels))

  data.table::fwrite(
    gene_group_df %>% dplyr::select(sample, gene_id, expr_mean, expr_log, dplyr::all_of(cluster_levels), n_cluster_present, group),
    file.path(table_dir, "gene_group_assignment.by_sample.tsv"),
    sep = "\t"
  )

  compute_cluster_effect_one_sample <- function(sample_name) {
    df_s <- gene_group_df %>% dplyr::filter(sample == sample_name)
    kw_groups <- c("absent", paste0(cluster_levels, "_only"))
    if (length(cluster_levels) == 2 && any(df_s$group == "both")) {
      kw_groups <- c(kw_groups, "both")
    } else if (any(df_s$group == "multi")) {
      kw_groups <- c(kw_groups, "multi")
    }
    kw_p <- safe_kw(df_s %>% dplyr::filter(group %in% kw_groups))

    res <- purrr::map_dfr(cluster_levels, function(cc) {
      group_only <- paste0(cc, "_only")
      x <- df_s %>% dplyr::filter(group == group_only) %>% dplyr::pull(expr_log)
      y <- df_s %>% dplyr::filter(group == "absent") %>% dplyr::pull(expr_log)
      delta <- if (length(x) > 0 && length(y) > 0) stats::median(x, na.rm = TRUE) - stats::median(y, na.rm = TRUE) else NA_real_
      p_expr <- safe_wilcox_less(x, y, min_n = opt$min_group_n)
      data.frame(
        sample = sample_name,
        cluster = cc,
        n_only = length(x[is.finite(x)]),
        n_absent = length(y[is.finite(y)]),
        delta_regulatory = delta,
        p_expr = p_expr,
        kw_p = kw_p,
        stringsAsFactors = FALSE
      )
    })

    res$q_expr <- stats::p.adjust(res$p_expr, method = "BH")
    res$R_score <- ifelse(!is.na(res$delta_regulatory) & res$delta_regulatory < 0 & !is.na(res$q_expr), abs(res$delta_regulatory) * safe_neglog10(res$q_expr), 0)
    res
  }

  cluster_effect_df <- purrr::map_dfr(sample_levels, compute_cluster_effect_one_sample) %>%
    dplyr::mutate(sample = factor(sample, levels = sample_levels), cluster = factor(cluster, levels = cluster_levels))

  data.table::fwrite(cluster_effect_df, file.path(table_dir, "cluster_regulatory_effect.by_sample.tsv"), sep = "\t")

  homer_root <- if (dir.exists(file.path(opt$motif_run_dir, "homer_out"))) file.path(opt$motif_run_dir, "homer_out") else opt$motif_run_dir
  known_result_files <- list.files(homer_root, pattern = paste0(escape_regex(opt$homer_known_filename), "$"), recursive = TRUE, full.names = TRUE)
  if (length(known_result_files) == 0) stop2("No ", opt$homer_known_filename, " files found under: ", homer_root)

  homer_raw <- purrr::map_dfr(known_result_files, parse_homer_known_results, known_filename = opt$homer_known_filename) %>%
    dplyr::filter(
      sample %in% sample_levels,
      cluster %in% cluster_levels,
      !is.na(tf),
      tf != "",
      is.finite(q_value)
    ) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_levels),
      cluster = factor(cluster, levels = cluster_levels),
      sig_bin = q_value < opt$motif_q_cutoff,
      bin_score = ifelse(sig_bin, safe_neglog10(q_value), 0)
    )

  if (nrow(homer_raw) == 0) stop2("No valid HOMER motif records remain after sample/cluster filtering")
  data.table::fwrite(homer_raw, file.path(table_dir, "homer_knownResults.parsed.tsv"), sep = "\t")

  motif_support_df <- homer_raw %>%
    dplyr::group_by(sample, cluster, tf) %>%
    dplyr::summarise(
      M_score = sum(bin_score, na.rm = TRUE),
      n_sig_bins = sum(sig_bin, na.rm = TRUE),
      min_q = ifelse(any(sig_bin), min(q_value[sig_bin], na.rm = TRUE), NA_real_),
      strongest_interval = ifelse(any(sig_bin), interval_dir[which.min(q_value)], NA_character_),
      .groups = "drop"
    )

  tf_per_sample <- homer_raw %>% dplyr::distinct(sample, tf)
  motif_complete_df <- tf_per_sample %>%
    tidyr::crossing(cluster = factor(cluster_levels, levels = cluster_levels)) %>%
    dplyr::left_join(motif_support_df, by = c("sample", "cluster", "tf")) %>%
    dplyr::mutate(
      M_score = dplyr::coalesce(M_score, 0),
      n_sig_bins = dplyr::coalesce(n_sig_bins, 0L),
      min_q = dplyr::coalesce(min_q, NA_real_),
      strongest_interval = dplyr::coalesce(strongest_interval, NA_character_)
    )

  motif_scored_df <- motif_complete_df %>%
    dplyr::group_by(sample, tf) %>%
    dplyr::arrange(cluster, .by_group = TRUE) %>%
    dplyr::mutate(
      max_other_M = purrr::map_dbl(seq_along(M_score), function(i) {
        other <- M_score[-i]
        if (length(other) == 0) return(0)
        suppressWarnings(max(other, na.rm = TRUE))
      }),
      best_other_cluster = purrr::map_chr(seq_along(M_score), function(i) {
        other_idx <- setdiff(seq_along(M_score), i)
        if (length(other_idx) == 0) return(NA_character_)
        other_scores <- M_score[other_idx]
        other_clusters <- as.character(cluster[other_idx])
        other_clusters[which.max(other_scores)]
      }),
      S_score = M_score - max_other_M,
      S_direction = dplyr::case_when(
        S_score > 0 ~ paste0("favor_", cluster),
        S_score < 0 & !is.na(best_other_cluster) ~ paste0("favor_", best_other_cluster),
        TRUE ~ "neutral"
      )
    ) %>%
    dplyr::ungroup()

  prediction_df <- motif_scored_df %>%
    dplyr::left_join(
      cluster_effect_df %>% dplyr::select(sample, cluster, delta_regulatory, p_expr, q_expr, R_score, n_only, n_absent, kw_p),
      by = c("sample", "cluster")
    ) %>%
    dplyr::mutate(
      prediction_score = R_score * M_score,
      prediction_score_signed = prediction_score * sign(S_score),
      priority_score = ifelse(S_score > 0, prediction_score, 0)
    ) %>%
    dplyr::arrange(sample, cluster, dplyr::desc(prediction_score))

  priority_prediction_df <- prediction_df %>%
    dplyr::filter(priority_score > 0) %>%
    dplyr::group_by(sample, cluster) %>%
    dplyr::arrange(dplyr::desc(priority_score), .by_group = TRUE) %>%
    dplyr::mutate(priority_rank = dplyr::row_number()) %>%
    dplyr::ungroup()

  data.table::fwrite(motif_support_df, file.path(table_dir, "motif_support_score.by_sample_cluster_tf.tsv"), sep = "\t")
  data.table::fwrite(prediction_df, file.path(table_dir, "TF_repressor_prediction_score.by_sample_cluster.tsv"), sep = "\t")
  data.table::fwrite(priority_prediction_df, file.path(table_dir, "TF_repressor_prediction_priority.by_sample_cluster.tsv"), sep = "\t")

  top_summary_df <- priority_prediction_df %>%
    dplyr::group_by(sample, cluster) %>%
    dplyr::slice_max(order_by = priority_score, n = opt$top_n_plot, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      sample, cluster, priority_rank, tf,
      delta_regulatory, p_expr, q_expr, R_score,
      M_score, S_score, S_direction,
      prediction_score, prediction_score_signed, priority_score,
      n_sig_bins, strongest_interval
    )
  data.table::fwrite(top_summary_df, file.path(table_dir, paste0("top", opt$top_n_plot, "_priority_repressors.by_sample_cluster.tsv")), sep = "\t")

  if (!isTRUE(opt$no_plots)) {
    invisible(purrr::walk(sample_levels, plot_signed_score_one_sample, prediction_df = prediction_df, cluster_levels = cluster_levels, plot_dir = plot_dir, top_n = opt$top_n_plot))
    invisible(purrr::walk(sample_levels, plot_priority_score_one_sample, priority_prediction_df = priority_prediction_df, cluster_levels = cluster_levels, plot_dir = plot_dir, top_n = opt$top_n_plot))
  }

  metadata_df <- data.frame(
    cluster_tsv = cluster_tsv,
    cluster_run_dir = opt$cluster_run_dir %||% NA_character_,
    best_k_selection_tsv = cluster_info$best_k_path %||% NA_character_,
    best_k = cluster_info$best_k %||% NA_integer_,
    fit_space = cluster_info$fit_space %||% NA_character_,
    cluster_feature_tsv = cluster_info$cluster_feature_tsv %||% NA_character_,
    cluster_column_used = cluster_col,
    motif_run_dir = opt$motif_run_dir,
    homer_root = homer_root,
    tss_anno_tsv = opt$tss_anno_tsv,
    expr_tsv = opt$expr_tsv,
    expr_id_col = opt$expr_id_col %||% NA_character_,
    gene_id_col_used = gene_id_col,
    promoter_flag_col = opt$promoter_flag_col,
    sample_order = paste(sample_levels, collapse = ","),
    cluster_levels = paste(cluster_levels, collapse = ","),
    motif_q_cutoff = opt$motif_q_cutoff,
    min_group_n = opt$min_group_n,
    top_n_plot = opt$top_n_plot,
    n_expr_rows = nrow(expr_long),
    n_expr_gene_sample = nrow(expr_by_sample),
    n_cluster_regions = nrow(cluster_assign),
    n_promoter_hits = nrow(promoter_hits),
    n_homer_rows = nrow(homer_raw),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(metadata_df, file.path(table_dir, "run_metadata.tsv"), sep = "\t")

  message("Done.")
  message("Tables: ", table_dir)
  if (!isTRUE(opt$no_plots)) message("Plots:  ", plot_dir)
}

main()