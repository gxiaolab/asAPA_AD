#!/usr/bin/env Rscript

# Usage:
# Rscript aqtl_overlap_analysis_tensorQTL_only.R \
#   FP.top.cis_qtl.txt IFG.top.cis_qtl.txt PG.top.cis_qtl.txt STG.top.cis_qtl.txt \
#   asapa_snps.txt control_snps.txt asapa_genes.txt control_genes.txt \
#   outdir asapa_snp_ratio_file.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(
    paste(
      "Expected 10 arguments:\n",
      "1 FP, 2 IFG, 3 PG, 4 STG,\n",
      "5 asAPA SNPs, 6 control SNPs, 7 asAPA genes,",
      "8 control genes, 9 output directory, 10 asAPA SNP ratio file"
    )
  )
}

FP                <- args[1]
IFG               <- args[2]
PG                <- args[3]
STG               <- args[4]
asapa_file        <- args[5]
control_file      <- args[6]
asapa_gene_file   <- args[7]
control_gene_file <- args[8]
outdir            <- args[9]
asapa_ratio_file  <- args[10]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

tensor_dir <- file.path(outdir, "tensorQTL_combined")
dir.create(tensor_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- brewer.pal(3, "Set1")[1:2]
names(group_colors) <- c("asAPA", "control")

# General helpers
get_column <- function(dt, preferred_names = character(), fallback_index = NULL,
                       label = "column") {
  hit <- preferred_names[preferred_names %in% names(dt)]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (!is.null(fallback_index) && ncol(dt) >= fallback_index) {
    message(
      "Using column ", fallback_index, " ('", names(dt)[fallback_index],
      "') as ", label, "."
    )
    return(names(dt)[fallback_index])
  }

  stop(
    "Could not identify ", label, ". Available columns are:\n",
    paste(names(dt), collapse = ", ")
  )
}

require_column_index <- function(dt, index, label) {
  if (ncol(dt) < index) {
    stop(
      label, " is expected in column ", index,
      ", but the input has only ", ncol(dt), " columns."
    )
  }

  col <- names(dt)[index]
  message("Using column ", index, " ('", col, "') as ", label, ".")
  col
}

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

parse_variant_id <- function(x) {
  # Supports, for example:
  #   1_680742_A_T
  #   chr1_680742_A_T
  #   chr1:680742
  x <- as.character(x)

  out_chr <- rep(NA_character_, length(x))
  out_pos <- rep(NA_real_, length(x))
  out_ref <- rep(NA_character_, length(x))
  out_alt <- rep(NA_character_, length(x))

  colon <- grepl("^chr[^:]+:[0-9]+$", x)
  out_chr[colon] <- sub(":.*$", "", x[colon])
  out_pos[colon] <- suppressWarnings(as.numeric(sub("^.*:", "", x[colon])))

  underscore <- !colon & grepl("^[^_]+_[0-9]+_[^_]+_[^_]+$", x)
  pieces <- strsplit(x[underscore], "_", fixed = TRUE)
  out_chr[underscore] <- vapply(
    pieces, function(z) normalize_chr(z[1]), character(1)
  )
  out_pos[underscore] <- suppressWarnings(vapply(
    pieces, function(z) as.numeric(z[2]), numeric(1)
  ))
  out_ref[underscore] <- toupper(vapply(pieces, function(z) z[3], character(1)))
  out_alt[underscore] <- toupper(vapply(pieces, function(z) z[4], character(1)))

  data.table(
    chr = out_chr,
    pos = out_pos,
    ref = out_ref,
    alt = out_alt,
    SNP = paste0(
      out_chr, ":",
      format(out_pos, scientific = FALSE, trim = TRUE)
    )
  )
}

extract_gene <- function(phenotype) {
  # ENSG..._UTRk|chr:start-end -> ENSG...
  sub("_UTR.*$", "", sub("\\|.*$", "", as.character(phenotype)))
}

parse_phenotype_region <- function(phenotype) {
  # Extracts chr:start-end from column 1 / phenotype_id.
  phenotype <- as.character(phenotype)
  matched <- regmatches(
    phenotype,
    regexpr("chr[^:|_]+:[0-9]+-[0-9]+", phenotype, perl = TRUE)
  )
  matched[matched == ""] <- NA_character_

  chr <- sub(":.*$", "", matched)
  coords <- sub("^.*:", "", matched)
  start <- suppressWarnings(as.numeric(sub("-.*$", "", coords)))
  end   <- suppressWarnings(as.numeric(sub("^.*-", "", coords)))

  data.table(region = matched, chr = chr, start = start, end = end)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_max_abs <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(abs(x))
}

safe_chisq <- function(tab, output_file) {
  if (!all(dim(tab) == c(2, 2))) {
    writeLines("A complete 2 x 2 table could not be formed.", output_file)
    return(NA_real_)
  }

  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    writeLines(
      "Chi-squared test could not be performed because a row or column total is zero.",
      output_file
    )
    return(NA_real_)
  }

  result <- tryCatch(
    chisq.test(tab),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    writeLines(
      paste("Chi-squared test failed:", conditionMessage(result)),
      output_file
    )
    return(NA_real_)
  }

  capture.output(result, file = output_file)
  result$p.value
}

format_p <- function(p) {
  if (is.na(p)) "NA" else formatC(p, format = "e", digits = 2)
}

save_table <- function(dt, filename) {
  fwrite(dt, filename, sep = "\t", quote = FALSE, na = "NA")
}

make_enrichment_plot <- function(summary_dt, title, y_label, chisq_p, filename) {
  summary_dt[, group := factor(group, levels = c("asAPA", "control"))]

  ymax <- max(summary_dt$proportion, na.rm = TRUE)
  if (!is.finite(ymax) || ymax == 0) ymax <- 0.05

  p <- ggplot(summary_dt, aes(group, proportion, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = paste0(n_positive, "/", n_total)),
      vjust = -0.4,
      size = 3.5
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, ymax * 1.28),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      subtitle = paste0("Chi-squared p = ", format_p(chisq_p)),
      x = NULL,
      y = y_label
    ) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(filename, p, width = 4.4, height = 4.5)
}

qq_frame <- function(pvals, group) {
  pvals <- suppressWarnings(as.numeric(pvals))
  pvals <- sort(pvals[is.finite(pvals) & pvals > 0 & pvals <= 1])

  if (length(pvals) == 0) {
    return(data.table())
  }

  data.table(
    expected = -log10(ppoints(length(pvals))),
    observed = -log10(pvals),
    group = group
  )
}

make_qq_plot <- function(asapa_p, control_p, title, filename,
                         perform_ks = TRUE, subtitle_note = NULL) {
  a <- qq_frame(asapa_p, "asAPA")
  c <- qq_frame(control_p, "control")
  qq <- rbindlist(list(a, c), use.names = TRUE, fill = TRUE)

  if (nrow(qq) == 0) {
    writeLines(
      "No valid nominal p-values were available for this QQ plot.",
      sub("\\.pdf$", ".txt", filename)
    )
    return(invisible(NULL))
  }

  a_clean <- suppressWarnings(as.numeric(asapa_p))
  a_clean <- a_clean[is.finite(a_clean) & a_clean > 0 & a_clean <= 1]
  c_clean <- suppressWarnings(as.numeric(control_p))
  c_clean <- c_clean[is.finite(c_clean) & c_clean > 0 & c_clean <= 1]

  ks_p <- if (perform_ks && length(a_clean) > 1 && length(c_clean) > 1) {
    suppressWarnings(ks.test(a_clean, c_clean)$p.value)
  } else {
    NA_real_
  }

  subtitle_text <- if (perform_ks) {
    paste0(
      "KS p = ", format_p(ks_p),
      "; n(asAPA) = ", length(a_clean),
      "; n(control) = ", length(c_clean)
    )
  } else {
    paste0(
      "n(asAPA) = ", length(a_clean),
      "; n(control) = ", length(c_clean)
    )
  }

  if (!is.null(subtitle_note) && nzchar(subtitle_note)) {
    subtitle_text <- paste0(subtitle_text, "
", subtitle_note)
  }

  max_axis <- max(c(qq$expected, qq$observed), na.rm = TRUE)

  p <- ggplot(qq, aes(expected, observed, color = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 0.65, size = 1.4) +
    #coord_equal(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
    scale_color_manual(values = group_colors) +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = expression(Expected~~-log[10](p)),
      y = expression(Observed~~-log[10](p)),
      color = NULL
    ) +
    theme_classic()

  ggsave(filename, p, width = 5.3, height = 5)
  invisible(qq)
}

load_candidate_snps <- function(file, group) {
  dt <- fread(file, header = FALSE)

  if (ncol(dt) < 2) {
    stop(group, " SNP file must contain chromosome in column 1 and position in column 2.")
  }

  dt <- dt[, .(
    chr = normalize_chr(V1),
    pos = suppressWarnings(as.numeric(V2))
  )]
  dt <- dt[!is.na(chr) & is.finite(pos)]
  dt[, SNP := paste0(chr, ":", format(pos, scientific = FALSE, trim = TRUE))]
  unique(dt[, .(chr, pos, SNP)])
}

load_candidate_genes <- function(file) {
  dt <- fread(file, header = FALSE)
  unique(as.character(dt[[1]]))
}

prepare_candidate_sets <- function(asapa_snps, control_snps,
                                   asapa_genes, control_genes) {
  shared_snps <- intersect(asapa_snps$SNP, control_snps$SNP)
  shared_genes <- intersect(asapa_genes, control_genes)

  if (length(shared_snps) > 0) {
    warning(
      length(shared_snps),
      " SNPs occur in both candidate files. They are excluded from both groups."
    )
    asapa_snps <- asapa_snps[!SNP %in% shared_snps]
    control_snps <- control_snps[!SNP %in% shared_snps]
  }

  if (length(shared_genes) > 0) {
    warning(
      length(shared_genes),
      " genes occur in both candidate files. They are excluded from both groups."
    )
    asapa_genes <- setdiff(asapa_genes, shared_genes)
    control_genes <- setdiff(control_genes, shared_genes)
  }

  list(
    asapa_snps = asapa_snps,
    control_snps = control_snps,
    asapa_genes = asapa_genes,
    control_genes = control_genes,
    shared_snps = shared_snps,
    shared_genes = shared_genes
  )
}

# Analysis functions
make_tissue_aqtl_count_plot <- function(qtl, output_file,
                                        output_table = NULL) {
  required_cols <- c("tissue", "variant", "phenotype", "significant")
  missing_cols <- setdiff(required_cols, names(qtl))

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for tissue aQTL plot: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  unique_aqtls <- qtl[
    !is.na(tissue) & !is.na(variant) & !is.na(phenotype),
    .(significant = any(significant == TRUE, na.rm = TRUE)),
    by = .(tissue, variant, phenotype)
  ]

  tissue_counts <- unique_aqtls[, .(
    total_aqtls = .N,
    significant_aqtls = sum(significant, na.rm = TRUE)
  ), by = tissue]

  tissue_counts[, nonsignificant_aqtls :=
                  total_aqtls - significant_aqtls]
  tissue_counts <- tissue_counts[match(tissue, c("FP", "IFG", "PG", "STG"))]
  tissue_counts[, tissue := factor(tissue, levels = c("FP", "IFG", "PG", "STG"))]

  if (!is.null(output_table)) save_table(tissue_counts, output_table)

  plot_long <- melt(
    tissue_counts,
    id.vars = "tissue",
    measure.vars = c("nonsignificant_aqtls", "significant_aqtls"),
    variable.name = "status",
    value.name = "count"
  )
  plot_long[, status := factor(
    status,
    levels = c("nonsignificant_aqtls", "significant_aqtls"),
    labels = c("Not significant", "Significant")
  )]

  p <- ggplot(plot_long, aes(tissue, count, fill = status)) +
    geom_col(width = 0.72) +
    geom_text(
      data = tissue_counts,
      aes(x = tissue, y = total_aqtls, label = total_aqtls),
      inherit.aes = FALSE,
      vjust = -0.35,
      size = 3.8
    ) +
    scale_fill_manual(values = c(
      "Not significant" = "#9ECAE1",
      "Significant" = "#3182BD"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = "Detected aQTLs across brain regions",
      subtitle = "Dark shading indicates q < 0.05",
      x = "Brain region",
      y = "Number of unique variant-UTR aQTLs",
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave(output_file, p, width = 7, height = 5.5)
  invisible(tissue_counts)
}

load_asapa_ratio <- function(file) {
  dt <- fread(
    file,
    header = FALSE,
    col.names = c(
      "chr", "pos", "ref", "alt",
      "ref_count", "alt_count", "total_count", "allelic_ratio"
    )
  )

  dt[, chr := normalize_chr(chr)]
  dt[, pos := suppressWarnings(as.numeric(pos))]
  dt[, ref := toupper(as.character(ref))]
  dt[, alt := toupper(as.character(alt))]
  dt[, ref_count := suppressWarnings(as.numeric(ref_count))]
  dt[, alt_count := suppressWarnings(as.numeric(alt_count))]
  dt[, total_count := suppressWarnings(as.numeric(total_count))]

  dt <- dt[
    !is.na(chr) & is.finite(pos) &
      !is.na(ref) & ref != "" & !is.na(alt) & alt != "" &
      is.finite(ref_count) & is.finite(alt_count) &
      ref_count >= 0 & alt_count >= 0
  ]

  dt[, SNP := paste0(chr, ":", format(pos, scientific = FALSE, trim = TRUE))]

  dt <- dt[, .(
    ref_count = sum(ref_count, na.rm = TRUE),
    alt_count = sum(alt_count, na.rm = TRUE)
  ), by = .(chr, pos, ref, alt, SNP)]

  dt[, total_count := ref_count + alt_count]
  dt <- dt[is.finite(total_count) & total_count > 0]
  dt[, allelic_ratio := ref_count / total_count]

  # tensorQTL slope is ALT-dosage oriented, so orient the allelic imbalance
  # toward ALT as well: positive = ALT favored; negative = REF favored.
  dt[, alt_imbalance := 1 - 2 * allelic_ratio]

  # Exclude ambiguous duplicate allele definitions at the same genomic position.
  duplicate_positions <- dt[, .N, by = SNP][N > 1, SNP]
  if (length(duplicate_positions) > 0) {
    warning(
      length(duplicate_positions),
      " SNP positions have multiple REF/ALT definitions in the ratio file; ",
      "they are excluded from signed effect-size analyses."
    )
    dt <- dt[!SNP %in% duplicate_positions]
  }

  dt
}

make_effect_imbalance_plot <- function(qtl, asapa_ratio,
                                       output_file,
                                       output_table = NULL,
                                       output_correlation = NULL) {
  required_cols <- c(
    "SNP", "ref", "alt", "slope", "qvalue", "nominal_p",
    "significant", "tissue", "phenotype"
  )
  missing_cols <- setdiff(required_cols, names(qtl))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for effect-size plot: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  strongest_qtl <- copy(
  qtl[
    significant == TRUE &
      !is.na(SNP) &
      is.finite(slope) &
      is.finite(qvalue) &
      is.finite(nominal_p)
  ]
)

setorder(
  strongest_qtl,
  SNP,
  qvalue,
  nominal_p
)

# Data are already sorted from strongest to weakest association.
# Keep the first row for each oriented variant without using .SD/by.
strongest_qtl <- strongest_qtl[
  !duplicated(paste(SNP, ref, alt, sep = "|"))
]

  plot_dt <- merge(
    asapa_ratio,
    strongest_qtl,
    by = c("SNP", "ref", "alt"),
    all = FALSE
  )

  if (!"slope" %in% names(plot_dt)) {
    stop(
      "The merged effect-size table does not contain 'slope'. Columns are: ",
      paste(names(plot_dt), collapse = ", ")
    )
  }

  plot_dt <- plot_dt[is.finite(alt_imbalance) & is.finite(slope)]

  if (nrow(plot_dt) < 3) {
    warning("Fewer than three matched SNPs were available for the effect-size plot.")
    return(invisible(NULL))
  }

  cor_test <- cor.test(
    plot_dt$alt_imbalance,
    plot_dt$slope,
    method = "spearman",
    exact = FALSE
  )

  annotation <- paste0(
    "Spearman rho = ", sprintf("%.3f", unname(cor_test$estimate)),
    "\nP = ", formatC(cor_test$p.value, format = "e", digits = 2),
    "\nn = ", nrow(plot_dt)
  )

  if (!is.null(output_table)) save_table(plot_dt, output_table)
  if (!is.null(output_correlation)) capture.output(cor_test, file = output_correlation)

  p <- ggplot(plot_dt, aes(alt_imbalance, slope)) +
    geom_point(alpha = 0.65, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate(
      "text", x = Inf, y = Inf, label = annotation,
      hjust = 1.08, vjust = 1.15, size = 4
    ) +
    labs(
      title = "Signed aQTL effect versus allelic imbalance",
      subtitle = "Strongest significant tensorQTL association retained per asAPA SNP",
      x = "ALT-oriented allelic imbalance (1 - 2 × reference ratio)",
      y = "tensorQTL slope"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(output_file, p, width = 6.5, height = 5.5)
  invisible(cor_test)
}

make_effect_imbalance_plot_by_tissue <- function(
    qtl,
    asapa_ratio,
    output_file,
    output_table = NULL,
    output_correlation = NULL
) {
  required_cols <- c(
    "SNP",
    "ref",
    "alt",
    "slope",
    "qvalue",
    "nominal_p",
    "significant",
    "tissue",
    "phenotype"
  )

  missing_cols <- setdiff(required_cols, names(qtl))

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for tissue-specific effect-size plot: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Retain the strongest significant association for each SNP
  # separately within each tissue.
  strongest_qtl_by_tissue <- copy(
  qtl[
    significant == TRUE &
      !is.na(SNP) &
      !is.na(tissue) &
      is.finite(slope) &
      is.finite(qvalue) &
      is.finite(nominal_p)
  ]
)

setorder(
  strongest_qtl_by_tissue,
  tissue,
  SNP,
  qvalue,
  nominal_p
)

# Data are already sorted from strongest to weakest association.
# Keep the first row for each oriented variant within each tissue.
strongest_qtl_by_tissue <- strongest_qtl_by_tissue[
  !duplicated(paste(tissue, SNP, ref, alt, sep = "|"))
]

  plot_dt <- merge(
    asapa_ratio,
    strongest_qtl_by_tissue,
    by = c("SNP", "ref", "alt"),
    all = FALSE
  )

  if (!"slope" %in% names(plot_dt)) {
    stop(
      "The merged tissue-specific effect-size table does not contain 'slope'. Columns are: ",
      paste(names(plot_dt), collapse = ", ")
    )
  }

  plot_dt <- plot_dt[
    is.finite(alt_imbalance) &
      is.finite(slope)
  ]

  if (nrow(plot_dt) < 3) {
    warning(
      "Fewer than three matched SNP-tissue associations were ",
      "available for the tissue-specific scatterplot."
    )
    return(invisible(NULL))
  }

  # Keep the intended brain-region order.
  plot_dt[
    ,
    tissue := factor(
      tissue,
      levels = c("FP", "IFG", "PG", "STG")
    )
  ]

  # Calculate one Spearman correlation per tissue.
  tissue_correlations <- plot_dt[
    ,
    {
      valid <- is.finite(alt_imbalance) &
        is.finite(slope)

      n_valid <- sum(valid)

      if (n_valid >= 3 &&
          uniqueN(alt_imbalance[valid]) > 1 &&
          uniqueN(slope[valid]) > 1) {

        test <- suppressWarnings(
          cor.test(
            alt_imbalance[valid],
            slope[valid],
            method = "spearman",
            exact = FALSE
          )
        )

        list(
          n = n_valid,
          rho = unname(test$estimate),
          pvalue = test$p.value
        )
      } else {
        list(
          n = n_valid,
          rho = NA_real_,
          pvalue = NA_real_
        )
      }
    },
    by = tissue
  ]

  tissue_correlations[
    ,
    label := paste0(
      "rho = ",
      ifelse(
        is.na(rho),
        "NA",
        sprintf("%.3f", rho)
      ),
      "\nP = ",
      ifelse(
        is.na(pvalue),
        "NA",
        formatC(
          pvalue,
          format = "e",
          digits = 2
        )
      ),
      "\nn = ",
      n
    )
  ]

  # Coordinates for panel annotations.
  annotation_dt <- plot_dt[
    ,
    .(
      x = max(alt_imbalance, na.rm = TRUE),
      y = max(slope, na.rm = TRUE)
    ),
    by = tissue
  ]

  annotation_dt <- merge(
    annotation_dt,
    tissue_correlations[
      ,
      .(tissue, label)
    ],
    by = "tissue",
    all.x = TRUE
  )

  if (!is.null(output_table)) {
    save_table(
      plot_dt,
      output_table
    )
  }

  if (!is.null(output_correlation)) {
    save_table(
      tissue_correlations,
      output_correlation
    )
  }

  p <- ggplot(
    plot_dt,
    aes(
      x = alt_imbalance,
      y = slope
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 2
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    geom_text(
      data = annotation_dt,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.4
    ) +
    facet_wrap(
      ~ tissue,
      ncol = 2
    ) +
    labs(
      title = paste0(
        "Signed aQTL effect versus allelic imbalance ",
        "by brain region"
      ),
      subtitle = paste0(
        "Strongest significant association retained ",
        "per SNP within each tissue"
      ),
      x = "ALT-oriented allelic imbalance (1 - 2 × reference ratio)",
      y = "tensorQTL slope"
    ) +
    theme_classic(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      )
    )

  ggsave(
    output_file,
    p,
    width = 9,
    height = 7
  )

  invisible(
    tissue_correlations
  )
}

make_regional_effect_imbalance_plot <- function(
    qtl,
    asapa_ratio,
    output_file,
    output_table = NULL,
    output_correlation = NULL
) {
  required_cols <- c(
    "phenotype", "tissue", "region_chr", "region_start", "region_end",
    "slope", "qvalue", "nominal_p", "significant"
  )
  missing_cols <- setdiff(required_cols, names(qtl))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for regional effect-size plot: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Keep the strongest significant aQTL association for each UTR in each
  # tissue. The aQTL variant need not be the asAPA SNP in this regional plot.
  regional_qtl <- copy(qtl[
    significant == TRUE &
      !is.na(phenotype) & !is.na(tissue) &
      !is.na(region_chr) & is.finite(region_start) & is.finite(region_end) &
      is.finite(slope) & is.finite(qvalue) & is.finite(nominal_p)
  ])

  setorder(regional_qtl, tissue, phenotype, qvalue, nominal_p)
  regional_qtl <- regional_qtl[
    !duplicated(paste(tissue, phenotype, sep = "|"))
  ]

  if (nrow(regional_qtl) == 0) {
    warning("No significant UTR-level aQTL associations were available.")
    return(invisible(NULL))
  }

  regions <- regional_qtl[, .(
    tissue,
    phenotype,
    region_chr,
    region_start,
    region_end,
    qtl_variant = variant,
    qtl_slope = slope,
    qtl_abs_slope = abs(slope),
    qvalue,
    nominal_p
  )]
  regions[, region_id := .I]

  ratios <- copy(asapa_ratio)
  ratios[, `:=`(start = pos, end = pos)]

  reg <- regions[, .(
    chr = region_chr,
    start = region_start,
    end = region_end,
    region_id,
    tissue,
    phenotype,
    qtl_variant,
    qtl_slope,
    qtl_abs_slope,
    qvalue,
    nominal_p
  )]

  setkey(reg, chr, start, end)
  setkey(ratios, chr, start, end)

  overlaps <- foverlaps(
    ratios,
    reg,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "start", "end"),
    type = "within",
    nomatch = 0L
  )

  if (nrow(overlaps) == 0) {
    warning("No asAPA ratio SNPs overlapped significant aQTL UTRs.")
    return(invisible(NULL))
  }

  # One independent plotting unit per UTR/tissue. Because the aQTL variant and
  # asAPA SNP may differ, compare magnitudes rather than signed allele effects.
  plot_dt <- overlaps[, .(
    median_abs_asapa_imbalance = median(abs(alt_imbalance), na.rm = TRUE),
    mean_abs_asapa_imbalance = mean(abs(alt_imbalance), na.rm = TRUE),
    n_asapa_snps = uniqueN(SNP),
    overlapping_asapa_snps = paste(sort(unique(SNP)), collapse = ","),
    qtl_variant = qtl_variant[1],
    qtl_slope = qtl_slope[1],
    qtl_abs_slope = qtl_abs_slope[1],
    qvalue = qvalue[1],
    nominal_p = nominal_p[1]
  ), by = .(region_id, tissue, phenotype)]

  plot_dt <- plot_dt[
    is.finite(median_abs_asapa_imbalance) & is.finite(qtl_abs_slope)
  ]

  if (nrow(plot_dt) < 3) {
    warning("Fewer than three significant UTRs were available for the regional plot.")
    return(invisible(NULL))
  }

  cor_test <- cor.test(
    plot_dt$median_abs_asapa_imbalance,
    plot_dt$qtl_abs_slope,
    method = "spearman",
    exact = FALSE
  )

  annotation <- paste0(
    "Spearman rho = ", sprintf("%.3f", unname(cor_test$estimate)),
    "\nP = ", formatC(cor_test$p.value, format = "e", digits = 2),
    "\nn = ", nrow(plot_dt), " UTRs"
  )

  if (!is.null(output_table)) save_table(plot_dt, output_table)
  if (!is.null(output_correlation)) {
    capture.output(cor_test, file = output_correlation)
  }

  p <- ggplot(
    plot_dt,
    aes(x = median_abs_asapa_imbalance, y = log2(qtl_abs_slope))
  ) +
    geom_point(aes(size = n_asapa_snps), alpha = 0.65) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate(
      "text", x = Inf, y = Inf, label = annotation,
      hjust = 1.08, vjust = 1.15, size = 4
    ) +
    scale_size_continuous(name = "asAPA SNPs in UTR") +
    labs(
      title = "Regional aQTL effect versus asAPA imbalance",
      subtitle = paste0(
        "One point per significant UTR and tissue; aQTL and asAPA SNPs ",
        "need not be identical"
      ),
      x = "Median absolute asAPA allelic imbalance in UTR",
      y = "Log2 of absolute tensorQTL slope"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(output_file, p, width = 7, height = 5.8)
  invisible(cor_test)
}

run_direct_snp_analysis <- function(qtl, asapa_snps, control_snps,
                                    analysis_dir, prefix) {
  snp_summary <- qtl[, .(
    significant = any(significant, na.rm = TRUE),
    min_p = safe_min(nominal_p)
  ), by = SNP]

  candidates <- unique(rbindlist(list(
    asapa_snps[, .(SNP, group = "asAPA")],
    control_snps[, .(SNP, group = "control")]
  )))

  overlap <- merge(candidates, snp_summary, by = "SNP", all.x = TRUE)
  overlap[, tested := !is.na(min_p)]
  overlap[is.na(significant), significant := FALSE]

  save_table(
    overlap,
    file.path(analysis_dir, paste0(prefix, "_direct_snp_overlap.tsv"))
  )

  # Use every candidate SNP as the denominator.
  tab <- table(
    factor(overlap$group, levels = c("asAPA", "control")),
    factor(overlap$significant, levels = c(TRUE, FALSE))
  )
  colnames(tab) <- c("Significant", "Not_significant")

  write.table(
    tab,
    file.path(analysis_dir, paste0(prefix, "_direct_snp_enrichment_table.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )

  chisq_p <- safe_chisq(
    tab,
    file.path(analysis_dir, paste0(prefix, "_direct_snp_chisq_test.txt"))
  )

  plot_dt <- overlap[, .(
    n_positive = sum(significant, na.rm = TRUE),
    n_total = .N,
    proportion = sum(significant, na.rm = TRUE) / .N,
    n_tested = sum(tested, na.rm = TRUE)
  ), by = group]

  save_table(
    plot_dt,
    file.path(analysis_dir, paste0(prefix, "_direct_snp_enrichment_counts.tsv"))
  )

  make_enrichment_plot(
    plot_dt,
    paste0(prefix, ": significant aQTL SNPs in total candidate set"),
    "Proportion of total candidate SNPs significant",
    chisq_p,
    file.path(analysis_dir, paste0(prefix, "_direct_snp_enrichment.pdf"))
  )

  tested <- overlap[tested == TRUE]

  make_qq_plot(
    tested[group == "asAPA", min_p],
    tested[group == "control", min_p],
    paste0(prefix, ": direct SNP QQ plot"),
    file.path(analysis_dir, paste0(prefix, "_direct_snp_qqplot.pdf"))
  )
}
map_candidates_to_regions <- function(qtl, candidates) {
  regions <- unique(qtl[
    !is.na(region_chr) & is.finite(region_start) & is.finite(region_end),
    .(
      region_chr,
      region_start,
      region_end,
      significant,
      nominal_p,
      phenotype,
      gene,
      tissue
    )
  ])

  if (nrow(regions) == 0) {
    return(data.table())
  }

  reg <- regions[, .(
    chr = region_chr,
    start = region_start,
    end = region_end,
    significant,
    nominal_p,
    phenotype,
    gene,
    tissue
  )]
  reg[, region_row := .I]

  cand <- copy(candidates)
  cand[, `:=`(start = pos, end = pos)]
  cand[, candidate_row := .I]

  setkey(reg, chr, start, end)
  setkey(cand, chr, start, end)

  hits <- foverlaps(
    cand,
    reg,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "start", "end"),
    type = "within",
    nomatch = 0L
  )

  if (nrow(hits) == 0) {
    return(data.table())
  }

  hits[, .(
    SNP,
    group,
    phenotype,
    gene,
    region_chr = chr,
    region_start = start,
    region_end = end,
    snp_position = i.start,
    significant,
    nominal_p,
    tissue
  )]
}

run_region_snp_analysis <- function(qtl, asapa_snps, control_snps,
                                    analysis_dir, prefix) {
  candidates <- rbindlist(list(
    asapa_snps[, .(chr, pos, SNP, group = "asAPA")],
    control_snps[, .(chr, pos, SNP, group = "control")]
  ))

  hits <- map_candidates_to_regions(qtl, candidates)

  if (nrow(hits) > 0) {
  region_summary <- hits[, .(
    in_any_qtl_region = TRUE,

    in_significant_qtl_region = any(
      significant == TRUE,
      na.rm = TRUE
    ),

    # Minimum p-value among every overlapping QTL UTR.
    min_region_p = safe_min(nominal_p),

    # Minimum p-value only among significant overlapping QTL UTRs.
    min_significant_region_p = safe_min(
      nominal_p[significant == TRUE]
    ),

    n_overlapping_qtl_rows = .N,

    n_overlapping_significant_qtl_rows = sum(
      significant == TRUE,
      na.rm = TRUE
    )
  ), by = .(SNP, group)]
} else {
  region_summary <- data.table(
    SNP = character(),
    group = character(),
    in_any_qtl_region = logical(),
    in_significant_qtl_region = logical(),
    min_region_p = numeric(),
    min_significant_region_p = numeric(),
    n_overlapping_qtl_rows = integer(),
    n_overlapping_significant_qtl_rows = integer()
  )
}

  summary <- merge(
    candidates[, .(SNP, group)],
    region_summary,
    by = c("SNP", "group"),
    all.x = TRUE
  )

  summary[is.na(in_any_qtl_region), in_any_qtl_region := FALSE]
  summary[is.na(in_significant_qtl_region), in_significant_qtl_region := FALSE]
  summary[is.na(n_overlapping_qtl_rows), n_overlapping_qtl_rows := 0L]
  summary[
    is.na(n_overlapping_significant_qtl_rows),
    n_overlapping_significant_qtl_rows := 0L
  ]

  save_table(
    summary,
    file.path(analysis_dir, paste0(prefix, "_utr_region_snp_overlap.tsv"))
  )

  save_table(
    hits,
    file.path(analysis_dir, paste0(prefix, "_utr_region_all_pairwise_hits.tsv"))
  )

  tab <- table(
    factor(summary$group, levels = c("asAPA", "control")),
    factor(
      summary$in_significant_qtl_region,
      levels = c(TRUE, FALSE)
    )
  )
  colnames(tab) <- c("In_significant_QTL_UTR", "Not_in_significant_QTL_UTR")

  write.table(
    tab,
    file.path(analysis_dir, paste0(prefix, "_utr_region_enrichment_table.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )

  chisq_p <- safe_chisq(
    tab,
    file.path(analysis_dir, paste0(prefix, "_utr_region_chisq_test.txt"))
  )

  plot_dt <- summary[, .(
    n_positive = sum(in_significant_qtl_region),
    n_total = .N,
    proportion = mean(in_significant_qtl_region)
  ), by = group]

  make_enrichment_plot(
    plot_dt,
    paste0(prefix, ": SNPs within significant QTL UTRs"),
    "Proportion of candidate SNPs",
    chisq_p,
    file.path(analysis_dir, paste0(prefix, "_utr_region_enrichment.pdf"))
  )

all_tested_region_snp_qq <- summary[
  in_any_qtl_region == TRUE &
    !is.na(min_region_p) &
    min_region_p > 0 &
    min_region_p <= 1
]

save_table(
  all_tested_region_snp_qq,
  file.path(
    analysis_dir,
    paste0(prefix, "_utr_region_snps_for_qqplot.tsv")
  )
)

cat(
  prefix, " candidate-SNP UTR-region QQ counts:\n",
  "  asAPA SNPs: ",
  nrow(all_tested_region_snp_qq[group == "asAPA"]), "\n",
  "  control SNPs: ",
  nrow(all_tested_region_snp_qq[group == "control"]), "\n",
  sep = ""
)

make_qq_plot(
  all_tested_region_snp_qq[group == "asAPA", min_region_p],
  all_tested_region_snp_qq[group == "control", min_region_p],
  paste0(prefix, ": candidate SNPs overlapping tested QTL UTRs"),
  file.path(
    analysis_dir,
    paste0(prefix, "_utr_region_qqplot.pdf")
  ),
  perform_ks = FALSE,
  subtitle_note = paste0(
    "One minimum p-value per candidate SNP; SNPs within the same UTR ",
    "may share p-values"
  )
)
all_tested_utr_qq <- hits[
  !is.na(nominal_p) & nominal_p > 0 & nominal_p <= 1,
  .(
    utr_pvalue = safe_min(nominal_p),
    n_overlapping_candidate_snps = uniqueN(SNP),
    overlapping_candidate_snps = paste(sort(unique(SNP)), collapse = ",")
  ),
  by = .(
    group,
    phenotype,
    gene,
    region_chr,
    region_start,
    region_end
  )
]

save_table(
  all_tested_utr_qq,
  file.path(analysis_dir, paste0(prefix, "_all_tested_utrs_for_qqplot.tsv"))
)

cat(
  prefix, " tested UTR QQ counts:\n",
  "  asAPA UTRs: ", nrow(all_tested_utr_qq[group == "asAPA"]), "\n",
  "  control UTRs: ", nrow(all_tested_utr_qq[group == "control"]), "\n",
  sep = ""
)

make_qq_plot(
  all_tested_utr_qq[group == "asAPA", utr_pvalue],
  all_tested_utr_qq[group == "control", utr_pvalue],
  paste0(prefix, ": unique tested QTL UTRs"),
  file.path(analysis_dir, paste0(prefix, "_unique_tested_utr_qqplot.pdf"))
)

all_tested_snp_utr_pairs <- hits[
  !is.na(nominal_p) & nominal_p > 0 & nominal_p <= 1,
  .(
    pair_pvalue = safe_min(nominal_p)
  ),
  by = .(
    group,
    SNP,
    tissue,
    phenotype,
    gene,
    region_chr,
    region_start,
    region_end
  )
]

save_table(
  all_tested_snp_utr_pairs,
  file.path(
    analysis_dir,
    paste0(prefix, "_all_tested_snp_utr_pairs_for_qqplot.tsv")
  )
)

make_qq_plot(
  all_tested_snp_utr_pairs[group == "asAPA", pair_pvalue],
  all_tested_snp_utr_pairs[group == "control", pair_pvalue],
  paste0(prefix, ": all tested SNP-UTR overlap pairs"),
  file.path(
    analysis_dir,
    paste0(prefix, "_all_tested_snp_utr_pair_qqplot.pdf")
  ),
  perform_ks = FALSE,
  subtitle_note = paste0(
    "Descriptive sensitivity analysis; repeated UTR p-values are ",
    "not independent"
  )
)

}

run_gene_analysis <- function(qtl, asapa_genes, control_genes,
                              analysis_dir, prefix) {

  # Collapse QTL results to one row per gene
  gene_summary <- qtl[, .(
    significant = any(significant, na.rm = TRUE),
    min_p = safe_min(nominal_p),
    min_q = safe_min(qvalue),
    max_abs_slope = safe_max_abs(slope)
  ), by = gene]

  # Build complete candidate gene sets
  candidates <- rbindlist(list(
    data.table(
      gene = unique(asapa_genes),
      group = "asAPA"
    ),
    data.table(
      gene = unique(control_genes),
      group = "control"
    )
  ))

  candidates[, group := factor(
    group,
    levels = c("asAPA", "control")
  )]

  overlap <- merge(
    candidates,
    gene_summary,
    by = "gene",
    all.x = TRUE
  )

  overlap[, tested := !is.na(min_p)]

  overlap[is.na(significant), significant := FALSE]

  save_table(
    overlap,
    file.path(
      analysis_dir,
      paste0(prefix, "_gene_overlap.tsv")
    )
  )
  # Gene enrichment using the entire candidate gene set
  enrichment_table <- table(
    factor(
      overlap$group,
      levels = c("asAPA", "control")
    ),
    factor(
      overlap$significant,
      levels = c(TRUE, FALSE)
    )
  )

  colnames(enrichment_table) <- c(
    "Significant",
    "Not_significant"
  )

  write.table(
    enrichment_table,
    file.path(
      analysis_dir,
      paste0(prefix, "_gene_enrichment_table.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  chisq_p <- safe_chisq(
    enrichment_table,
    file.path(
      analysis_dir,
      paste0(prefix, "_gene_chisq_test.txt")
    )
  )

  # Numerator:
  # Number of significant candidate genes
  #
  # Denominator:
  # Total number of candidate genes in that group
  plot_dt <- overlap[, .(
    n_positive = sum(significant, na.rm = TRUE),
    n_total = .N,
    proportion = sum(significant, na.rm = TRUE) / .N,
    n_tested = sum(tested, na.rm = TRUE)
  ), by = group]

  save_table(
    plot_dt,
    file.path(
      analysis_dir,
      paste0(prefix, "_gene_enrichment_plot_counts.tsv")
    )
  )

  make_enrichment_plot(
    plot_dt,
    paste0(prefix, ": significant genes in total candidate set"),
    "Proportion of total candidate genes significant",
    chisq_p,
    file.path(
      analysis_dir,
      paste0(prefix, "_gene_enrichment.pdf")
    )
  )

  # Gene QQ plot using all tested candidate genes
  tested_gene_dt <- overlap[
    tested == TRUE &
      !is.na(min_p) &
      min_p > 0 &
      min_p <= 1
  ]

  save_table(
    tested_gene_dt,
    file.path(
      analysis_dir,
      paste0(prefix, "_all_tested_candidate_genes_for_qqplot.tsv")
    )
  )

  make_qq_plot(
    tested_gene_dt[group == "asAPA", min_p],
    tested_gene_dt[group == "control", min_p],
    paste0(prefix, ": all tested candidate genes"),
    file.path(
      analysis_dir,
      paste0(prefix, "_all_tested_gene_qqplot.pdf")
    )
  )
}

# Load candidate SNP and gene sets
asapa_snps   <- load_candidate_snps(asapa_file, "asAPA")
control_snps <- load_candidate_snps(control_file, "control")
asapa_genes   <- load_candidate_genes(asapa_gene_file)
control_genes <- load_candidate_genes(control_gene_file)

sets <- prepare_candidate_sets(
  asapa_snps, control_snps, asapa_genes, control_genes
)

asapa_snps   <- sets$asapa_snps
control_snps <- sets$control_snps
asapa_genes   <- sets$asapa_genes
control_genes <- sets$control_genes

asapa_ratio <- load_asapa_ratio(asapa_ratio_file)

save_table(
  data.table(SNP = sets$shared_snps),
  file.path(outdir, "shared_snps_excluded.tsv")
)
save_table(
  data.table(gene = sets$shared_genes),
  file.path(outdir, "shared_genes_excluded.tsv")
)

# Load and combine tensorQTL tissue results
load_tensor <- function(file, tissue) {
  dt <- fread(file)

  # phenotype/UTR information is column 1.
  phenotype_col <- require_column_index(
    dt, 1, paste0(tissue, " phenotype/UTR column")
  )
  variant_col <- get_column(
    dt,
    c("variant_id", "Variant", "variant"),
    label = paste0(tissue, " variant column")
  )
  # Per the requested tensorQTL format, nominal p-value is column 13.
  p_col <- require_column_index(
    dt, 13, paste0(tissue, " nominal p-value column")
  )

  # Per the requested tensorQTL format, q-value is column 18.
  q_col <- require_column_index(dt, 18, paste0(tissue, " q-value column"))

  slope_col <- get_column(
    dt,
    c("slope", "beta", "effect", "Effect"),
    fallback_index = NULL,
    label = paste0(tissue, " effect-size column")
  )

  variant_parsed <- parse_variant_id(dt[[variant_col]])
  region_parsed <- parse_phenotype_region(dt[[phenotype_col]])

  dt[, phenotype := as.character(get(phenotype_col))]
  dt[, variant := as.character(get(variant_col))]
  dt[, nominal_p := suppressWarnings(as.numeric(get(p_col)))]
  dt[, qvalue := suppressWarnings(as.numeric(get(q_col)))]
  if (all(is.na(dt$qvalue))) {
    stop(tissue, " q-value column could not be converted to numeric values.")
  }
  dt[, slope := suppressWarnings(as.numeric(get(slope_col)))]
  dt[, tissue := tissue]
  dt[, SNP := variant_parsed$SNP]
  dt[, ref := variant_parsed$ref]
  dt[, alt := variant_parsed$alt]
  dt[, gene := extract_gene(phenotype)]
  dt[, region := region_parsed$region]
  dt[, region_chr := region_parsed$chr]
  dt[, region_start := region_parsed$start]
  dt[, region_end := region_parsed$end]
  dt[, significant := !is.na(qvalue) & qvalue < 0.05]

  dt
}

tensor_all <- rbindlist(
  list(
    load_tensor(FP, "FP"),
    load_tensor(IFG, "IFG"),
    load_tensor(PG, "PG"),
    load_tensor(STG, "STG")
  ),
  use.names = TRUE,
  fill = TRUE,
  idcol = FALSE
)

save_table(
  tensor_all,
  file.path(tensor_dir, "tensorQTL_all_tissues_combined.tsv")
)
save_table(
  tensor_all[significant == TRUE],
  file.path(tensor_dir, "tensorQTL_all_tissues_qvalue_lt_0.05.tsv")
)

run_direct_snp_analysis(
  tensor_all, asapa_snps, control_snps, tensor_dir, "tensorQTL"
)
run_region_snp_analysis(
  tensor_all, asapa_snps, control_snps, tensor_dir, "tensorQTL"
)
run_gene_analysis(
  tensor_all, asapa_genes, control_genes, tensor_dir, "tensorQTL"
)

make_tissue_aqtl_count_plot(
  qtl = tensor_all,
  output_file = file.path(
    tensor_dir,
    "tensorQTL_aQTL_counts_by_tissue.pdf"
  ),
  output_table = file.path(
    tensor_dir,
    "tensorQTL_aQTL_counts_by_tissue.tsv"
  )
)

make_effect_imbalance_plot(
  qtl = tensor_all,
  asapa_ratio = asapa_ratio,
  output_file = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance.pdf"
  ),
  output_table = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance.tsv"
  ),
  output_correlation = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance_correlation.txt"
  )
)

make_regional_effect_imbalance_plot(
  qtl = tensor_all,
  asapa_ratio = asapa_ratio,
  output_file = file.path(
    tensor_dir,
    "tensorQTL_regional_effect_vs_asAPA_imbalance.pdf"
  ),
  output_table = file.path(
    tensor_dir,
    "tensorQTL_regional_effect_vs_asAPA_imbalance.tsv"
  ),
  output_correlation = file.path(
    tensor_dir,
    "tensorQTL_regional_effect_vs_asAPA_imbalance_correlation.txt"
  )
)

make_effect_imbalance_plot_by_tissue(
  qtl = tensor_all,
  asapa_ratio = asapa_ratio,
  output_file = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance_by_tissue.pdf"
  ),
  output_table = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance_by_tissue.tsv"
  ),
  output_correlation = file.path(
    tensor_dir,
    "tensorQTL_effect_size_vs_asAPA_imbalance_by_tissue_correlations.tsv"
  )
)

# Complete
cat(
  "\nCompleted all tensorQTL analyses.\n",
  "tensorQTL output: ", tensor_dir, "\n",
  sep = ""
)
