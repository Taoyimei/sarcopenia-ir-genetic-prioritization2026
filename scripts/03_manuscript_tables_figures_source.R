## =============================================================================
## Manuscript-only figures and tables for the IR-sarcopenia project
##
## This is a post-processing script. It does not source or modify
## IR_sarcopenia.R and does not rerun the discovery analyses. It reads only the
## validated workbooks created by that pipeline and writes publication-facing
## figures plus one Excel workbook containing main and supplementary tables.
## =============================================================================


## 00. Packages and paths -----------------------------------------------------

required_packages <- c(
  "data.table", "dplyr", "ggplot2", "openxlsx", "patchwork",
  "scales", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Please install the following packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "."
result_dir <- file.path(project_dir, "results_extended")
output_dir <- file.path(project_dir, "manuscript_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  main = file.path(result_dir, "manuscript_analysis_tables.xlsx"),
  magic = file.path(result_dir, "MAGIC_IR_mechanistic_extension.xlsx"),
  snrna = file.path(result_dir, "GSE167186_snRNA_celltype_validation.xlsx"),
  clinical = file.path(result_dir, "clinical_translation_public_databases.xlsx"),
  sqtl = file.path(result_dir, "GTEx_v11_sQTL_TierA_analysis.xlsx"),
  reviewer7 = file.path(result_dir, "Reviewer7_external_validation_snRNA_supplement.xlsx")
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Required result workbook(s) missing:\n", paste(missing_files, collapse = "\n"))
}

candidate_order <- c("ABCC8", "MAPK1", "YWHAZ", "ZBTB7B", "SMAD3", "RXRA")
trait_labels <- c(
  ALM = "Appendicular lean mass",
  GRIP = "Grip strength",
  WALK = "Walking pace"
)

read_sheet <- function(file_key, sheet) {
  data.table::as.data.table(openxlsx::read.xlsx(input_files[[file_key]], sheet = sheet))
}

theme_manuscript <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#20262E"),
      plot.subtitle = ggplot2::element_text(colour = "#4B5563"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#E9EEF5", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

save_figure <- function(plot, stem, width, height) {
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_tmp <- tempfile(pattern = paste0(stem, "_"), fileext = ".png")
  ggplot2::ggsave(
    png_tmp, plot = plot, device = ragg::agg_png,
    width = width, height = height, units = "in", dpi = 600,
    background = "white"
  )
  if (!file.copy(png_tmp, png_file, overwrite = TRUE)) {
    stop("PNG output failed: ", png_file)
  }
  unlink(png_tmp)
  ggplot2::ggsave(
    pdf_file, plot = plot, device = grDevices::cairo_pdf,
    width = width, height = height, units = "in",
    bg = "white"
  )
  invisible(c(png_file, pdf_file))
}


## 01. Read validated analysis results ---------------------------------------

tier_a_coloc <- read_sheet("main", "04_coloc_TierA_prioritised_cand")
mr_evidence <- read_sheet("main", "07_Manuscript_Table2_evidence_h")
fusion <- read_sheet("main", "08_FUSION_independent_replicati")
susie <- read_sheet("main", "04b_coloc_SuSiE_TierA")
sn_localization <- read_sheet("snrna", "07_TierA_localization")
sn_de <- read_sheet("snrna", "08_TierA_pseudobulk_DE")
sqtl_mr <- read_sheet("sqtl", "06_targeted_sQTL_MR")
sqtl_summary <- read_sheet("sqtl", "07_layer_summary")
clinical <- read_sheet("clinical", "01_gene_evidence_matrix")

tier_a_mr <- mr_evidence[
  evidence_tier == "Tier_A" & gene_symbol %in% candidate_order
]
if (nrow(tier_a_mr) != 6L || uniqueN(tier_a_mr$gene_symbol) != 6L) {
  stop("Expected exactly six Tier A gene-trait results.")
}
tier_a_mr[, `:=`(
  gene_symbol = factor(gene_symbol, levels = rev(candidate_order)),
  trait_label = unname(trait_labels[trait])
)]


## 02. Figure 1: study workflow ----------------------------------------------

workflow_nodes <- data.table(
  x = 1:6,
  y = 1,
  label = c(
    "Candidate definition\n24 IR-related pathways\n771 genes",
    "Muscle cis-eQTL screen\nGTEx v8\n551 main gene-trait tests",
    "MR and colocalization\n84 FDR associations\n6 Tier A genes",
    "Genetic robustness\nFUSION replication\nSuSiE multi-signal coloc",
    "Biological resolution\nsnRNA-seq and GTEx v11 sQTL",
    "Translation context\nHPA, Open Targets,\nChEMBL and trials"
  ),
  stage = factor(
    c("Discovery", "Discovery", "Prioritization", "Validation", "Mechanism", "Translation"),
    levels = c("Discovery", "Prioritization", "Validation", "Mechanism", "Translation")
  )
)
workflow_edges <- data.table(
  x = workflow_nodes$x[-nrow(workflow_nodes)] + 0.39,
  xend = workflow_nodes$x[-1] - 0.39,
  y = 1,
  yend = 1
)

figure1 <- ggplot() +
  geom_segment(
    data = workflow_edges,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.7, colour = "#4B5563",
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  geom_tile(
    data = workflow_nodes,
    aes(x = x, y = y, fill = stage),
    width = 0.78, height = 0.70, colour = "#243447", linewidth = 0.5
  ) +
  geom_text(
    data = workflow_nodes,
    aes(x = x, y = y, label = label),
    size = 3.0, lineheight = 1.05, colour = "#17212B"
  ) +
  scale_fill_manual(values = c(
    Discovery = "#DCEAF7", Prioritization = "#F5DDA7",
    Validation = "#DCE6B5", Mechanism = "#F3D0CB",
    Translation = "#D8D5EA"
  )) +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.5, 1.5), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", colour = "#20262E"),
    plot.subtitle = element_text(colour = "#4B5563"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    plot.margin = margin(12, 12, 12, 12)
  )
save_figure(figure1, "Figure1_Study_workflow", 15, 4.2)


## 03. Figure 2: Tier A cis-MR forest ----------------------------------------

figure2 <- ggplot(
  tier_a_mr,
  aes(x = beta, y = gene_symbol, colour = trait_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280") +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    width = 0.16, linewidth = 0.7, orientation = "y"
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(
      x = ci_upper,
      label = paste0("  PP4=", sprintf("%.3f", PP.H4))
    ),
    hjust = 0, size = 3.0, colour = "#30343B", show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Appendicular lean mass" = "#2563A6",
    "Grip strength" = "#C46A1A",
    "Walking pace" = "#6B7D2A"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.28))) +
  labs(
    title = "Figure 2. Skeletal-muscle cis-MR estimates for Tier A genes",
    subtitle = "Points are MR estimates; bars are 95% confidence intervals; labels show primary-prior colocalization PP4",
    x = "MR effect on the corresponding sarcopenia-related trait (GWAS units)",
    y = NULL,
    colour = "Outcome"
  ) +
  guides(
    colour = guide_legend(
      title.position = "left",
      override.aes = list(linewidth = 1.2, size = 4.0)
    )
  ) +
  theme_manuscript(10) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10.5),
    legend.key.width = grid::unit(0.34, "in"),
    legend.key.height = grid::unit(0.18, "in"),
    legend.box.spacing = grid::unit(0.18, "in"),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
    plot.margin = margin(12, 12, 14, 12)
  )
save_figure(figure2, "Figure2_TierA_cisMR_forest", 8.5, 5.8)


## 04. Figure 3: cell localization, aging, and splicing ----------------------

susie_summary <- susie[, .(
  susie_max_PP4 = max(PP.H4.abf, na.rm = TRUE)
), by = gene_symbol]
sn_support <- sn_de[
  test_status == "Tested" & is.finite(candidate_FDR_global) &
    candidate_FDR_global < 0.05,
  .(snrna_supported = TRUE), by = .(gene_symbol = gene)
]
sqtl_support <- sqtl_summary[, .(
  gene_symbol,
  sqtl_supported = is.finite(sqtl_mr_fdr) & sqtl_mr_fdr < 0.05
)]

evidence_base <- data.table(gene_symbol = candidate_order)
evidence_base <- Reduce(
  function(x, y) merge(x, y, by = "gene_symbol", all.x = TRUE, sort = FALSE),
  list(
    evidence_base,
    fusion[, .(gene_symbol, fusion_replicated = replicated)],
    susie_summary,
    sn_support,
    sqtl_support,
    clinical[, .(gene_symbol, chembl_mechanism_count)]
  )
)
evidence_base[is.na(fusion_replicated), fusion_replicated := FALSE]
evidence_base[is.na(snrna_supported), snrna_supported := FALSE]
evidence_base[is.na(sqtl_supported), sqtl_supported := FALSE]
evidence_base[is.na(chembl_mechanism_count), chembl_mechanism_count := 0]

evidence_long <- rbindlist(list(
  evidence_base[, .(gene_symbol, evidence = "Primary cis-MR + coloc", score = 2L)],
  evidence_base[, .(
    gene_symbol, evidence = "FUSION replication",
    score = ifelse(fusion_replicated, 2L, 0L)
  )],
  evidence_base[, .(
    gene_symbol, evidence = "SuSiE shared signal",
    score = ifelse(is.finite(susie_max_PP4) & susie_max_PP4 >= 0.80, 2L, 0L)
  )],
  evidence_base[, .(
    gene_symbol, evidence = "snRNA age association",
    score = ifelse(snrna_supported, 1L, 0L)
  )],
  evidence_base[, .(
    gene_symbol, evidence = "sQTL-MR clue",
    score = ifelse(sqtl_supported, 1L, 0L)
  )],
  evidence_base[, .(
    gene_symbol, evidence = "ChEMBL mechanism",
    score = ifelse(chembl_mechanism_count > 0, 1L, 0L)
  )]
))
evidence_long[, `:=`(
  gene_symbol = factor(gene_symbol, levels = rev(candidate_order)),
  evidence = factor(evidence, levels = c(
    "Primary cis-MR + coloc", "FUSION replication",
    "SuSiE shared signal", "snRNA age association",
    "sQTL-MR clue", "ChEMBL mechanism"
  )),
  evidence_level = factor(
    score, levels = 0:2,
    labels = c("Not established", "Contextual", "Statistical support")
  )
)]

figure4 <- ggplot(evidence_long, aes(x = evidence, y = gene_symbol, fill = evidence_level)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = c("", "C", "S")[score + 1L]), size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c(
    "Not established" = "#ECEFF2",
    "Contextual" = "#E8C36A",
    "Statistical support" = "#4A86B8"
  )) +
  labs(x = NULL, y = NULL, fill = "Evidence level") +
  theme_manuscript(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )
save_figure(figure4, "Figure4_Multilayer_evidence_matrix", 10.5, 6.2)


## 05. Figure 3: cell localization, aging, and splicing ----------------------

cell_order <- c(
  "Fast skeletal fiber", "Slow skeletal fiber", "Satellite cell",
  "Fibro-adipogenic progenitor", "Endothelial cell",
  "Smooth muscle/pericyte", "Immune cell"
)
localization_plot_data <- copy(sn_localization)
localization_plot_data[, `:=`(
  gene = factor(gene, levels = rev(candidate_order)),
  cell_type = factor(cell_type, levels = cell_order)
)]

panel_a <- ggplot(
  localization_plot_data,
  aes(x = cell_type, y = gene, fill = sqrt(percent_detected))
) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient(
    low = "#F2F5F8", high = "#2563A6",
    labels = function(x) sprintf("%.0f", x^2)
  ) +
  labs(
    title = "A. Candidate-gene localization",
    subtitle = "GSE167186; fill represents nuclei with detected expression (%)",
    x = NULL, y = NULL, fill = "% detected"
  ) +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

fast_fiber_de <- sn_de[
  cell_type == "Fast skeletal fiber" & gene %in% candidate_order &
    test_status == "Tested"
]
fast_fiber_de[, `:=`(
  gene = factor(gene, levels = rev(candidate_order)),
  fdr_supported = candidate_FDR_global < 0.05
)]
panel_b <- ggplot(fast_fiber_de, aes(x = logFC_Old_vs_Young, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280") +
  geom_segment(
    aes(x = 0, xend = logFC_Old_vs_Young, yend = gene),
    colour = "#9CA3AF", linewidth = 0.7
  ) +
  geom_point(aes(fill = fdr_supported), shape = 21, size = 3.2, colour = "#243447") +
  scale_fill_manual(
    values = c(`FALSE` = "white", `TRUE` = "#C46A1A"),
    labels = c(`FALSE` = "Not FDR-supported", `TRUE` = "FDR < 0.05")
  ) +
  labs(
    title = "B. Fast-fiber pseudobulk expression",
    subtitle = "Old versus young donors; filled point indicates global candidate FDR < 0.05",
    x = "log2 fold change (old versus young)", y = NULL, fill = "FDR support"
  ) +
  theme_manuscript(9)

sqtl_plot_data <- copy(sqtl_mr)
sqtl_plot_data[, gene_symbol := factor(gene_symbol, levels = rev(c("YWHAZ", "SMAD3")))]
panel_c <- ggplot(sqtl_plot_data, aes(x = beta, y = gene_symbol)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280") +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    width = 0.14, linewidth = 0.7, orientation = "y"
  ) +
  geom_point(
    aes(fill = fdr_targeted < 0.05), shape = 21,
    size = 3.2, colour = "#243447"
  ) +
  scale_fill_manual(
    values = c(`FALSE` = "white", `TRUE` = "#4A86B8"),
    labels = c(`FALSE` = "Not FDR-supported", `TRUE` = "FDR < 0.05")
  ) +
  labs(
    title = "C. Targeted skeletal-muscle sQTL-MR",
    subtitle = "Wald-ratio estimates for the prespecified gene-trait pairs",
    x = "MR effect on the corresponding trait (GWAS units)", y = NULL,
    fill = "FDR support"
  ) +
  theme_manuscript(9)

figure3 <- panel_a / (panel_b | panel_c)
save_figure(figure3, "Figure3_Cell_and_splicing_context", 12, 9.2)


## 06. Main manuscript tables ------------------------------------------------

table1_sources <- data.frame(
  Analysis_layer = c(
    "Candidate genes", "Muscle eQTL", "Outcome GWAS", "Outcome GWAS",
    "Outcome GWAS", "Independent muscle eQTL",
    "snRNA-seq", "Muscle sQTL",
    "Translation annotation"
  ),
  Source = c(
    "MSigDB", "GTEx v8 Muscle_Skeletal", "GCST90000025",
    "UKB-b-10215", "UKB-b-4711",
    "FUSION skeletal muscle", "GSE167186",
    "GTEx v11 Muscle_Skeletal", "HPA; Open Targets; ChEMBL; ClinicalTrials.gov"
  ),
  Population_or_tissue = c(
    "24 prespecified IR-related pathways", "Human skeletal muscle",
    "Appendicular lean mass", "Right-hand grip strength", "Walking pace",
    "Human skeletal muscle", "Human vastus lateralis",
    "Human skeletal muscle",
    "Human protein, target, drug and trial records"
  ),
  Sample_size_or_scope = c(
    "771 unique genes", "706 samples", "450,243", "461,089", "459,915",
    "301 samples",
    "17 snRNA donors; approximately 143,000 nuclei",
    "Sample size not encoded in the supplied sGenes summary file",
    "Six Tier A genes"
  ),
  Primary_use = c(
    "Prespecified discovery space", "cis-eQTL instruments", "MR outcome",
    "MR outcome", "MR outcome", "Independent replication",
    "Cell localization and donor-level pseudobulk",
    "Splicing-layer MR", "Translation prioritization only"
  ),
  stringsAsFactors = FALSE
)

table2_tier_a <- tier_a_mr[, .(
  Gene = as.character(gene_symbol),
  Outcome = trait_label,
  Method = method,
  SNPs = nsnp,
  MR_beta = beta,
  MR_SE = se,
  CI_lower_95 = ci_lower,
  CI_upper_95 = ci_upper,
  P_value = pval,
  FDR_within_trait = fdr_within_trait,
  Minimum_F = min_F,
  Coloc_SNPs = n_coloc_snps,
  Coloc_PP3 = PP.H3,
  Coloc_PP4 = PP.H4
)]
table2_tier_a[, Gene := factor(Gene, levels = candidate_order)]
setorder(table2_tier_a, Gene)
table2_tier_a[, Gene := as.character(Gene)]

table3_multilayer <- data.table(gene_symbol = candidate_order)
table3_multilayer <- Reduce(
  function(x, y) merge(x, y, by = "gene_symbol", all.x = TRUE, sort = FALSE),
  list(
    table3_multilayer,
    tier_a_coloc[, .(gene_symbol, trait, primary_PP4 = PP.H4)],
    fusion[, .(
      gene_symbol, fusion_direction_concordant = direction_concordant,
      fusion_PP4 = fusion_PP.H4, fusion_replicated = replicated
    )],
    susie_summary,
    sn_support,
    sqtl_summary[, .(
      gene_symbol, sgene_qval, sqtl_mr_beta, sqtl_mr_pval,
      sqtl_mr_fdr, integrated_inference
    )],
    clinical[, .(
      gene_symbol, hpa_skeletal_muscle_protein_reported,
      max_relevant_ot_association, open_targets_supported_tractability_items,
      chembl_mechanism_count, highest_chembl_phase,
      directionally_concordant_mechanisms, translation_evidence_level
    )]
  )
)
table3_multilayer[is.na(snrna_supported), snrna_supported := FALSE]
table3_multilayer[, gene_order := match(gene_symbol, candidate_order)]
setorder(table3_multilayer, gene_order)
table3_multilayer[, gene_order := NULL]


## 07. Supplementary tables retained for reproducibility --------------------

main_traits <- c("ALM", "GRIP", "WALK")
filter_main_traits <- function(x) {
  trait_columns <- intersect(c("trait", "sarcopenia_trait", "outcome"), names(x))
  for (col in trait_columns) {
    x <- x[get(col) %in% main_traits | !get(col) %in% c("LOW_GRIP")]
  }
  x
}
magic_two_step_main <- read_sheet("magic", "07_two_step_evidence_chain")
magic_two_step_main <- magic_two_step_main[sarcopenia_trait %in% main_traits]


## 07b. BEGIN NEW CROSS-PHENOTYPE MAGMA MODULE ------------------------------
##
## Added for Experimental Gerontology revision. This module formally compares
## the primary MAGMA enrichment coefficients across ALM, grip strength, and
## walking pace, instead of inferring between-trait differences from one
## significant and two nonsignificant enrichment tests.

magma_dir <- file.path(result_dir, "MAGMA")
magma_traits <- c("ALM", "GRIP", "WALK")
magma_trait_labels <- c(
  ALM = "Appendicular lean mass",
  GRIP = "Right-hand grip strength",
  WALK = "Self-reported usual walking pace"
)

fit_binary_set_contrast <- function(y, is_set) {
  ok <- is.finite(y) & !is.na(is_set)
  y <- as.numeric(y[ok])
  is_set <- as.logical(is_set[ok])
  n_set <- sum(is_set)
  n_background <- sum(!is_set)
  if (n_set < 2L || n_background < 2L) {
    return(list(
      beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_,
      n_genes = length(y), n_set = n_set, n_background = n_background
    ))
  }
  beta_background <- mean(y[!is_set])
  beta <- mean(y[is_set]) - beta_background
  fitted <- beta_background + beta * as.numeric(is_set)
  residual <- y - fitted
  sigma2 <- sum(residual^2) / (length(y) - 2L)
  se <- sqrt(sigma2 * (1 / n_set + 1 / n_background))
  t_value <- beta / se
  p_value <- 2 * stats::pt(abs(t_value), df = length(y) - 2L, lower.tail = FALSE)
  list(
    beta = beta, se = se, t = t_value, p = p_value,
    n_genes = length(y), n_set = n_set, n_background = n_background
  )
}

run_magma_cross_phenotype_comparison <- function(
    magma_dir,
    traits,
    trait_labels,
    trait_n,
    magma_exe = file.path(project_dir, "magma.exe")) {
  summary_file <- file.path(magma_dir, "MAGMA_competitive_results.tsv")
  set_file <- file.path(magma_dir, "IR_gene_sets_entrez.txt")
  pair_table <- data.table::data.table(
    trait_1 = c("ALM", "ALM", "GRIP"),
    trait_2 = c("GRIP", "WALK", "WALK")
  )
  required_inputs <- c(
    magma_exe,
    set_file,
    file.path(magma_dir, paste0(traits, "_gene.genes.raw"))
  )
  if (any(!file.exists(required_inputs))) {
    stop(
      "Required MAGMA input(s) missing:\n",
      paste(required_inputs[!file.exists(required_inputs)], collapse = "\n")
    )
  }

  magma_cli_path <- function(path, mustWork = TRUE) {
    if (mustWork) {
      return(utils::shortPathName(normalizePath(
        path, winslash = "\\", mustWork = TRUE
      )))
    }
    parent <- utils::shortPathName(normalizePath(
      dirname(path), winslash = "\\", mustWork = TRUE
    ))
    file.path(parent, basename(path))
  }

  run_magma_gsa <- function(trait, out_prefix, include_file = NULL,
                            gene_info = FALSE) {
    args <- c(
      "--gene-results",
      magma_cli_path(file.path(magma_dir, paste0(trait, "_gene.genes.raw"))),
      "--set-annot", magma_cli_path(set_file), "col=2,1"
    )
    settings <- character()
    if (!is.null(include_file)) {
      settings <- c(
        settings,
        paste0("gene-include=", magma_cli_path(include_file))
      )
    }
    if (gene_info) settings <- c(settings, "gene-info")
    if (length(settings) > 0L) args <- c(args, "--settings", settings)
    args <- c(
      args,
      "--out", magma_cli_path(out_prefix, mustWork = FALSE)
    )
    magma_log <- system2(
      magma_cli_path(magma_exe),
      args = args,
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(magma_log, "status")
    if (is.null(status)) status <- 0L
    if (status != 0L) {
      stop(
        "MAGMA gene-set analysis failed for ", trait, ":\n",
        paste(tail(magma_log, 20L), collapse = "\n")
      )
    }
    invisible(out_prefix)
  }

  read_gsa <- function(prefix, trait) {
    file <- paste0(prefix, ".gsa.out.txt")
    x <- data.table::fread(file, skip = "VARIABLE", header = TRUE)
    x[, trait := trait]
    x[]
  }

  read_ir_full <- function(prefix, trait, excluded_chr = NA_integer_) {
    x <- read_gsa(prefix, trait)
    is_target <- x$VARIABLE == "IR_FULL_771"
    if ("FULL_NAME" %in% names(x)) {
      is_target <- is_target | x$FULL_NAME == "IR_FULL_771"
    }
    x <- x[is_target]
    if (nrow(x) != 1L) {
      stop("IR_FULL_771 result not uniquely identified for ", trait)
    }
    x[, .(
      trait,
      excluded_chr,
      NGENES,
      BETA,
      BETA_STD,
      SE,
      P
    )]
  }

  # Re-run all three full gene-set analyses with identical output settings.
  for (trait in traits) {
    message("MAGMA full gene-set analysis with gene residuals: ", trait)
    run_magma_gsa(
      trait,
      file.path(magma_dir, paste0(trait, "_IR_competitive")),
      gene_info = TRUE
    )
  }

  magma_summary <- data.table::rbindlist(lapply(traits, function(trait) {
    read_gsa(file.path(magma_dir, paste0(trait, "_IR_competitive")), trait)
  }), fill = TRUE)
  magma_summary[, analysis_family := data.table::fifelse(
    VARIABLE == "IR_FULL_771" | FULL_NAME == "IR_FULL_771",
    "primary_IR_total",
    "secondary_pathways"
  )]
  magma_summary[, FDR := stats::p.adjust(P, method = "BH"), by = analysis_family]
  data.table::setcolorder(
    magma_summary,
    c(
      "trait", "VARIABLE", "analysis_family", "NGENES", "BETA",
      "BETA_STD", "SE", "P", "FDR", "TYPE", "FULL_NAME"
    )
  )
  data.table::fwrite(magma_summary, summary_file, sep = "\t")

  beta_summary <- magma_summary[
    analysis_family == "primary_IR_total" & trait %in% traits,
    .(
      trait,
      trait_label = unname(trait_labels[trait]),
      sample_size = unname(trait_n[trait]),
      NGENES,
      BETA,
      BETA_STD,
      SE,
      P,
      FDR
    )
  ]
  beta_summary[, trait_order := match(trait, traits)]
  data.table::setorder(beta_summary, trait_order)
  beta_summary[, trait_order := NULL]

  naive_cross_trait <- data.table::rbindlist(lapply(seq_len(nrow(pair_table)), function(i) {
    left <- beta_summary[trait == pair_table$trait_1[i]]
    right <- beta_summary[trait == pair_table$trait_2[i]]
    if (nrow(left) != 1L || nrow(right) != 1L) {
      return(data.table::data.table(
        comparison = paste(pair_table$trait_1[i], "vs", pair_table$trait_2[i]),
        trait_1 = pair_table$trait_1[i],
        trait_2 = pair_table$trait_2[i],
        status = "not_run_missing_beta_summary",
        beta_1 = NA_real_, beta_2 = NA_real_, delta_beta = NA_real_,
        se_delta = NA_real_, z = NA_real_, p_value = NA_real_
      ))
    }
    delta_beta <- left$BETA - right$BETA
    se_delta <- sqrt(left$SE^2 + right$SE^2)
    z_value <- delta_beta / se_delta
    data.table::data.table(
      comparison = paste(pair_table$trait_1[i], "vs", pair_table$trait_2[i]),
      trait_1 = pair_table$trait_1[i],
      trait_2 = pair_table$trait_2[i],
      status = "complete_naive_zero_covariance",
      beta_1 = left$BETA,
      beta_2 = right$BETA,
      delta_beta = delta_beta,
      se_delta = se_delta,
      z = z_value,
      p_value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
    )
  }), fill = TRUE)
  naive_cross_trait[, fdr_BH := stats::p.adjust(p_value, method = "BH")]
  naive_cross_trait[, interpretation_boundary := paste0(
    "Naive comparison of published MAGMA coefficients. It assumes zero ",
    "cross-outcome covariance and is therefore not the primary ",
    "sample-overlap-aware result."
  )]

  covariance_sensitivity <- data.table::rbindlist(lapply(seq_len(nrow(pair_table)), function(i) {
    left <- beta_summary[trait == pair_table$trait_1[i]]
    right <- beta_summary[trait == pair_table$trait_2[i]]
    if (nrow(left) != 1L || nrow(right) != 1L) return(NULL)
    rho_grid <- seq(-0.50, 0.90, by = 0.10)
    out <- data.table::data.table(
      comparison = paste(pair_table$trait_1[i], "vs", pair_table$trait_2[i]),
      trait_1 = pair_table$trait_1[i],
      trait_2 = pair_table$trait_2[i],
      assumed_correlation = rho_grid
    )
    out[, delta_beta := left$BETA - right$BETA]
    out[, se_delta := sqrt(left$SE^2 + right$SE^2 -
                             2 * assumed_correlation * left$SE * right$SE)]
    out[, z := delta_beta / se_delta]
    out[, p_value := 2 * stats::pnorm(abs(z), lower.tail = FALSE)]
    out
  }), fill = TRUE)

  gene_sets <- data.table::fread(set_file, header = FALSE)
  data.table::setnames(gene_sets, c("gene_set", "GENE"))
  candidate_entrez <- unique(as.character(gene_sets[gene_set == "IR_FULL_771", GENE]))

  residual_tables <- list()
  for (trait in traits) {
    residual_file <- file.path(magma_dir, paste0(trait, "_IR_competitive.gsa.genes.out.txt"))
    if (!file.exists(residual_file)) {
      stop("MAGMA gene-residual output missing after re-run: ", residual_file)
    }
    residual_table <- data.table::fread(residual_file, skip = "GENE")
    residual_required <- c("GENE", "CHR", "ZRESID_BASE")
    missing_residual_columns <- setdiff(residual_required, names(residual_table))
    if (length(missing_residual_columns) > 0L) {
      stop(
        trait, " MAGMA gene-residual file is missing required column(s): ",
        paste(missing_residual_columns, collapse = ", ")
      )
    }
    residual_tables[[trait]] <- residual_table[, .(
      GENE = as.character(GENE),
      CHR = as.character(CHR),
      zresid = as.numeric(ZRESID_BASE)
    )]
  }

  common_genes <- Reduce(
    intersect,
    lapply(residual_tables, function(x) x[is.finite(zresid), GENE])
  )
  common_map <- residual_tables[[traits[1L]]][
    GENE %in% common_genes,
    .(GENE, CHR)
  ]
  common_map <- unique(common_map, by = "GENE")
  common_map[, CHR := as.integer(CHR)]
  data.table::setorder(common_map, CHR, GENE)
  for (trait in traits[-1L]) {
    check_map <- residual_tables[[trait]][
      GENE %in% common_genes,
      .(GENE, CHR_check = as.integer(CHR))
    ]
    check_map <- unique(check_map, by = "GENE")
    chromosome_check <- merge(common_map, check_map, by = "GENE")
    if (nrow(chromosome_check) != length(common_genes) ||
        any(chromosome_check$CHR != chromosome_check$CHR_check)) {
      stop("Cross-trait chromosome assignment mismatch for ", trait)
    }
  }
  blocks <- sort(unique(common_map$CHR))
  if (!identical(blocks, 1:22)) {
    stop("Expected autosomal chromosome blocks 1-22; observed: ",
         paste(blocks, collapse = ", "))
  }

  cross_trait_tmp <- tempfile(pattern = "magma_cross_trait_")
  dir.create(cross_trait_tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cross_trait_tmp, recursive = TRUE, force = TRUE), add = TRUE)

  common_include_file <- file.path(cross_trait_tmp, "common_genes.txt")
  writeLines(common_map$GENE, common_include_file)
  common_beta <- data.table::rbindlist(lapply(traits, function(trait) {
    message("MAGMA common-universe analysis: ", trait)
    prefix <- file.path(cross_trait_tmp, paste0(trait, "_common"))
    run_magma_gsa(trait, prefix, include_file = common_include_file)
    read_ir_full(prefix, trait)
  }))
  common_beta[, `:=`(
    trait_label = unname(trait_labels[trait]),
    sample_size = unname(trait_n[trait]),
    gene_universe = "three-trait common gene universe",
    n_total_genes = nrow(common_map),
    n_candidate_genes = sum(common_map$GENE %in% candidate_entrez)
  )]

  checkpoint_file <- file.path(output_dir, ".MAGMA_LOCO_checkpoint.tsv")
  checkpoint_required <- c(
    "trait", "excluded_chr", "NGENES", "BETA", "BETA_STD", "SE", "P",
    "n_common_genes", "n_common_candidate_genes"
  )
  if (file.exists(checkpoint_file)) {
    loco_checkpoint <- data.table::fread(checkpoint_file)
    if (length(setdiff(checkpoint_required, names(loco_checkpoint))) > 0L) {
      stop("Existing MAGMA LOCO checkpoint has an incompatible format.")
    }
    loco_checkpoint <- loco_checkpoint[
      trait %in% traits & excluded_chr %in% blocks &
        n_common_genes == nrow(common_map) &
        n_common_candidate_genes == sum(common_map$GENE %in% candidate_entrez)
    ]
    loco_checkpoint <- unique(loco_checkpoint, by = c("trait", "excluded_chr"))
  } else {
    loco_checkpoint <- data.table::data.table(
      trait = character(),
      excluded_chr = integer(),
      NGENES = integer(),
      BETA = numeric(),
      BETA_STD = numeric(),
      SE = numeric(),
      P = numeric(),
      n_common_genes = integer(),
      n_common_candidate_genes = integer()
    )
  }

  loco_estimates <- data.table::rbindlist(lapply(blocks, function(block) {
    completed <- loco_checkpoint[excluded_chr == block]
    missing_traits <- setdiff(traits, completed$trait)
    message(
      "MAGMA leave-one-chromosome analysis: chromosome ", block,
      " of 22; remaining traits = ", length(missing_traits)
    )
    include_file <- file.path(cross_trait_tmp, paste0("common_no_chr", block, ".txt"))
    writeLines(common_map[CHR != block, GENE], include_file)
    newly_completed <- data.table::rbindlist(lapply(missing_traits, function(trait) {
      prefix <- file.path(
        cross_trait_tmp,
        paste0(trait, "_common_no_chr", block)
      )
      run_magma_gsa(trait, prefix, include_file = include_file)
      read_ir_full(prefix, trait, excluded_chr = block)
    }), fill = TRUE)
    if (nrow(newly_completed) > 0L) {
      newly_completed[, `:=`(
        n_common_genes = nrow(common_map),
        n_common_candidate_genes = sum(common_map$GENE %in% candidate_entrez)
      )]
      data.table::fwrite(
        newly_completed,
        checkpoint_file,
        sep = "\t",
        append = file.exists(checkpoint_file),
        col.names = !file.exists(checkpoint_file)
      )
      loco_checkpoint <- data.table::rbindlist(
        list(loco_checkpoint, newly_completed),
        fill = TRUE
      )
    }
    data.table::rbindlist(list(completed, newly_completed), fill = TRUE)
  }))
  loco_estimates <- unique(loco_estimates, by = c("trait", "excluded_chr"))
  if (nrow(loco_estimates) != length(traits) * length(blocks)) {
    stop("MAGMA LOCO analysis is incomplete after checkpoint recovery.")
  }
  block_counts <- common_map[, .(
    n_genes_excluded = .N,
    n_candidate_genes_excluded = sum(GENE %in% candidate_entrez)
  ), by = .(excluded_chr = CHR)]
  loco_estimates <- merge(
    loco_estimates,
    block_counts,
    by = "excluded_chr",
    all.x = TRUE,
    sort = FALSE
  )
  loco_estimates[, `:=`(
    trait_label = unname(trait_labels[trait]),
    n_genes_retained = nrow(common_map) - n_genes_excluded,
    n_candidate_genes_retained =
      sum(common_map$GENE %in% candidate_entrez) - n_candidate_genes_excluded
  )]

  primary_loco <- data.table::rbindlist(lapply(seq_len(nrow(pair_table)), function(i) {
    trait_1 <- pair_table$trait_1[i]
    trait_2 <- pair_table$trait_2[i]
    full_1 <- common_beta[trait == trait_1]
    full_2 <- common_beta[trait == trait_2]
    loo <- merge(
      loco_estimates[trait == trait_1, .(excluded_chr, beta_1_loo = BETA)],
      loco_estimates[trait == trait_2, .(excluded_chr, beta_2_loo = BETA)],
      by = "excluded_chr"
    )
    loo[, delta_beta_loo := beta_1_loo - beta_2_loo]
    k <- nrow(loo)
    delta_beta <- full_1$BETA - full_2$BETA
    se_delta <- sqrt(
      (k - 1) / k *
        sum((loo$delta_beta_loo - mean(loo$delta_beta_loo))^2)
    )
    var_1 <- (k - 1) / k *
      sum((loo$beta_1_loo - mean(loo$beta_1_loo))^2)
    var_2 <- (k - 1) / k *
      sum((loo$beta_2_loo - mean(loo$beta_2_loo))^2)
    covariance <- (k - 1) / k * sum(
      (loo$beta_1_loo - mean(loo$beta_1_loo)) *
        (loo$beta_2_loo - mean(loo$beta_2_loo))
    )
    statistic <- delta_beta / se_delta
    critical <- stats::qt(0.975, df = k - 1L)
    data.table::data.table(
      comparison = paste(trait_1, "vs", trait_2),
      trait_1,
      trait_2,
      beta_1 = full_1$BETA,
      beta_2 = full_2$BETA,
      delta_beta,
      jackknife_covariance = covariance,
      jackknife_correlation = covariance / sqrt(var_1 * var_2),
      se_delta_jackknife = se_delta,
      t_statistic = statistic,
      degrees_of_freedom = k - 1L,
      p_value = 2 * stats::pt(abs(statistic), df = k - 1L, lower.tail = FALSE),
      p_value_normal_reference = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
      ci_95_lower = delta_beta - critical * se_delta,
      ci_95_upper = delta_beta + critical * se_delta,
      n_common_genes = nrow(common_map),
      n_candidate_genes = sum(common_map$GENE %in% candidate_entrez),
      n_blocks = k,
      status = "complete_primary_common_universe_MAGMA_LOCO"
    )
  }))
  primary_loco[, fdr_BH := stats::p.adjust(p_value, method = "BH")]
  primary_loco[, interpretation_boundary := paste0(
    "Primary two-sided cross-phenotype contrast of IR_FULL_771 MAGMA BETA ",
    "coefficients in the same gene universe. Cross-trait covariance and the ",
    "difference SE were estimated from 22 paired leave-one-chromosome-out ",
    "replicates. This is a summary-statistic block analysis and does not ",
    "identify the exact number of overlapping participants."
  )]

  residual_wide <- Reduce(function(left, right) {
    merge(left, right, by = c("GENE", "CHR"), all = FALSE)
  }, lapply(traits, function(trait) {
    x <- data.table::copy(residual_tables[[trait]][GENE %in% common_genes])
    data.table::setnames(x, "zresid", trait)
    x
  }))
  residual_wide[, is_candidate := GENE %in% candidate_entrez]

  transform_residuals <- function(x, scale_name) {
    if (scale_name == "raw_ZRESID_BASE") return(x)
    if (scale_name == "within_trait_z_standardized") {
      return(as.numeric(scale(x)))
    }
    if (scale_name == "within_trait_rank_normalized") {
      return(stats::qnorm((rank(x, ties.method = "average") - 0.5) / length(x)))
    }
    stop("Unknown residual scale: ", scale_name)
  }

  residual_scales <- c(
    "raw_ZRESID_BASE",
    "within_trait_z_standardized",
    "within_trait_rank_normalized"
  )
  residual_sensitivity <- data.table::rbindlist(lapply(residual_scales, function(scale_name) {
    transformed <- data.table::copy(residual_wide)
    for (trait in traits) {
      transformed[, (trait) := transform_residuals(get(trait), scale_name)]
    }
    data.table::rbindlist(lapply(seq_len(nrow(pair_table)), function(i) {
      trait_1 <- pair_table$trait_1[i]
      trait_2 <- pair_table$trait_2[i]
      transformed[, residual_difference := get(trait_1) - get(trait_2)]
      full_fit <- fit_binary_set_contrast(
        transformed$residual_difference,
        transformed$is_candidate
      )
      leave_one_block <- vapply(blocks, function(block) {
        fit_binary_set_contrast(
          transformed[CHR != block, residual_difference],
          transformed[CHR != block, is_candidate]
        )$beta
      }, numeric(1))
      k <- length(leave_one_block)
      se_delta <- sqrt(
        (k - 1) / k *
          sum((leave_one_block - mean(leave_one_block))^2)
      )
      statistic <- full_fit$beta / se_delta
      critical <- stats::qt(0.975, df = k - 1L)
      data.table::data.table(
        scale = scale_name,
        comparison = paste(trait_1, "vs", trait_2),
        trait_1,
        trait_2,
        residual_contrast_beta = full_fit$beta,
        se_delta_jackknife = se_delta,
        t_statistic = statistic,
        degrees_of_freedom = k - 1L,
        p_value = 2 * stats::pt(abs(statistic), df = k - 1L, lower.tail = FALSE),
        p_value_normal_reference = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
        ci_95_lower = full_fit$beta - critical * se_delta,
        ci_95_upper = full_fit$beta + critical * se_delta,
        n_common_genes = nrow(transformed),
        n_candidate_genes = sum(transformed$is_candidate),
        n_blocks = k
      )
    }))
  }))
  residual_sensitivity[, fdr_BH := stats::p.adjust(p_value, method = "BH"), by = scale]
  residual_sensitivity[, interpretation_boundary := paste0(
    "Sensitivity analysis based on paired residualized gene statistics. ",
    "Within-trait standardization and rank normalization address residual-scale ",
    "differences but are not formal SNP-heritability adjustment."
  )]

  residual_correlations <- data.table::rbindlist(lapply(seq_len(nrow(pair_table)), function(i) {
    trait_1 <- pair_table$trait_1[i]
    trait_2 <- pair_table$trait_2[i]
    background <- !residual_wide$is_candidate
    data.table::data.table(
      comparison = paste(trait_1, "vs", trait_2),
      trait_1,
      trait_2,
      pearson_all_genes = stats::cor(residual_wide[[trait_1]], residual_wide[[trait_2]]),
      pearson_background_genes = stats::cor(
        residual_wide[[trait_1]][background],
        residual_wide[[trait_2]][background]
      ),
      spearman_all_genes = stats::cor(
        residual_wide[[trait_1]], residual_wide[[trait_2]], method = "spearman"
      ),
      spearman_background_genes = stats::cor(
        residual_wide[[trait_1]][background],
        residual_wide[[trait_2]][background], method = "spearman"
      ),
      n_common_genes = nrow(residual_wide),
      n_background_genes = sum(background)
    )
  }))

  input_qc <- data.table::rbindlist(lapply(traits, function(trait) {
    data.table::data.table(
      trait,
      trait_label = unname(trait_labels[trait]),
      sample_size = unname(trait_n[trait]),
      gene_result_file = paste0(trait, "_gene.genes.raw"),
      residual_file = paste0(trait, "_IR_competitive.gsa.genes.out.txt"),
      n_trait_genes = nrow(residual_tables[[trait]]),
      n_three_trait_common_genes = nrow(common_map),
      n_common_candidate_genes = sum(common_map$GENE %in% candidate_entrez),
      n_completed_LOCO_runs = sum(loco_estimates[["trait"]] == trait),
      status = "complete"
    )
  }))

  analysis_notes <- data.table::data.table(
    item = c(
      "Primary comparison",
      "Covariance handling",
      "Multiple testing",
      "Residual sensitivity",
      "Sample overlap boundary",
      "Power and heritability boundary",
      "Software",
      "Generated"
    ),
    value = c(
      paste(
        "Two-sided paired contrasts of IR_FULL_771 MAGMA BETA coefficients",
        "within the same", format(nrow(common_map), big.mark = ","),
        "gene universe using 22 leave-one-chromosome-out replicates."
      ),
      paste(
        "The jackknife covariance between each pair of trait-specific MAGMA",
        "coefficients is estimated from paired chromosome deletions."
      ),
      "Benjamini-Hochberg correction across the three prespecified pairwise contrasts.",
      paste(
        "Raw, within-trait Z-standardized, and within-trait rank-normalized",
        "gene-residual contrasts are reported as sensitivity analyses."
      ),
      paste(
        "The paired block analysis accommodates cross-trait covariance in summary",
        "statistics but does not identify the exact number of overlapping participants."
      ),
      paste(
        "Sample sizes and scale-standardized sensitivities are reported; these",
        "analyses do not constitute formal SNP-heritability adjustment."
      ),
      paste0("MAGMA v1.10; R ", getRversion()),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )

  list(
    beta_summary = beta_summary,
    primary_loco = primary_loco,
    common_beta = common_beta,
    loco_estimates = loco_estimates,
    residual_sensitivity = residual_sensitivity,
    residual_correlations = residual_correlations,
    naive_cross_trait = naive_cross_trait,
    covariance_sensitivity = covariance_sensitivity,
    input_qc = input_qc,
    analysis_notes = analysis_notes,
    checkpoint_file = checkpoint_file
  )
}

cross_trait_workbook_file <- file.path(
  output_dir,
  "MAGMA_cross_phenotype_comparison.xlsx"
)
magma_cross_trait <- if (
    identical(Sys.getenv("REUSE_MAGMA_CROSS_TRAIT"), "1") &&
      file.exists(cross_trait_workbook_file)) {
  list(
    beta_summary = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "01_single_trait")),
    primary_loco = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "02_primary_MAGMA_LOCO")),
    common_beta = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "03_common_universe_beta")),
    loco_estimates = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "04_LOCO_trait_estimates")),
    residual_sensitivity = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "05_residual_sensitivity")),
    residual_correlations = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "06_residual_correlations")),
    naive_cross_trait = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "07_zero_covariance")),
    covariance_sensitivity = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "08_assumed_rho")),
    input_qc = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "09_input_QC")),
    analysis_notes = data.table::as.data.table(openxlsx::read.xlsx(cross_trait_workbook_file, "10_analysis_notes")),
    checkpoint_file = NA_character_
  )
} else {
  run_magma_cross_phenotype_comparison(
    magma_dir,
    magma_traits,
    magma_trait_labels,
    trait_n = c(ALM = 450243, GRIP = 461089, WALK = 459915)
  )
}
if (dir.exists(dirname(cross_trait_workbook_file))) {
  openxlsx::write.xlsx(
    list(
      `01_single_trait` = as.data.frame(magma_cross_trait$beta_summary),
      `02_primary_MAGMA_LOCO` = as.data.frame(magma_cross_trait$primary_loco),
      `03_common_universe_beta` = as.data.frame(magma_cross_trait$common_beta),
      `04_LOCO_trait_estimates` = as.data.frame(magma_cross_trait$loco_estimates),
      `05_residual_sensitivity` = as.data.frame(magma_cross_trait$residual_sensitivity),
      `06_residual_correlations` = as.data.frame(magma_cross_trait$residual_correlations),
      `07_zero_covariance` = as.data.frame(magma_cross_trait$naive_cross_trait),
      `08_assumed_rho` = as.data.frame(magma_cross_trait$covariance_sensitivity),
      `09_input_QC` = as.data.frame(magma_cross_trait$input_qc),
      `10_analysis_notes` = as.data.frame(magma_cross_trait$analysis_notes)
    ),
    file = cross_trait_workbook_file,
    overwrite = TRUE,
    asTable = TRUE
  )
  if (file.exists(magma_cross_trait$checkpoint_file)) {
    unlink(magma_cross_trait$checkpoint_file, force = TRUE)
  }
}

## 07b. END NEW CROSS-PHENOTYPE MAGMA MODULE --------------------------------


## 07c. BEGIN NEW CANDIDATE-SET SPECIFICITY SENSITIVITY MODULE ---------------
##
## Added for Experimental Gerontology revision. This module evaluates whether
## the ALM enrichment of the IR-related 771-gene set is robust to pathway-set
## definition, broad pathway exclusion, leave-one-pathway-out analyses, and
## matched random gene-set comparisons.

candidate_specificity_config <- list(
  target_trait = "ALM",
  random_seed = 20260904L,
  random_iterations = 10000L,
  random_sets_to_test = c(
    "IR_FULL_771"
  )
)

classify_ir_pathway <- function(pathway_name) {
  broad_metabolic <- grepl(
    paste(c(
      "TYPE_II_DIABETES", "PPAR", "GLUCOSE_HOMEOSTASIS",
      "REGULATION_OF_GLUCOSE_METABOLIC"
    ), collapse = "|"),
    pathway_name,
    ignore.case = TRUE
  )
  downstream_contextual <- grepl(
    paste(c(
      "MTOR", "MTORC1", "FOXO", "AUTOPHAGY", "GATOR",
      "FLCN", "TSC1_2", "AMINO_ACIDS_REGULATE_MTORC1"
    ), collapse = "|"),
    pathway_name,
    ignore.case = TRUE
  )
  core_insulin <- grepl(
    paste(c(
      "INSULIN_SIGNALING", "INSULIN_RECEPTOR",
      "SIGNALLING_BY_INSULIN_RECEPTOR", "SIGNALING_BY_INSULIN_RECEPTOR",
      "CELLULAR_RESPONSE_TO_INSULIN", "RESPONSE_TO_INSULIN",
      "GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN", "SLC2A4", "GLUT4"
    ), collapse = "|"),
    pathway_name,
    ignore.case = TRUE
  )

  data.table::fifelse(
    pathway_name %in% c("IR_FULL_771", "IR_CANONICAL_UNION"),
    "aggregate_set",
    data.table::fifelse(
      broad_metabolic,
      "broad_metabolic_or_disease",
      data.table::fifelse(
        downstream_contextual,
        "downstream_or_contextual",
        data.table::fifelse(core_insulin, "core_insulin_response", "other")
      )
    )
  )
}

fit_residual_set_enrichment <- function(gene_table, gene_ids, set_name) {
  gene_ids <- unique(as.character(gene_ids))
  working <- data.table::copy(gene_table)
  working[, is_set := GENE %in% gene_ids]
  working <- working[is.finite(zresid) & !is.na(is_set)]
  n_set <- sum(working$is_set)
  n_background <- sum(!working$is_set)
  if (n_set < 2L || n_background < 2L) {
    return(data.table::data.table(
      set_name = set_name,
      status = "not_run_insufficient_genes",
      n_entrez_input = length(gene_ids),
      n_magma_eligible = n_set,
      n_background = n_background,
      beta = NA_real_,
      se = NA_real_,
      t = NA_real_,
      p_positive = NA_real_,
      p_two_sided = NA_real_
    ))
  }

  beta_background <- mean(working[is_set == FALSE, zresid])
  beta <- mean(working[is_set == TRUE, zresid]) - beta_background
  fitted <- beta_background + beta * as.numeric(working$is_set)
  residual <- working$zresid - fitted
  sigma2 <- sum(residual^2) / (nrow(working) - 2L)
  se <- sqrt(sigma2 * (1 / n_set + 1 / n_background))
  t_value <- beta / se
  data.table::data.table(
    set_name = set_name,
    status = "complete",
    n_entrez_input = length(gene_ids),
    n_magma_eligible = n_set,
    n_background = n_background,
    beta = beta,
    se = se,
    t = t_value,
    p_positive = stats::pt(t_value, df = nrow(working) - 2L, lower.tail = FALSE),
    p_two_sided = 2 * stats::pt(abs(t_value), df = nrow(working) - 2L, lower.tail = FALSE)
  )
}

make_quantile_bin <- function(x, n_bins = 5L) {
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)
  if (sum(ok) == 0L) return(out)
  ranks <- data.table::frank(x[ok], ties.method = "average")
  out[ok] <- pmin(n_bins, pmax(1L, ceiling(ranks / length(ranks) * n_bins)))
  out
}

sample_matched_gene_indices <- function(
    target_table,
    pool_table,
    primary_pools = NULL,
    secondary_pools = NULL,
    strata = NULL) {
  n_target <- nrow(target_table)
  selected <- integer()
  if (is.null(primary_pools)) {
    primary_pools <- split(pool_table$gene_index, pool_table$primary_stratum)
  }
  if (is.null(secondary_pools)) {
    secondary_pools <- split(pool_table$gene_index, pool_table$secondary_stratum)
  }
  if (is.null(strata)) {
    strata <- target_table[
      ,
      .N,
      by = .(primary_stratum, secondary_stratum)
    ]
  }

  for (i in seq_len(nrow(strata))) {
    primary <- strata$primary_stratum[i]
    secondary <- strata$secondary_stratum[i]
    needed <- strata$N[i]
    pool <- primary_pools[[primary]]
    pool <- setdiff(pool, selected)
    if (length(pool) < needed) {
      pool <- secondary_pools[[secondary]]
      pool <- setdiff(pool, selected)
    }
    if (length(pool) < needed) {
      pool <- setdiff(pool_table$GENE, selected)
    }
    chosen <- sample(pool, needed, replace = length(pool) < needed)
    selected <- c(selected, chosen)
  }

  selected <- unique(selected)
  if (length(selected) < n_target) {
    fill_pool <- setdiff(pool_table$gene_index, selected)
    selected <- c(
      selected,
      sample(fill_pool, n_target - length(selected), replace = FALSE)
    )
  }
  selected[seq_len(n_target)]
}

run_matched_random_sets <- function(
    gene_table,
    set_table,
    target_set_name,
    observed,
    n_iterations) {
  target_ids <- unique(
    as.character(set_table[set_name == target_set_name, GENE])
  )
  target_table <- gene_table[GENE %in% target_ids]
  pool_table <- gene_table[!GENE %in% target_ids]
  if (nrow(target_table) < 2L || nrow(pool_table) < nrow(target_table)) {
    return(list(
      summary = data.table::data.table(
        set_name = target_set_name,
        status = "not_run_insufficient_matched_pool"
      ),
      iterations = data.table::data.table()
    ))
  }

  primary_pools <- split(pool_table$gene_index, pool_table$primary_stratum)
  secondary_pools <- split(pool_table$gene_index, pool_table$secondary_stratum)
  strata <- target_table[
    ,
    .N,
    by = .(primary_stratum, secondary_stratum)
  ]
  y <- gene_table$zresid
  y_sum <- sum(y)
  y_sumsq <- sum(y^2)
  fit_index_set <- function(set_index) {
    set_index <- unique(set_index[!is.na(set_index)])
    n_set <- length(set_index)
    n_background <- length(y) - n_set
    if (n_set < 2L || n_background < 2L) {
      return(c(beta = NA_real_, p_positive = NA_real_))
    }
    set_sum <- sum(y[set_index])
    set_sumsq <- sum(y[set_index]^2)
    mean_set <- set_sum / n_set
    beta_background <- (y_sum - set_sum) / n_background
    beta <- mean_set - beta_background
    sse <- (set_sumsq - n_set * mean_set^2) +
      ((y_sumsq - set_sumsq) - n_background * beta_background^2)
    sigma2 <- sse / (length(y) - 2L)
    se <- sqrt(sigma2 * (1 / n_set + 1 / n_background))
    t_value <- beta / se
    c(
      beta = beta,
      p_positive = stats::pt(
        t_value,
        df = length(y) - 2L,
        lower.tail = FALSE
      )
    )
  }

  random_beta <- numeric(n_iterations)
  random_p <- numeric(n_iterations)
  for (iteration in seq_len(n_iterations)) {
    sampled_index <- sample_matched_gene_indices(
      target_table,
      pool_table,
      primary_pools = primary_pools,
      secondary_pools = secondary_pools,
      strata = strata
    )
    fit <- fit_index_set(sampled_index)
    random_beta[iteration] <- fit[["beta"]]
    random_p[iteration] <- fit[["p_positive"]]
  }

  observed_row <- observed[set_name == target_set_name]
  random_iterations <- data.table::data.table(
    set_name = target_set_name,
    iteration = seq_len(n_iterations),
    random_beta = random_beta,
    random_p_positive = random_p
  )
  random_summary <- data.table::data.table(
    set_name = target_set_name,
    status = "complete",
    n_iterations = n_iterations,
    observed_beta = observed_row$beta,
    observed_p_positive = observed_row$p_positive,
    random_beta_mean = mean(random_beta, na.rm = TRUE),
    random_beta_sd = stats::sd(random_beta, na.rm = TRUE),
    random_beta_q025 = stats::quantile(random_beta, 0.025, na.rm = TRUE),
    random_beta_q500 = stats::quantile(random_beta, 0.500, na.rm = TRUE),
    random_beta_q975 = stats::quantile(random_beta, 0.975, na.rm = TRUE),
    empirical_p_positive =
      (1L + sum(random_beta >= observed_row$beta, na.rm = TRUE)) /
      (1L + n_iterations),
    n_target_magma_eligible = nrow(target_table),
    matching_features = paste(
      "chromosome, NSNP quintile, and gene-length quintile;",
      "fallback to NSNP/gene-length quintile or genome background if needed"
    )
  )

  list(summary = random_summary, iterations = random_iterations)
}

run_magma_candidate_set_specificity <- function(
    magma_dir,
    output_dir,
    config) {
  magma_exe <- file.path(project_dir, "magma.exe")
  summary_file <- file.path(magma_dir, "MAGMA_competitive_results.tsv")
  set_file <- file.path(magma_dir, "IR_gene_sets_entrez.txt")
  gene_residual_file <- file.path(
    magma_dir,
    paste0(config$target_trait, "_IR_competitive.gsa.genes.out.txt")
  )
  candidate_archive_file <- file.path(
    project_dir,
    "\u8001\u6587\u4ef6",
    "results_refactored",
    "01_candidate_pathway_membership.csv"
  )
  candidate_archive_display <- paste0(
    "./\u8001\u6587\u4ef6/",
    "results_refactored/01_candidate_pathway_membership.csv"
  )

  input_qc <- data.table::data.table(
    item = c(
      "MAGMA competitive summary",
      "MAGMA set annotation",
      paste0(config$target_trait, " MAGMA gene-residual output"),
      "archived candidate pathway membership"
    ),
    path = c(summary_file, set_file, gene_residual_file, candidate_archive_display),
    check_path = c(summary_file, set_file, gene_residual_file, candidate_archive_file)
  )
  input_qc[, present := file.exists(check_path)]
  input_qc[, last_write_time := as.character(file.info(check_path)$mtime)]
  input_qc[, status := ifelse(present, "available", "missing")]
  input_qc[, check_path := NULL]

  if (!file.exists(summary_file) ||
      !file.exists(set_file) ||
      !file.exists(gene_residual_file)) {
    empty <- data.table::data.table(status = "missing_required_input")
    return(list(
      readme = empty,
      input_qc = input_qc,
      pathway_rules = empty,
      pathway_annotation = empty,
      set_sensitivity = empty,
      leave_one_pathway = empty,
      matched_random_summary = empty,
      matched_random_iterations = empty,
      existing_magma_summary = empty
    ))
  }

  magma_summary <- data.table::fread(summary_file)
  magma_sets <- data.table::fread(set_file, header = FALSE)
  data.table::setnames(magma_sets, c("set_name", "GENE"))
  magma_sets[, `:=`(
    set_name = as.character(set_name),
    GENE = as.character(GENE)
  )]

  gene_table <- data.table::fread(gene_residual_file, skip = "GENE")
  residual_columns <- c("GENE", "CHR", "START", "STOP", "NSNPS", "ZRESID_BASE")
  missing_residual_columns <- setdiff(residual_columns, names(gene_table))
  if (length(missing_residual_columns) > 0L) {
    stop(
      "MAGMA gene-residual file is missing required column(s): ",
      paste(missing_residual_columns, collapse = ", ")
    )
  }
  gene_table <- gene_table[
    ,
    .(
      GENE = as.character(GENE),
      CHR = as.character(CHR),
      START = as.numeric(START),
      STOP = as.numeric(STOP),
      NSNPS = as.numeric(NSNPS),
      zresid = as.numeric(ZRESID_BASE)
    )
  ][is.finite(zresid)]
  gene_table[, gene_length := pmax(1, STOP - START + 1)]
  gene_table[, nsnp_bin := make_quantile_bin(log1p(NSNPS), 5L)]
  gene_table[, length_bin := make_quantile_bin(log1p(gene_length), 5L)]
  gene_table[, primary_stratum := paste(CHR, nsnp_bin, length_bin, sep = "|")]
  gene_table[, secondary_stratum := paste(nsnp_bin, length_bin, sep = "|")]
  gene_table[, gene_index := .I]

  aggregate_sets <- c("IR_FULL_771", "IR_CANONICAL_UNION")
  pathway_names <- setdiff(sort(unique(magma_sets$set_name)), aggregate_sets)
  pathway_annotation <- data.table::data.table(set_name = pathway_names)
  pathway_annotation[, pathway_role := classify_ir_pathway(set_name)]
  pathway_annotation[, n_entrez_ids := magma_sets[
    set_name == pathway_annotation$set_name[.I],
    uniqueN(GENE)
  ], by = .I]
  pathway_annotation[, n_magma_eligible := magma_sets[
    set_name == pathway_annotation$set_name[.I] & GENE %in% gene_table$GENE,
    uniqueN(GENE)
  ], by = .I]
  pathway_annotation[, reviewer_flag := pathway_role %in% c(
    "broad_metabolic_or_disease",
    "downstream_or_contextual"
  )]
  if ("I" %in% names(pathway_annotation)) {
    pathway_annotation[, I := NULL]
  }

  full_genes <- unique(magma_sets[set_name == "IR_FULL_771", GENE])
  canonical_genes <- unique(magma_sets[set_name == "IR_CANONICAL_UNION", GENE])
  core_pathways <- pathway_annotation[
    pathway_role == "core_insulin_response",
    set_name
  ]
  broad_metabolic_pathways <- pathway_annotation[
    pathway_role == "broad_metabolic_or_disease",
    set_name
  ]
  broad_contextual_pathways <- pathway_annotation[
    reviewer_flag == TRUE,
    set_name
  ]
  core_genes <- unique(magma_sets[set_name %in% core_pathways, GENE])
  broad_metabolic_genes <- unique(
    magma_sets[set_name %in% broad_metabolic_pathways, GENE]
  )
  broad_contextual_genes <- unique(
    magma_sets[set_name %in% broad_contextual_pathways, GENE]
  )

  sensitivity_sets <- data.table::rbindlist(list(
    data.table::data.table(
      set_name = "IR_FULL_771",
      sensitivity_family = "primary_full_set",
      definition = "Original full IR-related 771-gene candidate set.",
      source_pathways = "all prespecified pathways",
      GENE = full_genes
    ),
    data.table::data.table(
      set_name = "IR_CANONICAL_UNION",
      sensitivity_family = "core_pathway_set",
      definition = "Union of genes from pathways originally labeled canonical.",
      source_pathways = "IR_CANONICAL_UNION",
      GENE = canonical_genes
    ),
    data.table::data.table(
      set_name = "IR_CORE_INSULIN_RESPONSE",
      sensitivity_family = "core_pathway_set",
      definition = paste0(
        "Union of insulin, insulin-receptor, insulin-response, and ",
        "GLUT4/SLC2A4 response pathways."
      ),
      source_pathways = paste(core_pathways, collapse = "; "),
      GENE = core_genes
    ),
    data.table::data.table(
      set_name = "IR_EXCLUDE_BROAD_METABOLIC",
      sensitivity_family = "broad_pathway_exclusion",
      definition = paste0(
        "Full set after conservatively removing any gene annotated to ",
        "type 2 diabetes, PPAR, glucose homeostasis, or regulation of ",
        "glucose metabolic process pathways."
      ),
      source_pathways = paste(broad_metabolic_pathways, collapse = "; "),
      GENE = setdiff(full_genes, broad_metabolic_genes)
    ),
    data.table::data.table(
      set_name = "IR_EXCLUDE_BROAD_METABOLIC_DOWNSTREAM",
      sensitivity_family = "broad_pathway_exclusion",
      definition = paste0(
        "Full set after conservatively removing any gene annotated to ",
        "broad metabolic, disease, FOXO, mTOR/mTORC1, or related ",
        "downstream contextual pathways flagged by the reviewer concern."
      ),
      source_pathways = paste(broad_contextual_pathways, collapse = "; "),
      GENE = setdiff(full_genes, broad_contextual_genes)
    )
  ), use.names = TRUE, fill = TRUE)

  set_metadata <- sensitivity_sets[
    ,
    .(
      sensitivity_family = unique(sensitivity_family)[1L],
      definition = unique(definition)[1L],
      source_pathways = unique(source_pathways)[1L],
      n_entrez_ids = uniqueN(GENE)
    ),
    by = set_name
  ]

  ## Reviewer 1 deep revision: run the alternative candidate sets through
  ## MAGMA so the sensitivity tests retain gene-gene correlation and the same
  ## technical-property corrections as the primary competitive analysis.
  leave_one_membership <- data.table::rbindlist(lapply(
    pathway_names,
    function(pathway) {
      removed_genes <- unique(magma_sets[set_name == pathway, GENE])
      remaining_genes <- setdiff(full_genes, removed_genes)
      data.table::data.table(
        set_name = paste0("IR_FULL_771_minus_", pathway),
        GENE = remaining_genes,
        removed_pathway = pathway,
        removed_pathway_role = classify_ir_pathway(pathway),
        n_removed_entrez_ids = length(removed_genes),
        n_remaining_entrez_ids = length(remaining_genes)
      )
    }
  ), use.names = TRUE, fill = TRUE)

  magma_sensitivity_membership <- data.table::rbindlist(list(
    sensitivity_sets[, .(set_name, GENE)],
    leave_one_membership[, .(set_name, GENE)]
  ), use.names = TRUE)
  magma_sensitivity_membership <- unique(
    magma_sensitivity_membership[!is.na(GENE)],
    by = c("set_name", "GENE")
  )

  sensitivity_tmp <- tempfile("candidate_magma_")
  dir.create(sensitivity_tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(sensitivity_tmp, recursive = TRUE, force = TRUE), add = TRUE)
  sensitivity_annot <- file.path(sensitivity_tmp, "candidate_sets.txt")
  sensitivity_prefix <- file.path(sensitivity_tmp, "ALM_candidate_sets")
  data.table::fwrite(
    magma_sensitivity_membership[, .(set_name, GENE)],
    sensitivity_annot,
    sep = " ", col.names = FALSE
  )
  candidate_magma_path <- function(path, mustWork = TRUE) {
    if (mustWork) {
      return(utils::shortPathName(normalizePath(
        path, winslash = "\\", mustWork = TRUE
      )))
    }
    parent <- utils::shortPathName(normalizePath(
      dirname(path), winslash = "\\", mustWork = TRUE
    ))
    file.path(parent, basename(path))
  }
  candidate_magma_log <- system2(
    candidate_magma_path(magma_exe),
    args = c(
      "--gene-results",
      candidate_magma_path(file.path(magma_dir, "ALM_gene.genes.raw")),
      "--set-annot", candidate_magma_path(sensitivity_annot), "col=2,1",
      "--out", candidate_magma_path(sensitivity_prefix, mustWork = FALSE)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  candidate_magma_status <- attr(candidate_magma_log, "status")
  if (is.null(candidate_magma_status)) candidate_magma_status <- 0L
  if (candidate_magma_status != 0L) {
    stop(
      "MAGMA candidate-set sensitivity analysis failed:\n",
      paste(tail(candidate_magma_log, 20L), collapse = "\n")
    )
  }
  magma_sensitivity <- data.table::fread(
    paste0(sensitivity_prefix, ".gsa.out.txt"),
    skip = "VARIABLE", header = TRUE
  )
  if (!"FULL_NAME" %in% names(magma_sensitivity)) {
    magma_sensitivity[, FULL_NAME := NA_character_]
  }
  magma_sensitivity[, set_name := data.table::fifelse(
    !is.na(FULL_NAME), FULL_NAME, VARIABLE
  )]
  magma_sensitivity <- magma_sensitivity[, .(
    set_name,
    status = "complete_MAGMA_competitive",
    n_magma_eligible = NGENES,
    beta = BETA,
    beta_std = BETA_STD,
    se = SE,
    t = BETA / SE,
    p_positive = P
  )]

  set_sensitivity <- merge(
    set_metadata,
    magma_sensitivity[set_name %in% set_metadata$set_name],
    by = "set_name",
    all.x = TRUE
  )
  set_sensitivity[, `:=`(
    n_entrez_input = n_entrez_ids,
    n_background = nrow(gene_table) - n_magma_eligible,
    fdr_positive_BH = stats::p.adjust(p_positive, method = "BH")
  )]
  data.table::setorder(set_sensitivity, sensitivity_family, p_positive)

  # The matched-random empirical null uses residualised gene Z scores. Compute
  # the observed statistic on that same scale rather than using the MAGMA beta.
  random_observed <- data.table::rbindlist(lapply(
    unique(sensitivity_sets$set_name),
    function(current_set) {
      fit_residual_set_enrichment(
        gene_table,
        sensitivity_sets[set_name == current_set, GENE],
        current_set
      )
    }
  ), use.names = TRUE, fill = TRUE)

  leave_one_pathway <- merge(
    unique(leave_one_membership[, .(
      set_name, removed_pathway, removed_pathway_role,
      n_removed_entrez_ids, n_remaining_entrez_ids
    )]),
    magma_sensitivity,
    by = "set_name",
    all.x = TRUE
  )
  leave_one_pathway[, `:=`(
    n_entrez_input = n_remaining_entrez_ids,
    n_background = nrow(gene_table) - n_magma_eligible,
    fdr_positive_BH = stats::p.adjust(p_positive, method = "BH")
  )]
  leave_one_pathway[, retained_nominal_ALM_enrichment :=
    is.finite(p_positive) & p_positive < 0.05 & beta > 0]
  data.table::setcolorder(
    leave_one_pathway,
    c(
      "removed_pathway", "removed_pathway_role",
      "retained_nominal_ALM_enrichment", "n_removed_entrez_ids",
      "n_remaining_entrez_ids", "n_magma_eligible", "beta", "se",
      "beta_std", "t", "p_positive", "fdr_positive_BH", "status"
    )
  )
  data.table::setorder(leave_one_pathway, -p_positive)

  random_targets <- intersect(
    config$random_sets_to_test,
    unique(sensitivity_sets$set_name)
  )
  set.seed(config$random_seed)
  random_results <- lapply(random_targets, function(current_set) {
    run_matched_random_sets(
      gene_table = gene_table,
      set_table = sensitivity_sets,
      target_set_name = current_set,
      observed = random_observed,
      n_iterations = config$random_iterations
    )
  })
  matched_random_summary <- data.table::rbindlist(
    lapply(random_results, `[[`, "summary"),
    use.names = TRUE,
    fill = TRUE
  )
  matched_random_iterations <- data.table::rbindlist(
    lapply(random_results, `[[`, "iterations"),
    use.names = TRUE,
    fill = TRUE
  )

  existing_magma_summary <- magma_summary[
    trait %in% magma_traits &
      (
        VARIABLE %in% aggregate_sets |
          FULL_NAME %in% pathway_names |
          VARIABLE %in% pathway_names
      )
  ]
  existing_magma_summary[, pathway_role := classify_ir_pathway(
    fifelse(!is.na(FULL_NAME), FULL_NAME, VARIABLE)
  )]

  pathway_rules <- data.table::data.table(
    rule = c(
      "candidate_space",
      "canonical_component",
      "expanded_component",
      "broad_metabolic_exclusion",
      "broad_metabolic_downstream_exclusion",
      "leave_one_pathway_out",
      "matched_random_sets",
      "outcome_independence",
      "timestamp_evidence"
    ),
    implementation = c(
      paste0(
        "Uses the frozen MAGMA set annotation in IR_gene_sets_entrez.txt; ",
        "duplicate genes across overlapping pathways are collapsed."
      ),
      "IR_CANONICAL_UNION from the original MAGMA set annotation.",
      paste0(
        "Pathway names retained in IR_gene_sets_entrez.txt after the ",
        "original MSigDB exact-name and keyword-selection procedure."
      ),
      paste(
        "Remove any full-set gene annotated to pathways classified as",
        "type 2 diabetes, PPAR, glucose homeostasis, or regulation of",
        "glucose metabolic process."
      ),
      paste(
        "Remove any full-set gene annotated to broad metabolic/disease",
        "pathways or FOXO, mTOR/mTORC1, autophagy, GATOR, FLCN, TSC1/2,",
        "or amino-acid mTORC1 contextual pathways."
      ),
      "Remove one individual pathway at a time from IR_FULL_771.",
      paste0(
        config$random_iterations,
        " random gene sets matched to each target set by chromosome, ",
        "NSNP quintile, and gene-length quintile."
      ),
      paste(
        "Candidate pathways and membership are constructed and frozen before",
        "the pipeline reads any sarcopenia-related GWAS; outcome, eQTL, MR,",
        "and colocalization results are not candidate-selection inputs."
      ),
      paste0(
        "Archived candidate membership file mtime: ",
        input_qc[item == "archived candidate pathway membership", last_write_time],
        "; MAGMA summary mtime: ",
        input_qc[item == "MAGMA competitive summary", last_write_time]
      )
    )
  )

  readme <- data.table::data.table(
    item = c(
      "analysis_goal",
      "target_trait",
      "random_seed",
      "random_iterations",
      "interpretation_boundary"
    ),
    value = c(
      paste0(
        "Assess whether ALM enrichment of the IR-related candidate set is ",
        "robust to candidate-set definition and pathway overlap."
      ),
      config$target_trait,
      as.character(config$random_seed),
      as.character(config$random_iterations),
      paste0(
        "Alternative-set and leave-one-pathway-out results are MAGMA competitive ",
        "tests. The matched-random empirical null is an additional residual-Z ",
        "analysis on a separate coefficient scale and is not compared directly ",
        "with MAGMA beta coefficients."
      )
    )
  )

  list(
    readme = readme,
    input_qc = input_qc,
    pathway_rules = pathway_rules,
    pathway_annotation = pathway_annotation,
    set_sensitivity = set_sensitivity,
    leave_one_pathway = leave_one_pathway,
    matched_random_summary = matched_random_summary,
    matched_random_iterations = matched_random_iterations,
    existing_magma_summary = existing_magma_summary
  )
}

candidate_set_specificity <- run_magma_candidate_set_specificity(
  magma_dir = magma_dir,
  output_dir = output_dir,
  config = candidate_specificity_config
)

candidate_specificity_workbook_file <- file.path(
  output_dir,
  "MAGMA_candidate_set_specificity_sensitivity.xlsx"
)

repair_ooxml_xlsx <- function(xlsx_file) {
  if (!requireNamespace("xml2", quietly = TRUE) ||
      !requireNamespace("zip", quietly = TRUE)) {
    stop("Workbook compatibility repair requires xml2 and zip.")
  }
  if (!file.exists(xlsx_file)) {
    stop("Workbook not found: ", xlsx_file)
  }

  col_to_num <- function(column_label) {
    letters <- strsplit(column_label, "", fixed = TRUE)[[1L]]
    sum((match(letters, LETTERS)) * 26 ^ rev(seq_along(letters) - 1L))
  }
  num_to_col <- function(column_number) {
    out <- character()
    while (column_number > 0L) {
      remainder <- (column_number - 1L) %% 26L
      out <- c(LETTERS[remainder + 1L], out)
      column_number <- (column_number - remainder - 1L) %/% 26L
    }
    paste(out, collapse = "")
  }
  normalise_ooxml_path <- function(path) {
    path <- gsub("\\\\", "/", path)
    path_parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
    resolved_parts <- character()
    for (path_part in path_parts) {
      if (!nzchar(path_part) || path_part == ".") next
      if (path_part == "..") {
        if (length(resolved_parts) > 0L) {
          resolved_parts <- resolved_parts[-length(resolved_parts)]
        }
      } else {
        resolved_parts <- c(resolved_parts, path_part)
      }
    }
    paste(resolved_parts, collapse = "/")
  }

  repair_dir <- tempfile("candidate_specificity_xlsx_repair_")
  repaired_file <- tempfile(fileext = ".xlsx")
  dir.create(repair_dir, recursive = TRUE)
  on.exit(unlink(repair_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(repaired_file, force = TRUE), add = TRUE)

  utils::unzip(xlsx_file, exdir = repair_dir)
  archive_files <- list.files(
    repair_dir,
    recursive = TRUE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  archive_files <- gsub("\\\\", "/", archive_files)

  relationship_files <- list.files(
    repair_dir,
    pattern = "\\.rels$",
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  removed_relationships <- 0L
  for (relationship_file in relationship_files) {
    relationship_path <- substring(
      gsub("\\\\", "/", relationship_file),
      nchar(gsub("\\\\", "/", repair_dir)) + 2L
    )
    relationship_base <- if (relationship_path == "_rels/.rels") {
      ""
    } else {
      sub("/_rels/[^/]+\\.rels$", "", relationship_path)
    }
    relationship_xml <- xml2::read_xml(relationship_file)
    relationship_nodes <- xml2::xml_find_all(
      relationship_xml,
      "//*[local-name()='Relationship']"
    )
    for (relationship_node in relationship_nodes) {
      if (identical(xml2::xml_attr(relationship_node, "TargetMode"), "External")) {
        next
      }
      relationship_target <- xml2::xml_attr(relationship_node, "Target")
      target_path <- if (startsWith(relationship_target, "/")) {
        sub("^/+", "", relationship_target)
      } else {
        normalise_ooxml_path(file.path(
          relationship_base,
          relationship_target
        ))
      }
      if (!target_path %in% archive_files) {
        xml2::xml_remove(relationship_node)
        removed_relationships <- removed_relationships + 1L
      }
    }
    xml2::write_xml(relationship_xml, relationship_file)
  }

  content_types_file <- file.path(repair_dir, "[Content_Types].xml")
  content_types_xml <- xml2::read_xml(content_types_file)
  override_nodes <- xml2::xml_find_all(
    content_types_xml,
    "//*[local-name()='Override']"
  )
  removed_content_types <- 0L
  for (override_node in override_nodes) {
    part_name <- sub("^/+", "", xml2::xml_attr(override_node, "PartName"))
    if (!part_name %in% archive_files) {
      xml2::xml_remove(override_node)
      removed_content_types <- removed_content_types + 1L
    }
  }
  xml2::write_xml(content_types_xml, content_types_file)

  sheet_files <- list.files(
    file.path(repair_dir, "xl", "worksheets"),
    pattern = "^sheet[0-9]+\\.xml$",
    full.names = TRUE
  )
  repaired_dimensions <- 0L
  for (sheet_file in sheet_files) {
    sheet_xml <- xml2::read_xml(sheet_file)
    cell_refs <- xml2::xml_attr(
      xml2::xml_find_all(sheet_xml, "//*[local-name()='c']"),
      "r"
    )
    cell_refs <- cell_refs[!is.na(cell_refs)]
    if (length(cell_refs) == 0L) next
    cell_columns <- gsub("[0-9]", "", cell_refs)
    cell_rows <- suppressWarnings(as.integer(gsub("[A-Z]", "", cell_refs)))
    max_col <- max(vapply(cell_columns, col_to_num, numeric(1)), na.rm = TRUE)
    max_row <- max(cell_rows, na.rm = TRUE)
    dimension_node <- xml2::xml_find_first(
      sheet_xml,
      "//*[local-name()='dimension']"
    )
    if (!is.na(xml2::xml_name(dimension_node))) {
      xml2::xml_set_attr(
        dimension_node,
        "ref",
        paste0("A1:", num_to_col(max_col), max_row)
      )
      repaired_dimensions <- repaired_dimensions + 1L
      xml2::write_xml(sheet_xml, sheet_file)
    }
  }

  top_level_entries <- list.files(
    repair_dir,
    all.files = TRUE,
    no.. = TRUE,
    full.names = FALSE
  )
  previous_directory <- setwd(repair_dir)
  tryCatch(
    zip::zipr(
      zipfile = repaired_file,
      files = top_level_entries,
      recurse = TRUE,
      include_directories = FALSE
    ),
    finally = setwd(previous_directory)
  )
  if (!file.exists(repaired_file) || file.info(repaired_file)$size == 0L) {
    stop("Failed to create the repaired workbook.")
  }
  if (!file.copy(repaired_file, xlsx_file, overwrite = TRUE)) {
    stop("Failed to replace the candidate-set specificity workbook.")
  }

  data.table::data.table(
    removed_relationships = removed_relationships,
    removed_content_types = removed_content_types,
    repaired_dimensions = repaired_dimensions
  )
}

candidate_specificity_tables <- list(
  `00_README` = as.data.frame(candidate_set_specificity$readme),
  `01_input_QC` = as.data.frame(candidate_set_specificity$input_qc),
  `02_pathway_rules` = as.data.frame(candidate_set_specificity$pathway_rules),
  `03_pathway_annotation` = as.data.frame(candidate_set_specificity$pathway_annotation),
  `04_set_sensitivity_ALM` = as.data.frame(candidate_set_specificity$set_sensitivity),
  `05_leave_one_pathway_ALM` = as.data.frame(candidate_set_specificity$leave_one_pathway),
  `06_matched_random_summary` = as.data.frame(candidate_set_specificity$matched_random_summary),
  `07_matched_random_iterations` = as.data.frame(candidate_set_specificity$matched_random_iterations),
  `08_existing_MAGMA_summary` = as.data.frame(candidate_set_specificity$existing_magma_summary)
)
candidate_specificity_workbook <- openxlsx::createWorkbook()
candidate_specificity_header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
for (sheet_name in names(candidate_specificity_tables)) {
  sheet_data <- candidate_specificity_tables[[sheet_name]]
  openxlsx::addWorksheet(
    candidate_specificity_workbook,
    sheet_name,
    gridLines = FALSE
  )
  openxlsx::writeData(candidate_specificity_workbook, sheet_name, sheet_data)
  if (ncol(sheet_data) > 0L) {
    openxlsx::addStyle(
      candidate_specificity_workbook,
      sheet_name,
      candidate_specificity_header_style,
      rows = 1L,
      cols = seq_len(ncol(sheet_data)),
      gridExpand = TRUE
    )
    openxlsx::freezePane(
      candidate_specificity_workbook,
      sheet_name,
      firstRow = TRUE
    )
    widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- head(as.character(sheet_data[[column_index]]), 500L)
      observed_lengths <- nchar(values)
      observed <- if (
        length(observed_lengths) == 0L || all(is.na(observed_lengths))
      ) 0L else max(observed_lengths, na.rm = TRUE)
      min(55, max(12, nchar(names(sheet_data)[column_index]) + 2L, observed + 2L))
    }, numeric(1))
    openxlsx::setColWidths(
      candidate_specificity_workbook,
      sheet_name,
      cols = seq_len(ncol(sheet_data)),
      widths = widths
    )
  }
}
openxlsx::saveWorkbook(
  candidate_specificity_workbook,
  candidate_specificity_workbook_file,
  overwrite = TRUE
)
candidate_specificity_workbook_repair <- repair_ooxml_xlsx(
  candidate_specificity_workbook_file
)
cross_trait_workbook_repair <- repair_ooxml_xlsx(
  cross_trait_workbook_file
)
if (!file.exists(candidate_specificity_workbook_file)) {
  stop("Candidate-set specificity workbook validation failed: file was not created.")
}

## 07c. END NEW CANDIDATE-SET SPECIFICITY SENSITIVITY MODULE -----------------

supplementary_tables <- list(
  `S1_IR_Pathways` = read_sheet("main", "01_candidate_pathway_membership"),
  `S2_IR_Genes` = read_sheet("main", "01_candidate_gene_layers"),
  `S3_MAGMA` = utils::read.delim(
    file.path(magma_dir, "MAGMA_competitive_results.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  `S4_All_main_cisMR` = filter_main_traits(read_sheet("main", "02_cisMR_all_results")),
  `S5_FDR_main_cisMR` = filter_main_traits(read_sheet("main", "02_cisMR_FDR_significant")),
  `S6_MR_Sensitivity` = filter_main_traits(read_sheet("main", "03_cisMR_sensitivity_multiSNP")),
  `S7_Coloc` = filter_main_traits(read_sheet("main", "04_coloc_all_priors")),
  `S8_FUSION` = fusion,
  `S9_SuSiE` = susie,
  `S10_snRNA_Localization` = sn_localization,
  `S11_snRNA_DE` = sn_de,
  `S12_sQTL` = sqtl_mr,
  `S13_MAGIC_main_outcomes` = magic_two_step_main,
  `S14_GEFOS_ALM` = read_sheet("reviewer7", "07_external_evidence_label"),
  `S15_MAGMA_CrossTrait_Beta` = magma_cross_trait$beta_summary,
  `S16_MAGMA_CrossTrait_Primary` = magma_cross_trait$primary_loco,
  `S17_MAGMA_CrossTrait_Residual` = magma_cross_trait$residual_sensitivity,
  `S18_MAGMA_CrossTrait_QC` = magma_cross_trait$input_qc
)


## 08. Figure and table catalog ---------------------------------------------

catalog <- data.frame(
  Item = c(
    "Figure 1", "Figure 2", "Figure 3", "Figure 4",
    "Table 1", "Table 2",
    paste("Supplementary Table", paste0("S", 1:18))
  ),
  Chinese_name = c(
    "研究设计与证据整合流程",
    "6个Tier A基因的骨骼肌cis-MR森林图",
    "Tier A基因的细胞定位、年龄相关表达与剪接层证据",
    "6个Tier A基因的多层证据矩阵",
    "研究数据来源及用途",
    "6个Tier A候选基因的MR及共定位结果",
    "24条预设通路及成员关系",
    "771个候选基因及通路层级",
    "MAGMA竞争性富集结果",
    "三个主要结局的全部551项cis-MR结果",
    "三个主要结局FDR显著的84项cis-MR结果",
    "MR方向性与敏感性分析",
    "共定位及先验敏感性结果",
    "FUSION骨骼肌eQTL来源评估",
    "SuSiE多信号共定位结果",
    "单核RNA细胞定位结果",
    "单核RNA伪批量差异表达结果",
    "GTEx v11骨骼肌sQTL-MR结果",
    "MAGIC IR机制边界分析",
    "GEFOS 2017 ALM外部验证",
    "MAGMA主富集系数汇总",
    "MAGMA跨表型共同基因空间留一染色体主检验",
    "MAGMA跨表型残差尺度敏感性分析",
    "MAGMA跨表型比较输入QC"
  ),
  English_name = c(
    "Study design and evidence-integration workflow",
    "Skeletal-muscle cis-MR estimates for Tier A genes",
    "Cellular and splicing-layer context of Tier A genes",
    "Multi-layer evidence profile of the six Tier A genes",
    "Data sources and analytical roles",
    "MR and colocalization results for the six Tier A genes",
    "Prespecified pathways and membership",
    "The 771 candidate genes and pathway layers",
    "MAGMA competitive gene-set results",
    "All main-outcome skeletal-muscle cis-MR results",
    "FDR-significant main-outcome skeletal-muscle cis-MR results",
    "MR directionality and sensitivity analyses",
    "Colocalization and prior-sensitivity results",
    "FUSION muscle eQTL-source assessment",
    "SuSiE multi-signal colocalization",
    "Single-nucleus cellular localization",
    "Single-nucleus pseudobulk differential expression",
    "GTEx v11 skeletal-muscle sQTL-MR results",
    "Exploratory MAGIC IR mechanism-boundary analysis",
    "GEFOS 2017 ALM external validation",
    "Primary MAGMA enrichment coefficients",
    "Primary common-universe leave-one-chromosome MAGMA contrasts",
    "Residual-scale sensitivity analyses for cross-phenotype contrasts",
    "Input QC for cross-phenotype MAGMA comparisons"
  ),
  Output_file_or_sheet = c(
    "Figure1_Study_workflow.png; Figure1_Study_workflow.pdf",
    "Figure2_TierA_cisMR_forest.png; Figure2_TierA_cisMR_forest.pdf",
    "Figure3_Cell_and_splicing_context.png; Figure3_Cell_and_splicing_context.pdf",
    "Figure4_Multilayer_evidence_matrix.png; Figure4_Multilayer_evidence_matrix.pdf",
    "Table1_DataSources", "Table2_TierA",
    names(supplementary_tables)
  ),
  Manuscript_location = c(
    rep("Main text", 6), rep("Supplementary material", 18)
  ),
  stringsAsFactors = FALSE
)


## 09. Write the single manuscript table workbook ---------------------------

output_tables <- c(
  list(
    `00_Catalog` = catalog,
    `Table1_DataSources` = table1_sources,
    `Table2_TierA` = as.data.frame(table2_tier_a)
  ),
  lapply(supplementary_tables, as.data.frame)
)

workbook_file <- file.path(output_dir, "Manuscript_tables_and_catalog.xlsx")
workbook <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
note_style <- openxlsx::createStyle(
  fgFill = "#DCEAF7", valign = "top", wrapText = TRUE
)

for (sheet_name in names(output_tables)) {
  sheet_data <- output_tables[[sheet_name]]
  openxlsx::addWorksheet(workbook, sheet_name, gridLines = FALSE)
  openxlsx::writeDataTable(
    workbook, sheet_name, sheet_data, tableStyle = "TableStyleMedium2"
  )
  if (ncol(sheet_data) > 0L) {
    openxlsx::addStyle(
      workbook, sheet_name, header_style,
      rows = 1L, cols = seq_len(ncol(sheet_data)), gridExpand = TRUE
    )
    openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
    widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- head(as.character(sheet_data[[column_index]]), 500L)
      observed_lengths <- nchar(values)
      observed <- if (
        length(observed_lengths) == 0L || all(is.na(observed_lengths))
      ) 0L else max(observed_lengths, na.rm = TRUE)
      min(45, max(11, nchar(names(sheet_data)[column_index]) + 2L, observed + 2L))
    }, numeric(1))
    openxlsx::setColWidths(
      workbook, sheet_name, cols = seq_len(ncol(sheet_data)), widths = widths
    )
  }
  if (sheet_name == "00_Catalog") {
    openxlsx::addStyle(
      workbook, sheet_name, note_style,
      rows = 2:(nrow(sheet_data) + 1L), cols = 2:3,
      gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setColWidths(
      workbook, sheet_name, cols = 1:5,
      widths = c(24, 48, 58, 62, 24)
    )
  }
}
openxlsx::saveWorkbook(workbook, workbook_file, overwrite = TRUE)
manuscript_workbook_repair <- repair_ooxml_xlsx(workbook_file)


## 10. Final validation -------------------------------------------------------

expected_files <- c(
  unlist(lapply(c(
    "Figure1_Study_workflow", "Figure2_TierA_cisMR_forest",
    "Figure3_Cell_and_splicing_context", "Figure4_Multilayer_evidence_matrix"
  ), function(stem) paste0(stem, c(".png", ".pdf")))),
  basename(workbook_file)
)
missing_outputs <- expected_files[
  !file.exists(file.path(output_dir, expected_files))
]
if (length(missing_outputs) > 0L) {
  stop("Manuscript output validation failed: ", paste(missing_outputs, collapse = ", "))
}
if (!identical(openxlsx::getSheetNames(workbook_file), names(output_tables))) {
  stop("Manuscript workbook validation failed: sheet names differ.")
}
if (!file.exists(cross_trait_workbook_file)) {
  stop("Cross-phenotype MAGMA workbook validation failed: file was not created.")
}
if (nrow(table2_tier_a) != 6L || nrow(table3_multilayer) != 6L) {
  stop("Manuscript table validation failed: Tier A candidate rows missing.")
}

cat(
  "Manuscript figures and tables created in:\n",
  output_dir, "\n",
  "Cross-phenotype MAGMA workbook:\n",
  cross_trait_workbook_file, "\n",
  paste(expected_files, collapse = "\n"), "\n"
)
