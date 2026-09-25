# [NEW: final submission outputs]
# Regenerates the 13 main supplementary tables, 13 supplementary figures,
# and copies the four main-text figures from one consistent analysis release.

if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# [UPDATED 2026-09-17] Resolve this package from the script, never from cwd.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
if (is.na(script_file) || basename(script_file) != "04_submission_outputs.R") {
  script_file <- NA_character_
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (is.character(candidate) && length(candidate) == 1L && !is.na(candidate) &&
        nzchar(candidate) && basename(candidate) == "04_submission_outputs.R") {
      script_file <- candidate
      break
    }
  }
}
if (is.na(script_file) || !file.exists(script_file)) stop("Cannot resolve submission script path.")
package_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
submission_dir <- file.path(package_dir, "submission_outputs")
supp_figure_dir <- file.path(submission_dir, "Supplementary_figures")
main_figure_dir <- file.path(submission_dir, "Main_figures")
dir.create(supp_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Standalone rendering uses only this package's frozen processed snapshot.
# A full workflow passes its current project directory explicitly.
live_root <- Sys.getenv("IR_SUBMISSION_SOURCE_DIR", unset = "")
search_roots <- if (nzchar(live_root)) {
  live_root <- normalizePath(live_root, winslash = "/", mustWork = TRUE)
  c(file.path(live_root, "results_public_release"),
    file.path(live_root, "manuscript_outputs_public_release"))
} else {
  c(file.path(package_dir, "processed_results"), file.path(package_dir, "figures"))
}
if (any(!dir.exists(search_roots))) stop("Missing submission source directory: ",
                                        paste(search_roots[!dir.exists(search_roots)], collapse = ", "))

find_input <- function(name, required = TRUE) {
  hits <- unique(unlist(lapply(search_roots, function(root) {
    direct <- file.path(root, name)
    nested <- file.path(root, "plots", name)
    c(direct[file.exists(direct)], nested[file.exists(nested)])
  })))
  if (!length(hits)) {
    if (required) stop("Submission source file not found: ", name)
    return(NA_character_)
  }
  hits[[1]]
}

read_sheet <- function(file, sheet) {
  as.data.table(openxlsx::read.xlsx(file, sheet = sheet, check.names = FALSE))
}

select_columns <- function(x, columns) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) x[, (missing) := NA]
  x[, ..columns]
}

analysis_file <- find_input("manuscript_analysis_tables.xlsx")
cross_file <- find_input("MAGMA_cross_phenotype_comparison.xlsx")
specificity_file <- find_input("MAGMA_candidate_set_specificity_sensitivity.xlsx")
magic_file <- find_input("MAGIC_IR_mechanistic_extension.xlsx")
snrna_file <- find_input("GSE167186_snRNA_celltype_validation.xlsx")
sqtl_file <- find_input("GTEx_v11_sQTL_TierA_analysis.xlsx")

magma <- read_sheet(specificity_file, "08_existing_MAGMA_summary")
magma[, FULL_NAME := fifelse(is.na(FULL_NAME) | !nzchar(FULL_NAME), VARIABLE, FULL_NAME)]
s1 <- select_columns(magma, c(
  "trait", "FULL_NAME", "analysis_family", "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR"
))

s2a <- select_columns(read_sheet(cross_file, "02_primary_MAGMA_LOCO"), c(
  "comparison", "trait_1", "trait_2", "beta_1", "beta_2", "delta_beta",
  "jackknife_covariance", "jackknife_correlation", "se_delta_jackknife", "t_statistic",
  "degrees_of_freedom", "p_value", "p_value_normal_reference", "ci_95_lower",
  "ci_95_upper", "n_common_genes", "n_candidate_genes", "n_blocks", "fdr_BH"
))
s2b <- select_columns(read_sheet(cross_file, "05_residual_sensitivity"), c(
  "scale", "comparison", "trait_1", "trait_2", "residual_contrast_beta",
  "se_delta_jackknife", "t_statistic", "degrees_of_freedom", "p_value",
  "p_value_normal_reference", "ci_95_lower", "ci_95_upper", "n_common_genes",
  "n_candidate_genes", "n_blocks", "fdr_BH"
))

s3a <- select_columns(read_sheet(specificity_file, "04_set_sensitivity_ALM"), c(
  "set_name", "sensitivity_family", "source_pathways", "n_entrez_ids", "n_magma_eligible",
  "beta", "beta_std", "se", "t", "p_positive", "fdr_positive_BH"
))
s3b <- select_columns(read_sheet(specificity_file, "05_leave_one_pathway_ALM"), c(
  "removed_pathway", "removed_pathway_role", "n_removed_entrez_ids", "n_remaining_entrez_ids",
  "n_magma_eligible", "beta", "se", "beta_std", "t", "p_positive", "fdr_positive_BH"
))
s3c <- select_columns(read_sheet(specificity_file, "06_matched_random_summary"), c(
  "set_name", "n_iterations", "observed_beta", "observed_p_positive", "random_beta_mean",
  "random_beta_sd", "random_beta_q025", "random_beta_q500", "random_beta_q975",
  "empirical_p_positive", "n_target_magma_eligible"
))

s4a <- select_columns(read_sheet(specificity_file, "03_pathway_annotation"), c(
  "set_name", "pathway_role", "n_entrez_ids", "n_magma_eligible"
))
setorder(s4a, set_name)
s4b <- select_columns(read_sheet(analysis_file, "01_candidate_pathway_membership"), c(
  "gs_name", "gene_symbol", "layer"
))
setorder(s4b, gs_name, gene_symbol)
s4c <- select_columns(read_sheet(analysis_file, "01_candidate_gene_layers"), c(
  "gene_symbol", "in_canonical", "in_expanded", "gene_layer", "n_pathways"
))
setorder(s4c, gene_symbol)

mr_columns <- c(
  "gene_ensembl", "gene_symbol", "gene_layer", "trait", "method", "nsnp", "beta", "se",
  "pval", "min_F", "mean_F", "fdr_within_trait", "fdr_within_trait_layer", "fdr_global",
  "ci_lower", "ci_upper"
)
s5 <- select_columns(read_sheet(analysis_file, "02_cisMR_all_results"), mr_columns)
s5 <- s5[trait %in% c("ALM", "GRIP", "WALK")]
s6 <- select_columns(read_sheet(analysis_file, "02_cisMR_FDR_significant"), mr_columns)
s6 <- s6[trait %in% c("ALM", "GRIP", "WALK")]

s7 <- select_columns(read_sheet(analysis_file, "04_coloc_all_priors"), c(
  "gene_symbol", "trait", "p1", "p2", "p12", "n_snps", "PP.H0", "PP.H1", "PP.H2",
  "PP.H3", "PP.H4", "gene_ensembl", "gene_layer"
))
s7 <- s7[trait %in% c("ALM", "GRIP", "WALK")]

s8a <- select_columns(read_sheet(analysis_file, "04c_coloc_SuSiE_TierA_QC"), c(
  "gene_symbol", "gene_ensembl", "trait", "chromosome", "region_start", "region_end",
  "coloc_abf_n_snps", "ld_reference_panel", "n_harmonised_snps", "n_ld_reference_snps",
  "n_ld_allele_matched_snps", "n_susie_input_snps", "susie_input_limit",
  "variant_selection_status", "susie_status", "n_eqtl_credible_sets",
  "n_outcome_credible_sets", "n_primary_signal_pairs", "max_primary_PP.H3.abf",
  "max_primary_PP.H4.abf", "top_primary_signal_pair"
))
s8b <- select_columns(read_sheet(analysis_file, "04d_coloc_SuSiE_all_priors"), c(
  "nsnps", "hit1", "hit2", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf",
  "PP.H4.abf", "idx1", "idx2", "signal_pair", "gene_symbol", "gene_ensembl",
  "gene_layer", "trait", "p1", "p2", "p12", "n_ld_snps", "PP.H3", "PP.H4"
))
s8c <- select_columns(read_sheet(analysis_file, "04e_coloc_SuSiE_credible_sets"), c(
  "gene_symbol", "gene_ensembl", "trait", "dataset", "credible_set", "susie_component",
  "cs_size", "lead_snp", "lead_pip", "coverage", "purity_min_abs_corr",
  "purity_mean_abs_corr", "min_pos", "max_pos", "member_snps", "member_snps_truncated"
))
s8d <- select_columns(read_sheet(analysis_file, "04f_coloc_SuSiE_CS_members"), c(
  "gene_symbol", "gene_ensembl", "trait", "dataset", "credible_set", "susie_component",
  "SNP", "pos", "beta", "se", "pval", "pip"
))

s9 <- select_columns(read_sheet(analysis_file, "08_FUSION_independent_replicati"), c(
  "gene_ensembl", "gene_symbol", "trait", "discovery_beta", "discovery_PP.H4",
  "n_fusion_region_snps", "n_fusion_strong_instruments", "n_fusion_clumped_instruments",
  "fusion_mr_nsnp", "fusion_mr_beta", "fusion_mr_se", "fusion_mr_pval", "fusion_PP.H3",
  "fusion_PP.H4", "fusion_mr_fdr"
))

gefos <- read_sheet(analysis_file, "09_GEFOS_external_evidence")
gefos_coloc <- read_sheet(analysis_file, "09_GEFOS_coloc_all_priors")[p12 == 5e-6]
gefos_coloc <- gefos_coloc[, .(gene_symbol, nsnps = n_snps)]
s10 <- merge(gefos, gefos_coloc, by = "gene_symbol", all.x = TRUE, sort = FALSE)
setorder(s10, target_order)
s10 <- select_columns(s10, c(
  "gene_symbol", "gene_ensembl", "discovery_outcome", "external_method", "external_nsnp",
  "external_beta", "external_se", "external_ci_lower", "external_ci_upper", "external_pval",
  "external_p_bonferroni", "external_fdr_4gene", "nominal_directional_support",
  "bonferroni_directional_support", "nsnps", "external_PP3", "external_PP4"
))
setnames(s10, c("external_PP3", "external_PP4"), c("external_coloc_PP3", "external_coloc_PP4"))

s11 <- read_sheet(magic_file, "07_two_step_evidence_chain")
s11 <- s11[sarcopenia_trait %in% c("ALM", "GRIP", "WALK")]
s11 <- select_columns(s11, c(
  "gene_symbol", "ir_trait", "sarcopenia_trait", "gene_ensembl", "beta_gene_ir", "se_gene_ir",
  "p_gene_ir", "fdr_gene_ir", "beta_ir_sarcopenia", "se_ir_sarcopenia",
  "p_ir_sarcopenia", "fdr_ir_sarcopenia", "beta_direct", "se_direct", "p_direct",
  "fdr_direct", "gene_ir_PP4", "indirect_beta", "indirect_se", "indirect_pval",
  "indirect_fdr", "indirect_direction_matches_direct"
))

sn_qc <- read_sheet(snrna_file, "09_pseudobulk_QC")
s12a <- sn_qc[, .(
  total_nuclei = sum(nuclei, na.rm = TRUE),
  contributing_donors = uniqueN(sample_name),
  pseudobulks_ge_20_nuclei = sum(nuclei >= 20, na.rm = TRUE),
  older_donors = uniqueN(sample_name[group == "Old"]),
  younger_donors = uniqueN(sample_name[group == "Young"])
), by = cell_type]
s12b <- select_columns(read_sheet(snrna_file, "01_sample_QC"), c(
  "sample_name", "geo_accession", "group", "nuclei_input", "nuclei_retained",
  "median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt"
))
s12c <- select_columns(sn_qc, c(
  "sample_name", "pseudobulk_id", "cell_type", "nuclei", "geo_accession", "group"
))
s12d <- select_columns(read_sheet(snrna_file, "07_TierA_localization"), c(
  "gene", "cell_type", "nuclei", "mean_log_normalized_expression", "percent_detected",
  "localization_rank"
))
s12e <- select_columns(read_sheet(snrna_file, "08_TierA_pseudobulk_DE"), c(
  "gene", "cell_type", "young_donors", "old_donors", "logFC_Old_vs_Young", "PValue",
  "candidate_FDR_global", "test_status"
))
s12f <- select_columns(read_sheet(snrna_file, "05_composition_tests"), c(
  "cell_type", "young_donors", "old_donors", "young_median_proportion",
  "old_median_proportion", "median_difference_Old_minus_Young", "PValue", "FDR"
))
cell_order <- c("Slow skeletal fiber", "Fibro-adipogenic progenitor", "Fast skeletal fiber",
                "Satellite cell", "Smooth muscle/pericyte", "Endothelial cell", "Immune cell")
gene_order <- c("ABCC8", "MAPK1", "YWHAZ", "ZBTB7B", "SMAD3", "RXRA")
s12a[, cell_order__ := match(cell_type, cell_order)]
setorder(s12a, cell_order__)
s12a[, cell_order__ := NULL]
setorder(s12b, sample_name)
s12c[, cell_order__ := match(cell_type, cell_order)]
setorder(s12c, sample_name, cell_order__)
s12c[, cell_order__ := NULL]
s12d[, `:=`(gene_order__ = match(gene, gene_order), cell_order__ = match(cell_type, cell_order))]
setorder(s12d, gene_order__, cell_order__)
s12d[, c("gene_order__", "cell_order__") := NULL]
s12e[, `:=`(gene_order__ = match(gene, gene_order), cell_order__ = match(cell_type, cell_order))]
setorder(s12e, gene_order__, cell_order__)
s12e[, c("gene_order__", "cell_order__") := NULL]
s12f[, cell_order__ := match(cell_type, cell_order)]
setorder(s12f, cell_order__)
s12f[, cell_order__ := NULL]

s13a <- select_columns(read_sheet(sqtl_file, "03_significant_sQTL_pairs"), c(
  "group_id", "phenotype_id", "variant_id", "start_distance", "af", "ma_samples", "ma_count",
  "pval_nominal", "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal",
  "pval_beta", "gene_symbol", "F_stat", "conventional_p_lt_5e_8", "empirical_pair_significant"
))
s13b <- select_columns(read_sheet(sqtl_file, "04_eligible_lead_sQTL"), c(
  "phenotype_id", "gene_id", "gene_name", "biotype", "gene_chr", "gene_start", "gene_end",
  "strand", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant_id",
  "tss_distance", "chr", "variant_pos", "ref", "alt", "num_alt_per_site",
  "rs_id_dbSNP157_GRCh38p14", "ma_samples", "ma_count", "af", "pval_nominal", "slope",
  "slope_se", "pval_perm", "pval_beta", "group_size", "qval", "pval_nominal_threshold",
  "gene_symbol", "gene_ensembl", "SNP", "effect_allele", "other_allele", "beta", "se",
  "eaf", "F_stat", "target_trait", "sgene_fdr_significant", "lead_pair_empirically_significant"
))
s13c <- select_columns(read_sheet(sqtl_file, "05_harmonisation"), c(
  "gene_symbol", "gene_ensembl", "phenotype_id", "trait", "SNP", "effect_allele",
  "other_allele", "exposure_beta", "exposure_se", "exposure_eaf", "outcome_beta_aligned",
  "outcome_se", "outcome_eaf_aligned", "af_difference", "harmonisation_action",
  "harmonisation_status"
))
s13d <- select_columns(read_sheet(sqtl_file, "06_targeted_sQTL_MR"), c(
  "gene_symbol", "gene_ensembl", "phenotype_id", "trait", "SNP", "method", "nsnp", "beta",
  "se", "exposure_beta", "exposure_se", "outcome_beta", "outcome_se", "F_stat",
  "exposure_effect_allele", "exposure_other_allele", "exposure_eaf", "outcome_eaf",
  "af_difference", "ci_lower", "ci_upper", "pval", "fdr_targeted", "splicing_direction"
))
s13e <- select_columns(read_sheet(sqtl_file, "07_layer_summary"), c(
  "gene_symbol", "target_trait", "gene_ensembl", "lead_phenotype_id", "lead_sQTL", "sgene_qval",
  "lead_sQTL_pval", "lead_sQTL_F", "sgene_fdr_significant", "n_reported_significant_pairs",
  "n_reported_splice_events", "minimum_nominal_p", "n_pairs_p_lt_5e_8", "sqtl_mr_beta",
  "sqtl_mr_se", "sqtl_mr_pval", "sqtl_mr_fdr"
))
s13a[, gene_order__ := match(gene_symbol, gene_order)]
setorder(s13a, gene_order__, phenotype_id, pval_nominal)
s13a[, gene_order__ := NULL]
s13b[, gene_order__ := match(gene_symbol, gene_order)]
setorder(s13b, gene_order__, phenotype_id)
s13b[, gene_order__ := NULL]
s13c[, gene_order__ := match(gene_symbol, gene_order)]
setorder(s13c, gene_order__, phenotype_id, trait)
s13c[, gene_order__ := NULL]
s13d[, gene_order__ := match(gene_symbol, gene_order)]
setorder(s13d, gene_order__, phenotype_id, trait)
s13d[, gene_order__ := NULL]
s13e[, gene_order__ := match(gene_symbol, gene_order)]
setorder(s13e, gene_order__)
s13e[, gene_order__ := NULL]

table_specs <- list(
  S10_MAGMA = list(title = "Supplementary Table S10. Corrected pathway-level MAGMA results.",
                  sections = list(list(label = "Results", data = s1))),
  S9_Cross_trait = list(title = "Supplementary Table S9. Corrected paired cross-phenotype MAGMA comparisons.",
                        sections = list(list(label = "A. Paired leave-one-chromosome-out contrasts", data = s2a),
                                        list(label = "B. Residual-scale sensitivity analyses", data = s2b))),
  S3_Candidate_sensitivity = list(title = "Supplementary Table S3. Candidate-set sensitivity analyses.",
                                  sections = list(list(label = "A. Prespecified candidate-set definitions", data = s3a),
                                                  list(label = "B. Leave-one-pathway-out analyses", data = s3b),
                                                  list(label = "C. Size-matched random gene sets", data = s3c))),
  S1_Candidate_membership = list(title = "Supplementary Table S1. Candidate pathways and gene membership.",
                                 sections = list(list(label = "A. Pathway-level annotation", data = s4a),
                                                 list(label = "B. Pathway-gene membership", data = s4b),
                                                 list(label = "C. Deduplicated candidate genes", data = s4c))),
  S11_All_cisMR = list(title = "Supplementary Table S11. Complete skeletal-muscle cis-MR results.", sections = list(list(label = "Results", data = s5))),
  S12_FDR_cisMR = list(title = "Supplementary Table S12. FDR-significant skeletal-muscle cis-MR results.", sections = list(list(label = "Results", data = s6))),
  S13_All_coloc = list(title = "Supplementary Table S13. Complete single-signal colocalization and prior sensitivity.", sections = list(list(label = "Results", data = s7))),
  S4_SuSiE = list(title = "Supplementary Table S4. Complete-region coloc-SuSiE results and LD audit.",
                  sections = list(list(label = "A. Locus-level QC and summary", data = s8a),
                                  list(label = "B. Signal-pair results across priors", data = s8b),
                                  list(label = "C. Credible-set summaries", data = s8c),
                                  list(label = "D. Credible-set members", data = s8d))),
  S14_FUSION = list(title = "Supplementary Table S14. Alternative skeletal-muscle eQTL-source analysis.", sections = list(list(label = "Results", data = s9))),
  S15_GEFOS = list(title = "Supplementary Table S15. GEFOS external-outcome follow-up.", sections = list(list(label = "Results", data = s10))),
  S16_MAGIC = list(title = "Supplementary Table S16. Exploratory MAGIC systemic-IR boundary analysis.", sections = list(list(label = "Results", data = s11))),
  S7_snRNA = list(title = "Supplementary Table S7. Single-nucleus localization and pseudobulk analyses.",
                   sections = list(list(label = "A. Cell-type donor and nuclei counts", data = s12a),
                                   list(label = "B. Donor-level QC", data = s12b),
                                   list(label = "C. Pseudobulk sample-cell-type QC", data = s12c),
                                   list(label = "D. Tier A cell-type localization", data = s12d),
                                   list(label = "E. Donor-aware pseudobulk contrasts", data = s12e),
                                   list(label = "F. Donor-level cell-type composition comparisons", data = s12f))),
  S8_sQTL = list(title = "Supplementary Table S8. GTEx v11 skeletal-muscle sQTL-MR results.",
                  sections = list(list(label = "A. Empirically significant sQTL pairs", data = s13a),
                                  list(label = "B. Eligible lead sQTLs", data = s13b),
                                  list(label = "C. Harmonization audit", data = s13c),
                                  list(label = "D. Targeted sQTL-MR", data = s13d),
                                  list(label = "E. Integrated summary", data = s13e)))
)

table_specs <- table_specs[order(as.integer(sub("^S([0-9]+).*", "\\1", names(table_specs))))]

title_style <- createStyle(fontColour = "#000000", textDecoration = "bold", fontSize = 12)
section_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#365F91", textDecoration = "bold")
header_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
                            halign = "center", valign = "center", wrapText = TRUE,
                            border = "TopBottomLeftRight", borderColour = "#D9E2F3")
body_style <- createStyle(fontColour = "#000000", valign = "top", border = "Bottom",
                          borderColour = "#E7E6E6")
number_style <- createStyle(numFmt = "0.00000E+00")
revision_style <- createStyle(fontColour = "#000000", textDecoration = "bold")

wb <- createWorkbook()
section_rows <- list()
for (sheet_name in names(table_specs)) {
  spec <- table_specs[[sheet_name]]
  addWorksheet(wb, sheet_name, gridLines = FALSE)
  writeData(wb, sheet_name, spec$title, startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, stack = TRUE)
  row <- 3L
  section_rows[[sheet_name]] <- list()
  for (section in spec$sections) {
    dat <- as.data.frame(section$data, check.names = FALSE)
    writeData(wb, sheet_name, section$label, startRow = row, startCol = 1, colNames = FALSE)
    mergeCells(wb, sheet_name, cols = 1:max(1, ncol(dat)), rows = row)
    addStyle(wb, sheet_name, section_style, rows = row, cols = 1:max(1, ncol(dat)), stack = TRUE)
    header_row <- row + 1L
    writeData(wb, sheet_name, dat, startRow = header_row, startCol = 1, withFilter = FALSE)
    addStyle(wb, sheet_name, header_style, rows = header_row, cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE)
    if (nrow(dat)) {
      body_rows <- (header_row + 1L):(header_row + nrow(dat))
      addStyle(wb, sheet_name, body_style, rows = body_rows, cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE)
      numeric_cols <- which(vapply(dat, is.numeric, logical(1)))
      if (length(numeric_cols)) addStyle(wb, sheet_name, number_style, rows = body_rows, cols = numeric_cols, gridExpand = TRUE, stack = TRUE)
    }
    section_rows[[sheet_name]][[section$label]] <- list(header = header_row, n = nrow(dat), columns = names(dat))
    row <- header_row + nrow(dat) + 2L
  }
  freezePane(wb, sheet_name, firstActiveRow = 5, firstActiveCol = 2)
  max_cols <- max(vapply(spec$sections, function(z) ncol(z$data), integer(1)))
  width_for_col <- function(col_index) {
    values <- unlist(lapply(spec$sections, function(z) {
      if (ncol(z$data) < col_index) return(character())
      c(names(z$data)[col_index], as.character(z$data[[col_index]]))
    }), use.names = FALSE)
    pmin(28, pmax(11, max(nchar(values), na.rm = TRUE) + 1))
  }
  setColWidths(wb, sheet_name, cols = seq_len(max_cols),
               widths = vapply(seq_len(max_cols), width_for_col, numeric(1)))
}

style_revision <- function(sheet, section_label, data, row_selector, columns) {
  meta <- section_rows[[sheet]][[section_label]]
  selected_rows <- which(row_selector)
  selected_cols <- match(columns, meta$columns)
  selected_cols <- selected_cols[!is.na(selected_cols)]
  if (length(selected_rows) && length(selected_cols)) {
    addStyle(wb, sheet, revision_style,
             rows = meta$header + selected_rows, cols = selected_cols,
             gridExpand = TRUE, stack = TRUE)
  }
}

# Black bold marks values recalculated during this revision. Identifiers and
# unchanged descriptive fields remain regular black text.
style_revision("S10_MAGMA", "Results", s1, s1$trait == "ALM",
               c("BETA", "BETA_STD", "SE", "P", "FDR"))
style_revision("S9_Cross_trait", "A. Paired leave-one-chromosome-out contrasts", s2a,
               s2a$comparison %in% c("ALM vs GRIP", "ALM vs WALK"),
               setdiff(names(s2a), c("comparison", "trait_1", "trait_2")))
style_revision("S9_Cross_trait", "B. Residual-scale sensitivity analyses", s2b,
               s2b$comparison %in% c("ALM vs GRIP", "ALM vs WALK"),
               setdiff(names(s2b), c("scale", "comparison", "trait_1", "trait_2")))
style_revision("S3_Candidate_sensitivity", "A. Prespecified candidate-set definitions", s3a,
               rep(TRUE, nrow(s3a)), c("beta", "beta_std", "se", "t", "p_positive", "fdr_positive_BH"))
style_revision("S3_Candidate_sensitivity", "B. Leave-one-pathway-out analyses", s3b,
               rep(TRUE, nrow(s3b)), c("beta", "beta_std", "se", "t", "p_positive", "fdr_positive_BH"))
style_revision("S3_Candidate_sensitivity", "C. Size-matched random gene sets", s3c,
               rep(TRUE, nrow(s3c)), setdiff(names(s3c), c("set_name", "n_iterations", "n_target_magma_eligible")))
style_revision("S15_GEFOS", "Results", s10, rep(TRUE, nrow(s10)),
               c("external_coloc_PP3", "external_coloc_PP4"))
style_revision("S7_snRNA", "A. Cell-type donor and nuclei counts", s12a,
               rep(TRUE, nrow(s12a)), setdiff(names(s12a), "cell_type"))
style_revision("S8_sQTL", "B. Eligible lead sQTLs", s13b,
               rep(TRUE, nrow(s13b)), c("pval_nominal", "slope", "slope_se", "pval_beta", "qval", "F_stat"))
style_revision("S8_sQTL", "D. Targeted sQTL-MR", s13d,
               rep(TRUE, nrow(s13d)), c("beta", "se", "pval", "fdr_targeted"))
style_revision("S8_sQTL", "E. Integrated summary", s13e,
               rep(TRUE, nrow(s13e)), c("sgene_qval", "lead_sQTL_pval", "lead_sQTL_F",
                                        "minimum_nominal_p", "sqtl_mr_beta", "sqtl_mr_se",
                                        "sqtl_mr_pval", "sqtl_mr_fdr"))

supp_table_file <- file.path(submission_dir, "Supplementary_tables_13_main.xlsx")
saveWorkbook(wb, supp_table_file, overwrite = TRUE)

repair_submission_xlsx <- function(xlsx_file) {
  if (!requireNamespace("xml2", quietly = TRUE) || !requireNamespace("zip", quietly = TRUE)) {
    stop("Workbook compatibility repair requires xml2 and zip.")
  }
  repair_dir <- tempfile("submission_xlsx_")
  repaired <- tempfile(fileext = ".xlsx")
  dir.create(repair_dir)
  on.exit(unlink(repair_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(repaired, force = TRUE), add = TRUE)
  utils::unzip(xlsx_file, exdir = repair_dir)
  archive_files <- gsub("\\\\", "/", list.files(repair_dir, recursive = TRUE, all.files = TRUE,
                                                     no.. = TRUE, include.dirs = FALSE))
  normalise_part <- function(path) {
    parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
    out <- character()
    for (part in parts) {
      if (!nzchar(part) || part == ".") next
      if (part == "..") out <- head(out, -1L) else out <- c(out, part)
    }
    paste(out, collapse = "/")
  }
  for (rel_file in list.files(repair_dir, pattern = "\\.rels$", recursive = TRUE,
                              full.names = TRUE, all.files = TRUE)) {
    rel_path <- substring(gsub("\\\\", "/", rel_file), nchar(gsub("\\\\", "/", repair_dir)) + 2L)
    rel_base <- if (rel_path == "_rels/.rels") "" else sub("/_rels/[^/]+\\.rels$", "", rel_path)
    rel_xml <- xml2::read_xml(rel_file)
    nodes <- xml2::xml_find_all(rel_xml, "//*[local-name()='Relationship']")
    for (node in nodes) {
      if (identical(xml2::xml_attr(node, "TargetMode"), "External")) next
      target <- xml2::xml_attr(node, "Target")
      part <- if (startsWith(target, "/")) sub("^/+", "", target) else normalise_part(file.path(rel_base, target))
      if (!part %in% archive_files) xml2::xml_remove(node)
    }
    xml2::write_xml(rel_xml, rel_file)
  }
  content_file <- file.path(repair_dir, "[Content_Types].xml")
  content_xml <- xml2::read_xml(content_file)
  for (node in xml2::xml_find_all(content_xml, "//*[local-name()='Override']")) {
    part <- sub("^/+", "", xml2::xml_attr(node, "PartName"))
    if (!part %in% archive_files) xml2::xml_remove(node)
  }
  xml2::write_xml(content_xml, content_file)
  col_number <- function(x) sum(match(strsplit(x, "")[[1]], LETTERS) * 26 ^ rev(seq_len(nchar(x)) - 1L))
  col_label <- function(x) {
    out <- character()
    while (x > 0L) {
      remainder <- (x - 1L) %% 26L
      out <- c(LETTERS[remainder + 1L], out)
      x <- (x - remainder - 1L) %/% 26L
    }
    paste(out, collapse = "")
  }
  for (sheet_file in list.files(file.path(repair_dir, "xl", "worksheets"),
                                pattern = "^sheet[0-9]+\\.xml$", full.names = TRUE)) {
    sheet_xml <- xml2::read_xml(sheet_file)
    refs <- xml2::xml_attr(xml2::xml_find_all(sheet_xml, "//*[local-name()='c']"), "r")
    refs <- refs[!is.na(refs)]
    if (!length(refs)) next
    cols <- gsub("[0-9]", "", refs)
    rows <- as.integer(gsub("[A-Z]", "", refs))
    dimension <- xml2::xml_find_first(sheet_xml, "//*[local-name()='dimension']")
    xml2::xml_set_attr(dimension, "ref", paste0("A1:", col_label(max(vapply(cols, col_number, numeric(1)))), max(rows)))
    xml2::write_xml(sheet_xml, sheet_file)
  }
  old_wd <- setwd(repair_dir)
  on.exit(setwd(old_wd), add = TRUE)
  zip::zipr(repaired, list.files(".", all.files = TRUE, no.. = TRUE), recurse = TRUE,
            include_directories = FALSE)
  setwd(old_wd)
  if (!file.copy(repaired, xlsx_file, overwrite = TRUE)) stop("Could not replace repaired workbook.")
}
repair_submission_xlsx(supp_table_file)

theme_submission <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", colour = "#111111"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      strip.background = element_rect(fill = "#E8EDF2", colour = "#AAB4BE"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_pair <- function(plot, stem, width, height) {
  png_file <- file.path(supp_figure_dir, paste0(stem, ".png"))
  pdf_file <- file.path(supp_figure_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(supp_figure_dir, paste0(stem, ".svg"))
  png_device <- function(filename, width, height, ...) {
    grDevices::png(filename, width = width, height = height, units = "in", res = 450)
  }
  ggsave(png_file, plot, width = width, height = height, dpi = 450,
         device = png_device, bg = "white")
  grDevices::cairo_pdf(pdf_file, width = width, height = height,
                       family = "Arial", bg = "white")
  print(plot)
  grDevices::dev.off()
  svglite::svglite(svg_file, width = width, height = height, bg = "white")
  print(plot)
  grDevices::dev.off()
}

cross_plot <- copy(s2a[comparison %in% c("ALM vs GRIP", "ALM vs WALK")])
cross_plot[, label := fifelse(comparison == "ALM vs GRIP", "ALM vs grip strength", "ALM vs walking pace")]
cross_plot[, label := factor(label, levels = rev(unique(label)))]
p_s1a <- ggplot(cross_plot, aes(delta_beta, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#6B7280") +
  geom_errorbar(aes(xmin = ci_95_lower, xmax = ci_95_upper), width = 0.15, colour = "#2F5D7C") +
  geom_point(size = 2.7, colour = "#B23A48") +
  labs(title = "a  Primary LOCO contrasts", x = expression(Delta*beta), y = NULL) +
  theme_submission()

resid <- copy(s2b[comparison %in% c("ALM vs GRIP", "ALM vs WALK")])
resid[, comparison_label := fifelse(comparison == "ALM vs GRIP", "ALM vs grip", "ALM vs walking pace")]
resid[, scale_label := factor(scale,
  levels = c("raw_ZRESID_BASE", "within_trait_z_standardized", "within_trait_rank_normalized"),
  labels = c("Raw residuals", "Z standardized", "Rank normalized"))]
p_s1b <- ggplot(resid, aes(scale_label, comparison_label, fill = residual_contrast_beta)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f\nq=%.3g", residual_contrast_beta, fdr_BH)), size = 3) +
  scale_fill_gradient2(low = "#3B78A8", mid = "white", high = "#B23A48", midpoint = 0) +
  labs(title = "b  Residual-scale sensitivity", x = NULL, y = NULL, fill = expression(Delta*beta)) +
  theme_submission() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_pair(p_s1a + p_s1b + plot_layout(widths = c(0.9, 1.35)), "Figure_S1_Cross_trait_MAGMA", 10, 4.2)

sets <- copy(s3a)
sets[, set_label := gsub("_", " ", set_name)]
sets[, set_label := factor(set_label, levels = rev(set_label))]
p_s2a <- ggplot(sets, aes(beta, set_label, colour = sensitivity_family)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#6B7280") +
  geom_errorbar(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se), width = 0.14) +
  geom_point(size = 2.3) +
  scale_colour_manual(values = c("primary_full_set" = "#2F5D7C", "core_pathway_set" = "#4F7F52",
                                 "broad_pathway_exclusion" = "#B23A48")) +
  labs(title = "a  Candidate-set definitions", x = "MAGMA beta (95% CI)", y = NULL, colour = NULL) +
  theme_submission(8) + theme(legend.position = "none")

random <- copy(s3c)
p_s2c <- ggplot(random, aes(y = 1)) +
  geom_segment(aes(x = random_beta_q025, xend = random_beta_q975, yend = 1), linewidth = 2.8, colour = "#AAB4BE") +
  geom_point(aes(x = random_beta_q500), shape = 21, fill = "white", size = 3) +
  geom_point(aes(x = observed_beta), colour = "#B23A48", size = 3.2) +
  annotate("text", x = random$observed_beta, y = 1.12,
           label = paste0("Observed; empirical P=", format(random$empirical_p_positive, digits = 3)),
           colour = "#B23A48", size = 3, hjust = 0.5) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.30))) +
  labs(title = "c  Size-matched random gene sets", x = "MAGMA beta") +
  theme_submission(8) + theme(legend.position = "none")

lopo <- copy(s3b)
lopo[, label := gsub("^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_", "", removed_pathway)]
lopo[, label := gsub("_", " ", label)]
lopo[, label := vapply(label, function(x) paste(strwrap(x, width = 34), collapse = "\n"), character(1))]
lopo[, display_order := .I]
lopo_left <- lopo[display_order <= ceiling(nrow(lopo) / 2)]
lopo_right <- lopo[display_order > ceiling(nrow(lopo) / 2)]
make_lopo_panel <- function(dat, panel_title = NULL, show_legend = FALSE) {
  dat <- copy(dat)
  dat[, label := factor(label, levels = rev(label))]
  ggplot(dat, aes(beta, label, colour = removed_pathway_role)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#6B7280") +
    geom_errorbar(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se), width = 0.12) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("core_insulin_response" = "#2F5D7C",
                                   "downstream_or_contextual" = "#4F7F52",
                                   "broad_metabolic_or_disease" = "#B23A48"),
                        breaks = c("core_insulin_response", "downstream_or_contextual",
                                   "broad_metabolic_or_disease"),
                        labels = c("Core", "Contextual", "Broad")) +
    coord_cartesian(xlim = c(-0.015, 0.32)) +
    labs(title = panel_title, x = "MAGMA beta (95% CI)", y = NULL, colour = "Pathway class") +
    theme_submission(7.2) +
    theme(axis.text.y = element_text(size = 6.1), legend.position = if (show_legend) "bottom" else "none")
}
p_s2b_left <- make_lopo_panel(lopo_left, "b  Leave-one-pathway-out estimates", FALSE)
p_s2b_right <- make_lopo_panel(lopo_right, NULL, TRUE)
save_pair((p_s2a | p_s2c) / (p_s2b_left | p_s2b_right) + plot_layout(heights = c(0.9, 1.45)),
          "Figure_S2_Candidate_set_robustness", 11, 7.2)

gefos_plot <- copy(s10[!is.na(external_beta)])
gefos_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
p_s13a <- ggplot(gefos_plot, aes(external_beta, gene_symbol)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#6B7280") +
  geom_errorbar(aes(xmin = external_ci_lower, xmax = external_ci_upper), width = 0.15, colour = "#2F5D7C") +
  geom_point(aes(shape = external_p_bonferroni < 0.05), size = 2.8, colour = "#B23A48") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19), guide = "none") +
  labs(title = "a  External ALM cis-MR", x = "Effect estimate (95% CI)", y = NULL) +
  theme_submission()
p_s13b <- ggplot(gefos_plot, aes(external_coloc_PP4, gene_symbol)) +
  geom_vline(xintercept = 0.8, linetype = 2, colour = "#6B7280") +
  geom_segment(aes(x = 0, xend = external_coloc_PP4, yend = gene_symbol), colour = "#AAB4BE") +
  geom_point(size = 2.8, colour = "#2F5D7C") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "b  External-outcome colocalization", x = "PP4", y = NULL) +
  theme_submission()
save_pair(p_s13a + p_s13b, "Figure_S13_GEFOS_external_followup", 9, 4.4)

supp_sources <- c(
  Figure_S3_cisMR_forest = "Figure_cisMR_forest",
  Figure_S4_coloc_PP3_PP4 = "Figure_coloc_PP3_PP4",
  Figure_S5_ABCC8_ALM = "coloc_multisignal_locus_ABCC8_ALM",
  Figure_S6_MAPK1_ALM = "coloc_multisignal_locus_MAPK1_ALM",
  Figure_S7_YWHAZ_ALM = "coloc_multisignal_locus_YWHAZ_ALM",
  Figure_S8_ZBTB7B_ALM = "coloc_multisignal_locus_ZBTB7B_ALM",
  Figure_S9_SMAD3_GRIP = "coloc_multisignal_locus_SMAD3_GRIP",
  Figure_S10_RXRA_WALK = "coloc_multisignal_locus_RXRA_WALK",
  Figure_S11_snRNA_UMAP = "GSE167186_snRNA_celltype_UMAP",
  Figure_S12_TierA_celltype_expression = "GSE167186_TierA_celltype_expression"
)
for (out_stem in names(supp_sources)) {
  for (ext in c("png", "pdf")) {
    src <- find_input(paste0(supp_sources[[out_stem]], ".", ext), required = (ext == "png"))
    if (!is.na(src)) {
      file.copy(src, file.path(supp_figure_dir, paste0(out_stem, ".", ext)), overwrite = TRUE)
    } else {
      stale_output <- file.path(supp_figure_dir, paste0(out_stem, ".", ext))
      if (file.exists(stale_output)) unlink(stale_output)
    }
  }
}

for (stem in c("Figure1_Study_workflow", "Figure2_TierA_cisMR_forest",
               "Figure3_Cell_and_splicing_context", "Figure4_Multilayer_evidence_matrix")) {
  for (ext in c("png", "pdf")) {
    src <- find_input(paste0(stem, ".", ext))
    file.copy(src, file.path(main_figure_dir, paste0(stem, ".", ext)), overwrite = TRUE)
  }
}

supp_png_stems <- paste0("Figure_S", 1:13, c(
  "_Cross_trait_MAGMA", "_Candidate_set_robustness", "_cisMR_forest", "_coloc_PP3_PP4",
  "_ABCC8_ALM", "_MAPK1_ALM", "_YWHAZ_ALM", "_ZBTB7B_ALM", "_SMAD3_GRIP",
  "_RXRA_WALK", "_snRNA_UMAP", "_TierA_celltype_expression", "_GEFOS_external_followup"
))
expected_png <- c(
  file.path(supp_figure_dir, paste0(supp_png_stems, ".png")),
  file.path(main_figure_dir, paste0(c("Figure1_Study_workflow", "Figure2_TierA_cisMR_forest",
                                     "Figure3_Cell_and_splicing_context", "Figure4_Multilayer_evidence_matrix"), ".png")),
  supp_table_file
)
missing_outputs <- expected_png[!file.exists(expected_png)]
if (length(missing_outputs)) stop("Submission output validation failed: ", paste(basename(missing_outputs), collapse = ", "))
if (!identical(getSheetNames(supp_table_file), names(table_specs))) stop("Supplementary workbook sheet validation failed.")

message("Submission-ready tables and figures created in: ", submission_dir)
