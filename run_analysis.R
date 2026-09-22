## PUBLIC RELEASE -- synchronized 2026-09-17
## Run from R: source("run_analysis.R", encoding = "UTF-8")
## Or: Rscript run_analysis.R
## Check prerequisites without running analyses: Rscript run_analysis.R --check
##
## [NEW: release notes]
## This single file reuses IR_sarcopenia.R, the GEFOS extension, and the manuscript
## output script. Original files are not sourced or modified. Changed blocks are
## marked [NEW] or [REVISED]. Place this script beside the original input folders.
## Outputs: results_public_release/ and manuscript_outputs_public_release/.
## Input data and licensed reference resources must be obtained separately.
## MAGMA/PLINK intermediate input files, caches and provenance logs are retained
## where required by the inherited analysis. Final result tables are XLSX.
## Clinical database queries need network access and may change over time.
##
## Reviewer 3: embedded archived membership (24 pathways, 771 unique symbols);
## formal MAGMA competitive alternative-set/leave-one-pathway-out tests; the
## matched-random residual-Z analysis is supplementary and uses a different scale.
## Exclusions remove ALL genes annotated to excluded pathways, including shared
## genes. BH correction is across five set definitions, or 24 deletion tests.
## Archived timestamps and today's checksum do NOT prove selection preceded every
## historical outcome analysis. The fixed input establishes reproducible reruns.
## Reviewer 5: complete LD-matched regions; stop above the 8000-variant memory
## guard instead of selecting variants by association P value. Report all priors,
## PP3, PP4, credible sets, member variants and non-interpretable fits separately.
## Reviewer 6: GTEx v8 positions are GRCh38; primary GWAS positions are GRCh37.
## Main MR/coloc join by rsID with allele checks, not by equal coordinate values.
## Explicit hg38->hg19 liftOver is used in the inherited MAGIC extension.
## Primary MR/eQTL import keeps the first valid occurrence of each rsID;
## the MAGMA branch re-reads raw outcome rows and resolves duplicate rsIDs
## against 1000G EUR chromosome, position and unordered allele pairs. GEFOS
## sorts by P then descending N before deduplication. These distinct policies
## are reported explicitly, not described as one universal min-P rule.
## Regional plots place both datasets at the matched GWAS GRCh37 position.
## MR uses harmonise_data(action=3); coloc retains matching/swapped alleles and
## removes palindromic variants when either MAF >0.42 (see implementation).
## Reviewer 7: reuse GTEx instruments; four external ALM targets including missing
## targets in multiplicity correction; external single-signal coloc at all priors.
## This external coloc implementation explicitly reuses the primary sdY assumptions.
## Its PP4 values differ slightly from the earlier supplement without a located
## generating script. Use the regenerated tables consistently when reporting it.
## External support is limited, not definitive replication. snRNA age comparisons
## are exploratory donor-level, group-only edgeR models (17 male donors); sex is
## invariant, batch/comorbidity/activity terms are unavailable and not adjusted.
## Low detection does not establish absence of expression.
## Figures 1/3/4 and Figure 2 legend use the current revised manuscript code.
## No renv, automatic package installation or repository operations are used.
##
## Archived membership SHA256: d5bb47e1b9e64a823e43d1942c928f353299898cbe81b27582ecb31d118271f0
## Source IR_sarcopenia.R SHA256: 08599b4f0f437b5eb8c8e6beaee7f2fb864ef03383d2276b3b25314f2971f564
## Source IR_sarcopenia_with_GEFOS_ALM_validation_revised.R SHA256: 40cb8b74f4b8a30d5987e7670740af9a81ff5238b9e0337ce3bdfb4c26aaeb21
## Source IR_sarcopenia_manuscript_outputs_revised.R SHA256: 1ea16086ee80f05d582296ea014b633972781f12724ee94083212619c4c036ba


## =============================================================================
## IR-related skeletal-muscle genes and sarcopenia traits
## Refactored analysis pipeline
##
## Main evidence chain:
## candidate genes -> skeletal-muscle cis-eQTL -> cis-MR -> colocalisation
##
## Important interpretation:
## - FDR-significant cis-MR results are "MR-associated genes".
## - Genes with strong colocalisation support are "prioritised candidate genes".
## - Cross-gene IVW and cross-gene Cochran Q/I2 are intentionally not used as
##   causal or MR-heterogeneity tests in this script.
## =============================================================================


# [NEW] Windows R needs a UTF-8 character locale for non-ASCII project paths.
if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

## 00. Packages ---------------------------------------------------------------

required_packages <- c(
  "data.table", "dplyr", "stringr", "msigdbr",
  "org.Hs.eg.db", "AnnotationDbi", "AnnotationFilter",
  "EnsDb.Hsapiens.v86", "ensembldb",
  "GenomicRanges", "IRanges", "Rsamtools", "TwoSampleMR", "ieugwasr",
  "coloc", "susieR", "openxlsx", "clusterProfiler", "enrichplot",
  "ggplot2", "ggrepel")


missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]


if (length(missing_packages) > 0L) {
  stop(
    "Please install the following packages before running the script: ",
    paste(missing_packages, collapse = ", "))}


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(AnnotationFilter)
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  library(GenomicRanges)
  library(TwoSampleMR)
  library(coloc)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)})


## 01. Config -----------------------------------------------------------------

get_script_dir <- function() {
  command_file <- sub("^--file=", "", grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )[1])
  if (!is.na(command_file) && nzchar(command_file)) {
    return(normalizePath(dirname(command_file), winslash = "/", mustWork = TRUE))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    return(normalizePath(dirname(frame_files[[length(frame_files)]]),
                         winslash = "/", mustWork = TRUE))
  }

  normalizePath(".", winslash = "/", mustWork = TRUE)
}

resolve_project_dir <- function() {
  external_data_dir <- Sys.getenv("SARC_IR_DATA_DIR", unset = "")
  if (nzchar(external_data_dir)) {
    return(normalizePath(external_data_dir, winslash = "/", mustWork = TRUE))
  }
  get_script_dir()
}

config <- list(
  project_dir = resolve_project_dir(), # [REVISED] resolve inputs beside this file, or SARC_IR_DATA_DIR when set

  files = list(
    outcomes = c(
      ALM      = "outcome/ebi-a-GCST90000025_ALM.vcf",
      GRIP     = "outcome/ukb-b-10215_right_grip.vcf",
      WALK     = "outcome/ukb-b-4711_walkpace.vcf",
      LOW_GRIP = "outcome/ebi-a-GCST90007529_low_grip.vcf"
    ),
    gtex_muscle = "eQTL/Muscle_Skeletal.tsv.gz",
    fusion_eqtl = "qtltools_nominal.tsv.gz",
    fusion_variants = "genotype-variant_information.tsv.gz",
    fusion_sample_size = "sample_size.tsv",
    eqtlgen = "Significant_cis_eQTLs.gz",
    external_alm = "appendicularleanmass.results.metal_.txt.gz"
  ),

  outcome_n = c(
    ALM      = 450243,
    GRIP     = 461089,
    WALK     = 459915,
    LOW_GRIP = 256523
  ),

  outcome_type = c(
    ALM = "quant", GRIP = "quant", WALK = "quant", LOW_GRIP = "cc"
  ),

  outcome_cases = c(
    ALM = NA_real_, GRIP = NA_real_, WALK = NA_real_, LOW_GRIP = 20335
  ),

  gtex = list(
    sample_size = 706,
    cis_window_bp = 1000000L,
    minimum_region_rows = 20L
  ),

  fusion = list(
    sample_size = 301L
  ),

  instruments = list(
    p_threshold = 5e-8,
    minimum_F = 10,
    clump_kb = 1000,
    clump_r2 = 0.001,
    ancestry = "EUR",
    clump_mode = "local",
    plink_binary = "plink/plink/plink.exe",
    ld_reference = "1000G_EUR/1000G.EUR.QC"
  ),

  multiple_testing = list(
    alpha = 0.05,
    primary_scope = "within_trait" # also reports global FDR
  ),

  coloc = list(
    run = TRUE,
    minimum_common_snps = 50L,
    minimum_maf = 0.01,
    p1 = 1e-4,
    p2 = 1e-4,
    p12_values = c(1e-6, 1e-5, 1e-4),
    primary_p12 = 1e-5,
    tier_a_pp4 = 0.80,
    tier_b_pp4 = 0.60,
    tier_c_pp4 = 0.30,
    run_susie = TRUE,
    # Reviewer 1 revision: all six Tier A loci contain fewer than 8,000
    # LD-reference allele-matched variants. The higher ceiling retains the
    # complete matched region and avoids association-based SNP truncation.
    susie_max_snps = 8000L),
  optional = list(
    run_fusion_replication = TRUE,
    run_eqtlgen_comparison = TRUE,
    run_go_enrichment = TRUE,
    make_locus_plots = FALSE,
    coordinates_harmonised_for_plots = FALSE),
  # [NEW: Reviewer 7] Fixed discovery-derived external targets.
  external_alm = list(
    run = TRUE,
    dataset_label = "GEFOS_2017_ALM",
    sample_size = 28330L,
    target_genes = c("ABCC8", "MAPK1", "YWHAZ", "ZBTB7B"),
    nominal_alpha = 0.05,
    bonferroni_alpha = 0.05 / 4,
    require_primary_tier_a = TRUE),
  random_seed = 20260719L)

config$paths <- list(
  results = file.path(config$project_dir, "results_public_release"),
  plots = file.path(config$project_dir, "results_public_release", "plots"),
  cache = file.path(config$project_dir, "eqtl_cache_gtex_v8_local"))

config$files$outcomes[] <- file.path(config$project_dir, config$files$outcomes)
for (file_name in setdiff(names(config$files), "outcomes")) {
  config$files[[file_name]] <- file.path(
    config$project_dir, config$files[[file_name]]
  )
}
config$instruments$plink_binary <- file.path(
  config$project_dir, config$instruments$plink_binary
)
config$instruments$ld_reference <- file.path(
  config$project_dir, config$instruments$ld_reference
)

# [NEW: MAGMA reproducibility] Locate the shared duplicate-rsID resolver from
# either the package-root entrypoint or its archived scripts/ copy.
release_root <- get_script_dir()
if (!file.exists(file.path(release_root, "scripts", "06_magma_variant_resolution.R"))) {
  release_root <- dirname(release_root)
}
magma_variant_resolution_script <- file.path(
  release_root, "scripts", "06_magma_variant_resolution.R"
)

# [NEW: public-release prerequisite check; does not execute analyses]
setwd(config$project_dir)
locus_to_gene_script <- file.path(
  release_root, "scripts", "07_locus_to_gene_analysis.R"
)
if ("--check" %in% commandArgs(trailingOnly = TRUE)) {
  extra_packages <- c("rtracklayer", "DESeq2", "edgeR", "Seurat", "SeuratObject",
                      "Matrix", "patchwork", "scales", "ragg", "httr2", "jsonlite",
                      "xml2", "zip", "arrow")
  absent_packages <- extra_packages[!vapply(extra_packages, requireNamespace,
                                           logical(1), quietly = TRUE)]
  inputs <- c(unname(config$files$outcomes), config$files$gtex_muscle,
    paste0(config$files$gtex_muscle, ".tbi"), config$files$fusion_eqtl,
    config$files$fusion_variants, config$files$external_alm,
    config$instruments$plink_binary, "magma.exe", "NCBI37.3.gene.loc",
    "HOMA_IR.vcf.gz", "ISI_adjBMI_EUR.gz", "FI_EUR.gz", "hg38ToHg19.over.chain.gz",
    "GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv.gz",
    "GSE111016_series_matrix.txt.gz", "GSE167186_counts.csv.gz",
    "GSE167186-GPL20301_series_matrix.txt.gz", "GSE167186_RAW.tar",
    "GSE167186_SimplifiedMetadataSheet.xlsx", "GTEx_Analysis_v11_sQTL.tar",
    "GTEx_Analysis_v11_sQTL_groups.tar", locus_to_gene_script,
    magma_variant_resolution_script,
    paste0(config$instruments$ld_reference, ".", rep(1:22, each = 3L),
           rep(c(".bed", ".bim", ".fam"), 22L)))
  absent_inputs <- inputs[!file.exists(inputs)]
  if (length(absent_packages) || length(absent_inputs)) stop(paste(
    "Prerequisites missing. Packages:", paste(absent_packages, collapse = ", "),
    "Files:", paste(absent_inputs, collapse = "\n"), sep = "\n"))
  message("Prerequisite check passed. No analyses were run; network services were not tested.")
  quit(save = "no", status = 0L)
}

if (!file.exists(magma_variant_resolution_script)) {
  stop("MAGMA duplicate-rsID resolver not found: ", magma_variant_resolution_script)
}
sys.source(magma_variant_resolution_script, envir = environment())

dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(config$paths$plots, recursive = TRUE, showWarnings = FALSE)
dir.create(config$paths$cache, recursive = TRUE, showWarnings = FALSE)

set.seed(config$random_seed)

## Never paste a JWT into this script. Put OPENGWAS_JWT=your_new_token in the
## user-level .Renviron file, then restart R.
if (config$instruments$clump_mode == "remote" &&
    !nzchar(Sys.getenv("OPENGWAS_JWT"))) {
  stop(
    "Remote LD clumping requires OPENGWAS_JWT. ",
    "Set it as an environment variable; do not paste it into this script."
  )
}


## 02. General helper functions ----------------------------------------------

log_step <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

result_path <- function(filename) {
  file.path(config$paths$results, filename)
}

plot_path <- function(filename) {
  file.path(config$paths$plots, filename)
}

result_tables <- new.env(parent = emptyenv())

# [REVISED: Reviewer 1 comment 6] Audit duplicate rsID handling without
# changing the retained analysis records.
variant_dedup_audit <- data.table::data.table()

deduplicate_variant_rows <- function(dt, source, scope) {
  dt <- data.table::as.data.table(dt)
  n_input <- nrow(dt)
  n_unique <- data.table::uniqueN(dt$SNP)
  n_removed <- n_input - n_unique
  variant_dedup_audit <<- data.table::rbindlist(list(
    variant_dedup_audit,
    data.table::data.table(
      source = source,
      scope = scope,
      n_input_rows = n_input,
      n_unique_rsids = n_unique,
      n_removed_duplicate_rows = n_removed,
      retention_rule = "first valid occurrence in source order"
    )
  ), fill = TRUE)
  dt[!duplicated(SNP)]
}

write_result <- function(x, filename) {
  x <- data.table::as.data.table(x)
  if (ncol(x) == 0L) {
    x <- data.table::data.table(status = "no_results")
  }
  sheet_name <- substr(tools::file_path_sans_ext(filename), 1L, 31L)
  base_name <- sheet_name
  suffix <- 1L
  while (exists(sheet_name, envir = result_tables, inherits = FALSE)) {
    suffix <- suffix + 1L
    sheet_name <- paste0(substr(base_name, 1L, 27L), "_", suffix)
  }
  assign(sheet_name, as.data.frame(x), envir = result_tables)
}

assert_files_exist <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(label, " file(s) not found:\n", paste(missing, collapse = "\n"))
  }
}

empty_dt <- function() data.table::data.table()

bind_nonempty <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (length(x) == 0L) return(empty_dt())
  data.table::rbindlist(x, fill = TRUE, use.names = TRUE)
}

strip_ensembl_version <- function(x) sub("\\..*$", "", as.character(x))

is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT", "TA", "CG", "GC")
}

normal_pvalue <- function(beta, se) {
  2 * stats::pnorm(-abs(beta / se))
}

theme_manuscript <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}


## 03. Data import and standardisation functions -----------------------------

vcf_to_gwas <- function(path, trait, fixed_n, deduplicate = TRUE,
                         audit_scope = "outcome_GWAS") {
  log_step("Reading outcome: ", trait, " <- ", path)
  dt <- data.table::fread(
    path,
    skip = "#CHROM",
    select = c(1L, 2L, 3L, 4L, 5L, 9L, 10L),
    na.strings = c(".", "NA"),
    showProgress = TRUE
  )
  data.table::setnames(
    dt,
    c("chr", "pos", "SNP", "other_allele", "effect_allele",
      "FORMAT", "sample")
  )

  dt[, `:=`(
    beta = NA_real_,
    se = NA_real_,
    eaf = NA_real_,
    N = as.numeric(fixed_n)
  )]

  format_templates <- unique(dt$FORMAT)
  format_templates <- format_templates[!is.na(format_templates)]
  for (format_template in format_templates) {
    row_index <- which(dt$FORMAT == format_template)
    format_fields <- strsplit(format_template, ":", fixed = TRUE)[[1]]
    format_values <- data.table::tstrsplit(
      dt$sample[row_index], ":", fixed = TRUE
    )

    get_format_value <- function(field) {
      field_index <- match(field, format_fields)
      if (is.na(field_index)) return(rep(NA_real_, length(row_index)))
      as.numeric(format_values[[field_index]])
    }

    dt[row_index, `:=`(
      beta = get_format_value("ES"),
      se = get_format_value("SE"),
      eaf = get_format_value("AF")
    )]
    if ("SS" %in% format_fields) {
      dt[row_index, N := get_format_value("SS")]
    }
  }
  dt[, c("FORMAT", "sample") := NULL]

  dt[!is.finite(N) | is.na(N), N := fixed_n]
  dt[, pval := normal_pvalue(beta, se)]
  dt[, maf := pmin(eaf, 1 - eaf)]
  dt <- dt[
    grepl("^rs", SNP) & is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(eaf) & eaf > 0 & eaf < 1 &
      !is.na(effect_allele) & !is.na(other_allele) &
      !grepl(",", effect_allele, fixed = TRUE)
  ]
  if (isTRUE(deduplicate)) {
    deduplicate_variant_rows(dt, source = trait, scope = audit_scope)
  } else {
    dt
  }
}

gtex_columns <- c(
  "variant", "r2", "pvalue", "molecular_trait_object_id",
  "molecular_trait_id", "maf", "gene_id", "median_tpm", "beta", "se",
  "an", "ac", "chromosome", "position", "ref", "alt", "type", "rsid"
)

read_local_eqtl_region <- function(target_gene_id, chromosome, start, end) {
  query <- GenomicRanges::GRanges(
    as.character(chromosome),
    IRanges::IRanges(
      max(1L, as.integer(start - config$gtex$cis_window_bp)),
      as.integer(end + config$gtex$cis_window_bp)
    )
  )
  lines <- Rsamtools::scanTabix(config$files$gtex_muscle, param = query)[[1]]
  if (length(lines) == 0L) return(empty_dt())
  lines <- lines[grepl(target_gene_id, lines, fixed = TRUE)]
  if (length(lines) == 0L) return(empty_dt())

  dt <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    sep = "\t",
    showProgress = FALSE
  )
  data.table::setnames(dt, gtex_columns)
  dt <- dt[
    strip_ensembl_version(gene_id) == strip_ensembl_version(target_gene_id)
  ]
  if (nrow(dt) == 0L) return(empty_dt())

  dt[, .(
    SNP = rsid,
    beta = as.numeric(beta),
    se = as.numeric(se),
    pval = as.numeric(pvalue),
    eaf = as.numeric(ac) / as.numeric(an),
    maf = as.numeric(maf),
    effect_allele = toupper(as.character(alt)),
    other_allele = toupper(as.character(ref)),
    pos = as.integer(position),
    chr = as.character(chromosome),
    gene_ensembl = strip_ensembl_version(target_gene_id)
  )]
}

standardise_eqtl <- function(dt, gene_id) {
  if (is.null(dt) || nrow(dt) == 0L) return(empty_dt())
  dt <- data.table::as.data.table(dt)

  aliases <- list(
    SNP = c("SNP", "rsid", "snp", "variant_id"),
    pval = c("pval", "pvalue", "p_value"),
    eaf = c("eaf", "maf", "allele_frequency"),
    effect_allele = c("effect_allele", "alt", "assessed_allele"),
    other_allele = c("other_allele", "ref"),
    pos = c("pos", "position", "variant_position")
  )

  for (target in names(aliases)) {
    source <- aliases[[target]][aliases[[target]] %in% names(dt)][1]
    if (!is.na(source) && source != target && !target %in% names(dt)) {
      data.table::setnames(dt, source, target)
    }
  }

  required <- c(
    "SNP", "beta", "se", "pval", "eaf",
    "effect_allele", "other_allele", "pos"
  )
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    warning("Skipping ", gene_id, ": eQTL columns missing: ",
            paste(missing, collapse = ", "))
    return(empty_dt())
  }

  dt[, gene_ensembl := strip_ensembl_version(gene_id)]
  dt[, `:=`(
    SNP = as.character(SNP),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pval = as.numeric(pval),
    eaf = as.numeric(eaf),
    maf = pmin(as.numeric(eaf), 1 - as.numeric(eaf)),
    pos = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele))
  )]

  dt <- dt[
    grepl("^rs", SNP) & is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(pval) & is.finite(eaf) & eaf > 0 & eaf < 1
  ]
  deduplicate_variant_rows(dt, source = gene_id, scope = "GTEx_v8_eQTL_region")
}

load_gene_eqtl <- function(gene_id, chromosome, start, end) {
  cache_file <- file.path(config$paths$cache, paste0(gene_id, "_muscle.rds"))

  if (file.exists(cache_file)) {
    return(standardise_eqtl(readRDS(cache_file), gene_id))
  }

  raw <- read_local_eqtl_region(gene_id, chromosome, start, end)

  eqtl <- standardise_eqtl(raw, gene_id)
  saveRDS(eqtl, cache_file)
  eqtl
}

format_exposure <- function(
  eqtl, gene_symbol,
  exposure_n = config$gtex$sample_size
) {
  if (nrow(eqtl) == 0L) return(NULL)
  eqtl <- copy(eqtl)
  eqtl[, F_stat := (beta / se)^2]
  eqtl <- eqtl[
    pval < config$instruments$p_threshold &
      F_stat >= config$instruments$minimum_F
  ]
  if (nrow(eqtl) == 0L) return(NULL)

  out <- tryCatch(
    TwoSampleMR::format_data(
      as.data.frame(eqtl),
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      chr_col = "chr"
    ),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out$exposure <- gene_symbol
  out$id.exposure <- gene_symbol
  out$samplesize.exposure <- exposure_n
  out[!duplicated(out$SNP), ]
}

clump_exposure <- function(exposure) {
  if (is.null(exposure) || nrow(exposure) <= 1L) return(exposure)

  if (config$instruments$clump_mode == "local") {
    chromosome <- unique(as.character(exposure$chr.exposure))
    if (length(chromosome) != 1L || !chromosome %in% as.character(1:22)) {
      warning("Skipping exposure: local LD reference covers autosomes 1-22 only.")
      return(NULL)
    }
    bfile <- paste0(config$instruments$ld_reference, ".", chromosome)
    if (!file.exists(config$instruments$plink_binary) ||
        !file.exists(paste0(bfile, ".bed"))) {
      stop("Local clumping requires valid plink_binary and ld_reference config.")
    }
    clumped <- tryCatch(
      ieugwasr::ld_clump(
        data.frame(
          rsid = exposure$SNP,
          pval = exposure$pval.exposure,
          id = exposure$id.exposure
        ),
        clump_kb = config$instruments$clump_kb,
        clump_r2 = config$instruments$clump_r2,
        clump_p = config$instruments$p_threshold,
        bfile = bfile,
        plink_bin = config$instruments$plink_binary
      ),
      error = function(e) {
        warning(
          "Skipping exposure: no usable variant remained after local LD clumping. ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(clumped) || nrow(clumped) == 0L) return(NULL)
    return(exposure[exposure$SNP %in% clumped$rsid, , drop = FALSE])
  }

  tryCatch(
    TwoSampleMR::clump_data(
      exposure,
      clump_kb = config$instruments$clump_kb,
      clump_r2 = config$instruments$clump_r2,
      clump_p1 = config$instruments$p_threshold,
      pop = config$instruments$ancestry
    ),
    error = function(e) {
      warning("LD clumping failed: ", conditionMessage(e),
              ". The gene is skipped; no top-SNP fallback was used.")
      NULL
    }
  )
}

format_outcome <- function(gwas, snps, trait) {
  sub <- gwas[SNP %in% snps]
  if (nrow(sub) == 0L) return(NULL)
  out <- tryCatch(
    TwoSampleMR::format_data(
      as.data.frame(sub),
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval"
    ),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out$outcome <- trait
  out$id.outcome <- trait
  out$samplesize.outcome <- config$outcome_n[[trait]]
  out
}


## 04. Candidate IR gene set -------------------------------------------------

assert_files_exist(
  c(config$files$gtex_muscle, paste0(config$files$gtex_muscle, ".tbi")),
  "GTEx skeletal-muscle eQTL"
)

# [REVISED: Reviewer 3] Fixed membership embedded for single-file distribution.
# The following rules document the original selection algorithm; reruns use the
# archived members below, rather than querying the current MSigDB release.
##   canonical_names <- c(
##     "KEGG_INSULIN_SIGNALING_PATHWAY",
##     "KEGG_TYPE_II_DIABETES_MELLITUS",
##     "KEGG_PPAR_SIGNALING_PATHWAY",
##     "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE",
##     "REACTOME_SIGNALLING_BY_INSULIN_RECEPTOR",
##     "REACTOME_PI3K_CASCADE",
##     "REACTOME_FOXO_MEDIATED_TRANSCRIPTION",
##     "REACTOME_MTOR_SIGNALLING",
##     "REACTOME_MTORC1_MEDIATED_SIGNALLING",
##     "GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY",
##     "GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS",
##     "GOBP_RESPONSE_TO_INSULIN_STIMULUS",
##     "GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS",
##     "GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS",
##     "GOBP_GLUCOSE_HOMEOSTASIS",
##     "GOBP_PI3K_SIGNALING")
##   
## 
##   expansion_pattern <- paste(
##     c(
##       "INSULIN_SIGNALING", "INSULIN_RECEPTOR", "INSULIN_RESISTANCE",
##       "GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN", "GLUCOSE_HOMEOSTASIS",
##       "SLC2A4", "GLUT4", "IRS1", "IRS2", "FOXO_MEDIATED", "MTORC1"
##     ),
##     collapse = "|")
##   
build_candidate_gene_set <- function() {
  pathways <- data.table::fread(text = "gs_name,gene_symbol,layer
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ADIPOQ,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ADIPOR1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,AHSG,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,AKT1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,AKT2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ALAS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,AP3S1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,APC,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,APPL1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,APPL2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ATP2B1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,BAIAP2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,BCAR1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,BCAR3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,BGLAP,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,BLVRB,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,C1QTNF12,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,C2CD5,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CAPN10,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CAV2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CEACAM1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CFLAR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,COL6A1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CPEB1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CPEB2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CSRP3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CTSD,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CUL3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,CUL7,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,DENND4C,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,DNAI1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ECHDC3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,EEF2K,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,EIF4EBP2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ENPP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,EPRS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ERFE,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ERRFI1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FBP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FBXW8,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FER,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FFAR3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FOXC2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FOXO1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FOXO4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,FUT7,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,G6PC1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GAB1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GCK,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GCLC,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GHR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GHRHR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GHSR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GKAP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GOT1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GPLD1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GPR21,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GRB10,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GRB14,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GRB2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GRB7,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GSK3A,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GSK3B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,GSTP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,HDAC5,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,HDAC9,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,HMGCS2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,HRAS,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IDE,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IGF1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IGF1R,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IGF2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IGFBP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IL1B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INHBB,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INPP5K,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INS,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INSIG1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INSIG2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INSR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,INSRR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IRS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IRS2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,IRS4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,KANK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,KAT2B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,KBTBD2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,LEP,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,LONP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,LPIN1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,LPIN2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,LPIN3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MAPK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MAPK3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MAPKAP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MIR103A1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MIR107,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MIR1271,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MIR15B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MIR195,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MSTN,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MTCL2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MTOR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MYO1C,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MYO5A,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,MZB1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NAMPT,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NCK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NCL,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NCOA1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NCOA2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NCOA5,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NDEL1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NR1H4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,NUCKS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,OGT,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,OSBPL8,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,OTOP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PAK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PARP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PCK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PCK2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PCSK9,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PDE3B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PDK2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PDK4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PDPK1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PHIP,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PID1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIK3C2A,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIK3CA,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIK3R1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIK3R2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIK3R3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIP4K2A,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIP4K2B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PIP4K2C,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PKLR,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PKM,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PLA2G1B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PNPLA3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,POU4F2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PPARG,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKAA1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKCB,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKCD,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKCI,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKCQ,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKCZ,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PRKDC,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PTPN1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PTPN11,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PTPN2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PTPRE,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,PTPRJ,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RAB10,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RAB13,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RAB31,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RAB8A,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RAF1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RARRES2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RB1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RBM4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RBX1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RELA,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RHOQ,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RPS6KB1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,RPS6KB2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SELENOS,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SERPINA12,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SESN3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SGCB,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SH2B2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SHC1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SIK2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SIRT1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC22A12,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC25A33,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC27A1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC27A4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC2A4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC2A8,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC39A14,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SLC9A1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SMARCC1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SNX5,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SOCS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SOCS3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SOCS7,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SORBS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SORL1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SOS1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SOS2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SRD5A1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SREBF1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SRSF5,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,STAT1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,STXBP4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,SYAP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TBC1D4,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TNS2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TRARG1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TRIB3,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TRIM72,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,TSC2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,UCP2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,USF1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,USO1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,VAMP2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,VWA2,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,WDTC1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,XBP1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,YWHAG,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ZBTB7B,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ZDHHC7,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ZFP36L1,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ZNF106,canonical
GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,ZNF592,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,SGCB,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,SLC27A1,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,SLC27A4,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,SLC2A4,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,TRARG1,canonical
GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS,ZDHHC7,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ADIPOR1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,AHSG,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,AKT1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,AKT2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,AP3S1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,APC,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,APPL1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,BAIAP2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,BCAR1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,BCAR3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,BLVRB,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,C1QTNF12,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,C2CD5,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CAV2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,COL6A1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CSRP3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CTSD,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CUL7,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,DNAI1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,EIF4EBP2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ENPP1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ERFE,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FBXW8,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FER,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FFAR3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FOXC2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FOXO1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FOXO4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FUT7,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GAB1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GKAP1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GPLD1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GPR21,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB10,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB14,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB7,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GSK3A,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GSK3B,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,HRAS,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IDE,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IGF1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IGF1R,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IGF2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IGFBP1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IL1B,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INPP5K,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INS,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INSR,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INSRR,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IRS1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IRS2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IRS4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,KANK1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,LEP,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,LONP1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MAPK1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MAPK3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MAPKAP1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR103A1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR107,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR1271,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR15B,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR195,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MSTN,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MTCL2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MTOR,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MZB1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NAMPT,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCK1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCL,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCOA5,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NDEL1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NR1H4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NUCKS1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,OGT,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,OSBPL8,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PAK1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PDK2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PDK4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PDPK1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PHIP,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PID1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIK3C2A,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIK3CA,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIK3R1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIK3R2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIK3R3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2A,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2B,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2C,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKAA1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCB,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCD,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCQ,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCZ,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN11,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPRE,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPRJ,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RAF1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RARRES2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RBM4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RBX1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RELA,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RHOQ,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RPS6KB1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RPS6KB2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SERPINA12,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SESN3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SH2B2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SHC1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SIK2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SIRT1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SLC27A4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SLC2A8,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SLC39A14,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SMARCC1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SNX5,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOCS1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOCS3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOCS7,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SORBS1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SORL1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOS1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOS2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SREBF1,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,STXBP4,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TNS2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TRIB3,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TRIM72,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TSC2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,VWA2,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ZBTB7B,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ZNF106,canonical
GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ZNF592,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ACACB,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ACADM,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ACTN3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ADCY10,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ADIPOQ,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ADIPOR1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,AKT1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,AKT2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ARPP19,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,BCKDK,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,C1QTNF1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,C1QTNF12,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,C1QTNF3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,CLK2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,CRY1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,DDB1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,DGAT2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,DGKQ,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,DYRK2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,EP300,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,EPM2AIP1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ERFE,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,FBP1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,FOXK1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,FOXK2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,FOXO1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GCG,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GCK,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GNB3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GNMT,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GPLD1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,GSK3A,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,HMGB1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IGF1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IGF2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IGFBP3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IGFBP4,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,INS,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,INSR,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IRS1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,IRS2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,KAT2A,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,KAT2B,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,LCMT1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,LEP,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,LEPR,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MIR103A1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MIR107,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MIR210,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MLYCD,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MST1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,MTCL2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NCOA2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NFE2L1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NLN,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NNMT,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NR0B1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,NR3C1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,OGT,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PDK1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PDK2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PDK3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PDK4,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PGP,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PHKA1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PHKG2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PMAIP1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPARA,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPARGC1A,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP1CA,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP1R3B,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP1R3E,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP1R3G,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP4R3A,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PPP4R3B,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PRKACA,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PRKAG1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PRKAG2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PRKAG3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PRKN,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PTH,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,PTPN2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,RANBP2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,RORA,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,RORC,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SDHAF3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SELENOS,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SERPINA12,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SESN2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SIK1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SIRT1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SIRT6,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SIRT7,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SLC35B4,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SLC45A3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,SORBS1,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,TCF7L2,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,TFF3,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,TIGAR,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,TP53,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,USP7,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,WDR5,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ZFP92,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ZMPSTE24,canonical
GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS,ZNF692,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,ACACA,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,ACACB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,AKT1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,AKT2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,AKT3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,ARAF,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,BAD,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,BRAF,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALM1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALM2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALM3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALML3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALML5,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CALML6,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CBL,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CBLB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CBLC,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CRK,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,CRKL,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,EIF4E,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,EIF4E1B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,EIF4E2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,EIF4EBP1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,ELK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,EXOC7,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FASN,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FBP1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FBP2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FLOT1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FLOT2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,FOXO1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,G6PC1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,G6PC2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,GCK,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,GRB2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,GSK3B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,GYS1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,GYS2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,HK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,HK2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,HK3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,HRAS,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,IKBKB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,INPP5D,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,INPP5K,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,INS,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,INSR,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,IRS1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,IRS2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,IRS4,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,KRAS,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,LIPE,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAP2K1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAP2K2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAPK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAPK10,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAPK3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAPK8,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MAPK9,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MKNK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MKNK2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,MTOR,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,NRAS,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PCK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PCK2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PDE3A,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PDE3B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PDPK1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PHKA1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PHKA2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PHKB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PHKG1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PHKG2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3CA,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3CB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3CD,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3CG,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3R1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3R2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3R3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PIK3R5,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PKLR,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPARGC1A,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1CA,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1CB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1CC,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1R3A,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1R3B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1R3C,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PPP1R3D,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAA1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAA2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAB1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAB2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKACA,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKACB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKACG,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAG1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAG2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAG3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAR1A,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAR1B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAR2A,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKAR2B,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKCI,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKCZ,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PRKX,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PTPN1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PTPRF,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PYGB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PYGL,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,PYGM,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RAF1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RAPGEF1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RHEB,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RHOQ,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RPS6,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RPS6KB1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RPS6KB2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,RPTOR,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SH2B2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SHC1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SHC2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SHC3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SHC4,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SLC2A4,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOCS1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOCS2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOCS3,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOCS4,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SORBS1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOS1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SOS2,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,SREBF1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,TRIP10,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,TSC1,canonical
KEGG_INSULIN_SIGNALING_PATHWAY,TSC2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACAA1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACADL,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACADM,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACOX1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACOX2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACOX3,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACSL1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACSL3,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACSL4,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACSL5,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ACSL6,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ADIPOQ,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ANGPTL4,canonical
KEGG_PPAR_SIGNALING_PATHWAY,APOA1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,APOA2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,APOA5,canonical
KEGG_PPAR_SIGNALING_PATHWAY,APOC3,canonical
KEGG_PPAR_SIGNALING_PATHWAY,AQP7,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CD36,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CPT1A,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CPT1B,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CPT1C,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CPT2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CYP27A1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CYP4A11,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CYP4A22,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CYP7A1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,CYP8B1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,DBI,canonical
KEGG_PPAR_SIGNALING_PATHWAY,EHHADH,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP3,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP4,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP5,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP6,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FABP7,canonical
KEGG_PPAR_SIGNALING_PATHWAY,FADS2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,GK,canonical
KEGG_PPAR_SIGNALING_PATHWAY,GK2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,HMGCS2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ILK,canonical
KEGG_PPAR_SIGNALING_PATHWAY,LPL,canonical
KEGG_PPAR_SIGNALING_PATHWAY,ME1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,MMP1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,NR1H3,canonical
KEGG_PPAR_SIGNALING_PATHWAY,OLR1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PCK1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PCK2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PDPK1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PLIN1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PLTP,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PPARA,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PPARD,canonical
KEGG_PPAR_SIGNALING_PATHWAY,PPARG,canonical
KEGG_PPAR_SIGNALING_PATHWAY,RXRA,canonical
KEGG_PPAR_SIGNALING_PATHWAY,RXRB,canonical
KEGG_PPAR_SIGNALING_PATHWAY,RXRG,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SCD,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SCD5,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SCP2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SLC27A1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SLC27A2,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SLC27A4,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SLC27A5,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SLC27A6,canonical
KEGG_PPAR_SIGNALING_PATHWAY,SORBS1,canonical
KEGG_PPAR_SIGNALING_PATHWAY,UBC,canonical
KEGG_PPAR_SIGNALING_PATHWAY,UCP1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,ABCC8,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,ADIPOQ,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1A,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1B,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1C,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1D,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1E,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,CACNA1G,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,GCK,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,HK1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,HK2,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,HK3,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,IKBKB,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,INS,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,INSR,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,IRS1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,IRS2,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,IRS4,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,KCNJ11,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAFA,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAPK1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAPK10,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAPK3,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAPK8,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MAPK9,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,MTOR,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PDX1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3CA,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3CB,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3CD,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3CG,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3R1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3R2,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3R3,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PIK3R5,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PKLR,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PKM,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PRKCD,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PRKCE,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,PRKCZ,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SLC2A2,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SLC2A4,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SOCS1,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SOCS2,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SOCS3,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,SOCS4,canonical
KEGG_TYPE_II_DIABETES_MELLITUS,TNF,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,ABCA6,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,AGRP,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,AKT1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,AKT2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,AKT3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,ATXN3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,BBC3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,BCL2L11,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,BCL6,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,BTG1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CAT,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CAV1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CCNG2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CDKN1A,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CDKN1B,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CITED2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,CREBBP,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,DDIT3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,EP300,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FASLG,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FBXO32,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FOXG1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FOXO1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FOXO3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FOXO4,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,FOXO6,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,G6PC1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,GADD45A,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,GCK,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,HDAC1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,HDAC2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,IGFBP1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,INS,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,KAT2B,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,KLF4,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,MSTN,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,NFYA,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,NFYB,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,NFYC,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,NPY,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,NR3C1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,PCBP4,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,PCK1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,PINK1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,PLXNA4,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,POMC,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,PPARGC1A,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,RBL2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,RETN,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SFN,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SIN3A,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SIRT1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SIRT3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SMAD2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SMAD3,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SMAD4,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SOD2,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,SREBF1,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,STK11,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,TRIM63,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,TXN,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,TXNIP,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,YWHAB,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,YWHAG,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,YWHAQ,canonical
REACTOME_FOXO_MEDIATED_TRANSCRIPTION,YWHAZ,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,AKT2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF10,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF16,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF17,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF18,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF19,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF20,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF22,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF23,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF3,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF4,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF5,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF6,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF7,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF8,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGF9,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGFR1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGFR2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGFR3,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FGFR4,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FLT3,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FLT3LG,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,FRS2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,GAB1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,GAB2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,GRB10,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,GRB2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,HRAS,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,INS,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,INSR,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,IRS1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,IRS2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,KL,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,KLB,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,KRAS,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,MAPK1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,MAPK3,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,NRAS,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PDE3B,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PDPK1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3C3,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3CA,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3CB,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3R1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3R2,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PIK3R4,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,PTPN11,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,SHC1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,SOS1,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,THEM4,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,TLR9,canonical
REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE,TRIB3,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,AKT1S1,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,EEF2K,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,EIF4B,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,EIF4E,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,EIF4EBP1,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,EIF4G1,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,FKBP1A,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,LAMTOR1,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,LAMTOR2,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,LAMTOR3,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,LAMTOR4,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,LAMTOR5,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,MLST8,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,MTOR,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RHEB,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RPS6,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RPS6KB1,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RPTOR,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RRAGA,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RRAGB,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RRAGC,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,RRAGD,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,SLC38A9,canonical
REACTOME_MTORC1_MEDIATED_SIGNALLING,YWHAB,canonical
REACTOME_MTOR_SIGNALLING,AKT1,canonical
REACTOME_MTOR_SIGNALLING,AKT1S1,canonical
REACTOME_MTOR_SIGNALLING,AKT2,canonical
REACTOME_MTOR_SIGNALLING,CAB39,canonical
REACTOME_MTOR_SIGNALLING,CAB39L,canonical
REACTOME_MTOR_SIGNALLING,EEF2K,canonical
REACTOME_MTOR_SIGNALLING,EIF4B,canonical
REACTOME_MTOR_SIGNALLING,EIF4E,canonical
REACTOME_MTOR_SIGNALLING,EIF4EBP1,canonical
REACTOME_MTOR_SIGNALLING,EIF4G1,canonical
REACTOME_MTOR_SIGNALLING,FKBP1A,canonical
REACTOME_MTOR_SIGNALLING,LAMTOR1,canonical
REACTOME_MTOR_SIGNALLING,LAMTOR2,canonical
REACTOME_MTOR_SIGNALLING,LAMTOR3,canonical
REACTOME_MTOR_SIGNALLING,LAMTOR4,canonical
REACTOME_MTOR_SIGNALLING,LAMTOR5,canonical
REACTOME_MTOR_SIGNALLING,MLST8,canonical
REACTOME_MTOR_SIGNALLING,MTOR,canonical
REACTOME_MTOR_SIGNALLING,PPM1A,canonical
REACTOME_MTOR_SIGNALLING,PRKAA1,canonical
REACTOME_MTOR_SIGNALLING,PRKAA2,canonical
REACTOME_MTOR_SIGNALLING,PRKAB1,canonical
REACTOME_MTOR_SIGNALLING,PRKAB2,canonical
REACTOME_MTOR_SIGNALLING,PRKAG1,canonical
REACTOME_MTOR_SIGNALLING,PRKAG2,canonical
REACTOME_MTOR_SIGNALLING,PRKAG3,canonical
REACTOME_MTOR_SIGNALLING,RHEB,canonical
REACTOME_MTOR_SIGNALLING,RPS6,canonical
REACTOME_MTOR_SIGNALLING,RPS6KB1,canonical
REACTOME_MTOR_SIGNALLING,RPTOR,canonical
REACTOME_MTOR_SIGNALLING,RRAGA,canonical
REACTOME_MTOR_SIGNALLING,RRAGB,canonical
REACTOME_MTOR_SIGNALLING,RRAGC,canonical
REACTOME_MTOR_SIGNALLING,RRAGD,canonical
REACTOME_MTOR_SIGNALLING,SLC38A9,canonical
REACTOME_MTOR_SIGNALLING,STK11,canonical
REACTOME_MTOR_SIGNALLING,STRADA,canonical
REACTOME_MTOR_SIGNALLING,STRADB,canonical
REACTOME_MTOR_SIGNALLING,TSC1,canonical
REACTOME_MTOR_SIGNALLING,TSC2,canonical
REACTOME_MTOR_SIGNALLING,YWHAB,canonical
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ABCA12,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ABCC8,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ADCY5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ADCY8,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ADRA2A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,AGER,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ANO1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,BAD,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,BAIAP3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,BRSK2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,C1QTNF12,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,C2CD2L,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CARTPT,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CASR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CCDC186,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CDK16,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CFTR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CLTRN,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CMA1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,COL1A1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CRH,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CYBA,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,CYP7A1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,DYNLL1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,EFNA5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ENDOG,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ENY2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,EPHA5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ERN1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FIS1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FKBP1B,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXA2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXA3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXK1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXK2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXO1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FOXO3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,FUT1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GAS6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GATA4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GCG,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GCK,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GCKR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GCLC,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GCLM,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GHRHR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GHRL,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GJB6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPER1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPLD1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPR27,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPR68,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPRC6A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,GPX1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HIF1A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HK1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HK2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HK3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HKDC1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,HLA-DRB1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ICAM1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,IGF1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,IGF1R,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,IRS2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,JAGN1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,KAT5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,KCNB1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,KCNK16,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,KLF15,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,KLF7,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,LIN28A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,LRP5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,LRRC8A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,LRRC8D,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR103A1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR146A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR15A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR16-1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320B1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320B2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320C1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320C2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320D1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320D2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR320E,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIR337,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MIRLET7G,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MLXIPL,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,MPC2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NADK,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NCF1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NDUFAF2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NGFR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NKX6-1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NOX4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NPTX1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NR1D1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NR1H4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,NUCKS1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,OGT,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,OPRK1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,OSBP,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,OXCT1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PAX2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PCK1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PCK2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PDE8B,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PDX1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PHPT1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PIH1D1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PIK3CA,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PIK3R1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PIK3R2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PIM3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PLA2G6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PPARD,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PPARGC1A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PPP3CB,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PRKAA1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PRKAA2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PRKACA,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PRKCE,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PRKN,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PTPRN,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,PTPRN2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RAB11B,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RAB11FIP2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RAB11FIP5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RAC1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RACK1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RAF1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RBM4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RFX6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,RPTOR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SELENOT,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SERPINF1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SIDT2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SIN3A,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SIRT1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC12A6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC12A7,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC29A1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC2A2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC39A14,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SLC9B2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SMAD3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SMAD4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SMARCA4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SMARCB1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SOX4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SRF,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SRI,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,STX4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,STXBP3,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,STXBP4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,SYBU,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TGFB1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TRA2B,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TREM2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TRPA1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TRPM4,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TRPM5,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,TUNAR,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,UBTF,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,UCP2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,UNC13B,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,USF1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,USF2,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,VSNL1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,XBP1,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ZBED6,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ZBTB20,expanded
GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS,ZNF236,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ADIPOR1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,AHSG,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,BLVRB,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,C1QTNF12,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CTSD,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,CUL7,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ENPP1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ERFE,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FBXW8,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FFAR3,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,FUT7,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GKAP1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GPR21,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB10,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB14,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GRB7,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,GSK3A,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IGF2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IL1B,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INPP5K,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,INS,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,IRS1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,KANK1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,LEP,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,LONP1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MAPKAP1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR103A1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR107,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR1271,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR15B,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MIR195,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MSTN,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MTOR,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,MZB1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCK1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCL,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NCOA5,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NR1H4,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,NUCKS1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,OGT,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,OSBPL8,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PAK1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PID1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2A,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2B,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PIP4K2C,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKAA1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCB,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCD,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCQ,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PRKCZ,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN11,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPN2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPRE,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,PTPRJ,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RBM4,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RBX1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RELA,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RPS6KB1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,RPS6KB2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SERPINA12,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SESN3,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SIK2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SIRT1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SLC27A4,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SNX5,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOCS1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SOCS3,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SORBS1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,SORL1,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TNS2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TRIB3,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TRIM72,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,TSC2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,VWA2,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ZBTB7B,expanded
GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,ZNF592,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,AKT1S1,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,ATG101,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,ATG13,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,ATG14,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,BECN1,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,DEPTOR,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,MLST8,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,MTOR,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,NRBF2,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,PIK3C3,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,PIK3R4,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,RB1CC1,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,RPTOR,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,TELO2,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,TTI1,expanded
KEGG_MEDICUS_REFERENCE_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION_MTORC1_PI3KC3_C1,ULK1,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,AKT1S1,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,DEPTOR,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,FLCN,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,MLST8,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,MTOR,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,RPTOR,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,RRAGA,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,RRAGB,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,RRAGC,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,RRAGD,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,SLC3A2,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,SLC7A5,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,TELO2,expanded
KEGG_MEDICUS_REFERENCE_FLCN_MTORC1_SIGNALING_PATHWAY,TTI1,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,AKT1S1,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,DEPDC5,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,DEPTOR,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,MLST8,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,MTOR,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,NPRL2,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,NPRL3,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,RPTOR,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,RRAGA,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,RRAGB,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,RRAGC,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,RRAGD,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,TELO2,expanded
KEGG_MEDICUS_REFERENCE_GATOR1_MTORC1_SIGNALING_PATHWAY,TTI1,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,AKT1S1,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,DEPTOR,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,MLST8,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,MTOR,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,RHEB,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,RPTOR,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,TBC1D7,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,TELO2,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,TSC1,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,TSC2,expanded
KEGG_MEDICUS_REFERENCE_TSC1_2_MTORC1_SIGNALING_PATHWAY,TTI1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0B,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0C,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0D1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0D2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0E1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V0E2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1A,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1B1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1B2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1C1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1C2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1D,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1E1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1E2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1F,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1G1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1G2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1G3,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ATP6V1H,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,CASTOR1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,CASTOR2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,DEPDC5,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,FLCN,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,FNIP1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,FNIP2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,ITFG2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,KICS2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,KPTN,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,LAMTOR1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,LAMTOR2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,LAMTOR3,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,LAMTOR4,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,LAMTOR5,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,MIOS,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,MLST8,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,MTOR,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,NPRL2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,NPRL3,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RHEB,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RPTOR,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RRAGA,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RRAGB,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RRAGC,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,RRAGD,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SAMTOR,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SEC13,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SEH1L,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SESN1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SESN2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SH3BP4,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SLC38A9,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,SZT2,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,TCIRG1,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,WDR24,expanded
REACTOME_AMINO_ACIDS_REGULATE_MTORC1,WDR59,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,BTG1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,CAV1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,CCNG2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,CDKN1A,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,CDKN1B,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,FOXG1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,FOXO1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,FOXO3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,FOXO4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,GADD45A,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,KLF4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,MSTN,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,PCBP4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,RBL2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,SMAD2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,SMAD3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES,SMAD4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,BBC3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,BCL2L11,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,BCL6,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,CITED2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,CREBBP,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,DDIT3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,EP300,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,FASLG,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,FOXO1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,FOXO3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,FOXO4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,NFYA,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,NFYB,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,NFYC,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,PINK1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES,STK11,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,ABCA6,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,AGRP,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,ATXN3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,CAT,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,FBXO32,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,FOXO1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,FOXO3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,FOXO4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,FOXO6,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,G6PC1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,GCK,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,HDAC1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,HDAC2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,IGFBP1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,INS,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,NPY,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,NR3C1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,PCK1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,PLXNA4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,POMC,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,PPARGC1A,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,RETN,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SIN3A,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SIRT3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SMAD2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SMAD3,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SMAD4,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SOD2,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,SREBF1,expanded
REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES,TRIM63,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6AP1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0A1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0A2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0A4,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0B,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0C,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0D1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0D2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0E1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V0E2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1A,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1B1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1B2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1C1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1C2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1D,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1E1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1E2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1F,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1G1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1G2,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1G3,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,ATP6V1H,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,CTSD,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,IDE,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,INS,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,INSR,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,PTPN1,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,PTPRF,expanded
REACTOME_INSULIN_RECEPTOR_RECYCLING,TCIRG1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,AKT2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6AP1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0A1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0A2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0A4,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0B,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0C,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0D1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0D2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0E1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V0E2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1A,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1B1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1B2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1C1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1C2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1D,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1E1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1E2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1F,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1G1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1G2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1G3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,ATP6V1H,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,CTSD,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF10,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF16,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF17,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF18,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF19,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF20,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF22,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF23,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF4,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF5,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF6,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF7,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF8,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGF9,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGFR1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGFR2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGFR3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FGFR4,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FLT3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FLT3LG,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,FRS2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,GAB1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,GAB2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,GRB10,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,GRB2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,HRAS,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,IDE,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,INS,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,INSR,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,IRS1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,IRS2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,KL,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,KLB,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,KRAS,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,MAPK1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,MAPK3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,NRAS,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PDE3B,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PDPK1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3C3,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3CA,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3CB,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3R1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3R2,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PIK3R4,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PTPN1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PTPN11,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,PTPRF,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,SHC1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,SOS1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,TCIRG1,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,THEM4,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,TLR9,expanded
REACTOME_SIGNALING_BY_INSULIN_RECEPTOR,TRIB3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,ACTB,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,ACTG1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,AKT1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,AKT2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,ASPSCR1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,C2CD5,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,CALM1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC4,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC5,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC6,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC7,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,EXOC8,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,KIF3A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,KIF3B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,KIFAP3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,LNPEP,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,MYH9,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,MYO1C,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,MYO5A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAA2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAB1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAB2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAG1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAG2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,PRKAG3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB10,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB11A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB13,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB14,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB4A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAB8A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RAC1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RALA,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RALGAPA2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RALGAPB,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,RHOQ,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,SFN,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,SLC2A4,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,SNAP23,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,STX4,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,STXBP3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TBC1D1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TBC1D4,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA1A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA1B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA1C,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA3C,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA3D,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA3E,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA4A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA4B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBA8,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBAL3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB1,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB2A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB2B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB3,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB4A,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB4B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB6,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB8,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,TUBB8B,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,VAMP2,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAB,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAE,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAG,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAH,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAQ,expanded
REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE,YWHAZ,expanded", sep = ",")
  stopifnot(nrow(pathways) == 1514L, uniqueN(pathways$gs_name) == 24L,
            uniqueN(pathways$gene_symbol) == 771L,
            all(pathways$layer %in% c("canonical", "expanded")))
  membership <- pathways |>
    dplyr::group_by(gene_symbol) |>
    dplyr::summarise(
      in_canonical = any(layer == "canonical"),
      in_expanded = any(layer == "expanded"),
      gene_layer = ifelse(in_canonical, "canonical", "expanded_only"),
      n_pathways = dplyr::n_distinct(gs_name), .groups = "drop")
  list(pathways = pathways, membership = membership)
}

gene_set <- build_candidate_gene_set()
candidate_genes <- gene_set$membership$gene_symbol

core_genes <- c(
  "INSR", "IRS1", "IRS2", "PIK3CA", "PIK3R1", "AKT1", "AKT2",
  "SLC2A4", "FOXO1", "FOXO3", "MTOR", "RPTOR", "RPS6KB1",
  "PRKAA1", "PRKAA2", "PPARG", "FBXO32", "TRIM63" )


gene_set_qc <- data.table(
  metric = c(
    "candidate_genes", "canonical_genes", "expanded_only_genes",
    "pathways", "core_genes_recovered", "core_gene_recovery_rate"
  ),
  value = c(
    length(candidate_genes),
    sum(gene_set$membership$gene_layer == "canonical"),
    sum(gene_set$membership$gene_layer == "expanded_only"),
    dplyr::n_distinct(gene_set$pathways$gs_name),
    sum(core_genes %in% candidate_genes),
    mean(core_genes %in% candidate_genes)))


write_result(gene_set$pathways, "01_candidate_pathway_membership.csv")
write_result(gene_set$membership, "01_candidate_gene_layers.csv")
write_result(gene_set_qc, "01_candidate_gene_set_qc.csv")

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = candidate_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
) |>
  dplyr::filter(!is.na(ENSEMBL)) |>
  dplyr::mutate(ENSEMBL = strip_ensembl_version(ENSEMBL)) |>
  dplyr::distinct(SYMBOL, ENSEMBL)

gene_ranges <- ensembldb::genes(
  EnsDb.Hsapiens.v86,
  filter = AnnotationFilter::GeneIdFilter(unique(gene_map$ENSEMBL)),
  return.type = "GRanges")


gene_positions <- as.data.frame(gene_ranges) |>
  dplyr::transmute(
    gene_ensembl = strip_ensembl_version(gene_id),
    chr = as.character(seqnames),
    start = as.integer(start),
    end = as.integer(end),
    midpoint = as.integer(round((start + end) / 2))
  ) |>
  dplyr::inner_join(
    gene_map |> dplyr::rename(gene_symbol = SYMBOL, gene_ensembl = ENSEMBL),
    by = "gene_ensembl"
  ) |>
  dplyr::left_join(gene_set$membership, by = "gene_symbol") |>
  dplyr::filter(chr %in% c(as.character(1:22), "X")) |>
  dplyr::distinct(gene_ensembl, .keep_all = TRUE)

write_result(gene_positions, "01_candidate_gene_coordinates.csv")
log_step("Candidate genes: ", length(candidate_genes),
         "; mapped Ensembl genes: ", nrow(gene_positions))


## 05. Read outcome GWAS ------------------------------------------------------

assert_files_exist(config$files$outcomes, "Outcome GWAS")

gwas_list <- lapply(names(config$files$outcomes), function(trait) {
  vcf_to_gwas(
    config$files$outcomes[[trait]],
    trait,
    config$outcome_n[[trait]]
  )
})
names(gwas_list) <- names(config$files$outcomes)


## 06. Skeletal-muscle cis-MR -------------------------------------------------

run_cis_mr <- function(gene_positions, gwas_list) {
  mr_rows <- list()
  harmonised <- list()
  attrition <- list()
  pb <- utils::txtProgressBar(min = 0, max = nrow(gene_positions), style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(nrow(gene_positions))) {
    utils::setTxtProgressBar(pb, i)
    gp <- gene_positions[i, ]

    eqtl <- load_gene_eqtl(
      gp$gene_ensembl, gp$chr, gp$start, gp$end
    )
    n_region <- nrow(eqtl)
    exposure <- format_exposure(eqtl, gp$gene_symbol)
    n_strong <- if (is.null(exposure)) 0L else nrow(exposure)
    exposure <- clump_exposure(exposure)
    n_clumped <- if (is.null(exposure)) 0L else nrow(exposure)

    attrition[[gp$gene_ensembl]] <- data.table(
      gene_ensembl = gp$gene_ensembl,
      gene_symbol = gp$gene_symbol,
      gene_layer = gp$gene_layer,
      n_eqtl_region = n_region,
      n_strong_F_eligible = n_strong,
      n_after_clumping = n_clumped
    )
    if (is.null(exposure) || nrow(exposure) == 0L) next

    for (trait in names(gwas_list)) {
      outcome <- format_outcome(gwas_list[[trait]], exposure$SNP, trait)
      if (is.null(outcome)) next

      ## action = 3 is deliberately conservative because some eQTL endpoints
      ## return minor-allele frequency rather than effect-allele frequency.
      dat <- tryCatch(
        TwoSampleMR::harmonise_data(exposure, outcome, action = 3),
        error = function(e) NULL
      )
      if (is.null(dat)) next
      dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
      if (nrow(dat) == 0L) next

      method <- if (nrow(dat) == 1L) "mr_wald_ratio" else "mr_ivw"
      estimate <- tryCatch(
        TwoSampleMR::mr(dat, method_list = method),
        error = function(e) NULL
      )
      if (is.null(estimate) || nrow(estimate) == 0L) next

      key <- paste(gp$gene_ensembl, trait, sep = "__")
      harmonised[[key]] <- dat
      mr_rows[[key]] <- data.table(
        gene_ensembl = gp$gene_ensembl,
        gene_symbol = gp$gene_symbol,
        gene_layer = gp$gene_layer,
        trait = trait,
        method = estimate$method[1],
        nsnp = estimate$nsnp[1],
        beta = estimate$b[1],
        se = estimate$se[1],
        pval = estimate$pval[1],
        min_F = min((dat$beta.exposure / dat$se.exposure)^2, na.rm = TRUE),
        mean_F = mean((dat$beta.exposure / dat$se.exposure)^2, na.rm = TRUE)
      )
    }
  }

  list(
    mr = bind_nonempty(mr_rows),
    harmonised = harmonised,
    attrition = bind_nonempty(attrition)
  )
}

log_step("Running skeletal-muscle cis-MR")
cis_mr <- run_cis_mr(gene_positions, gwas_list)
mr_main <- cis_mr$mr

# [REVISED: Reviewer 1 comment 6] Preserve source-level duplicate counts and
# the exact deterministic rule used before MR/colocalization joins.
write_result(variant_dedup_audit, "06_variant_dedup_audit.csv")

if (nrow(mr_main) == 0L) stop("No analysable cis-MR results were produced.")

mr_main[, fdr_within_trait := p.adjust(pval, method = "BH"), by = trait]
mr_main[, fdr_within_trait_layer := p.adjust(pval, method = "BH"),
        by = .(trait, gene_layer)]
mr_main[, fdr_global := p.adjust(pval, method = "BH")]
mr_main[, ci_lower := beta - 1.96 * se]
mr_main[, ci_upper := beta + 1.96 * se]
mr_main[, mr_associated := fdr_within_trait < config$multiple_testing$alpha]
mr_main[, mr_associated_layered :=
          fdr_within_trait_layer < config$multiple_testing$alpha]

write_result(mr_main, "02_cisMR_all_results.csv")
write_result(mr_main[mr_associated == TRUE], "02_cisMR_FDR_significant.csv")
write_result(cis_mr$attrition, "02_cisMR_gene_attrition.csv")


## 07. MR sensitivity analyses ------------------------------------------------

run_mr_sensitivity <- function(mr_main, harmonised) {
  sensitivity <- list()
  heterogeneity <- list()
  pleiotropy <- list()
  steiger <- list()

  significant_pairs <- mr_main[mr_associated == TRUE]
  for (i in seq_len(nrow(significant_pairs))) {
    row <- significant_pairs[i]
    key <- paste(row$gene_ensembl, row$trait, sep = "__")
    dat <- harmonised[[key]]
    if (is.null(dat) || nrow(dat) == 0L) next

    direction <- tryCatch(
      TwoSampleMR::steiger_filtering(dat),
      error = function(e) NULL
    )
    if (!is.null(direction)) {
      steiger[[key]] <- data.table(
        gene_ensembl = row$gene_ensembl,
        gene_symbol = row$gene_symbol,
        trait = row$trait,
        nsnp = nrow(dat),
        n_direction_correct = sum(direction$steiger_dir %in% TRUE, na.rm = TRUE),
        proportion_direction_correct = mean(direction$steiger_dir %in% TRUE, na.rm = TRUE)
      )
    }

    if (nrow(dat) >= 2L) {
      h <- tryCatch(TwoSampleMR::mr_heterogeneity(dat), error = function(e) NULL)
      if (!is.null(h)) {
        h$gene_ensembl <- row$gene_ensembl
        h$gene_symbol <- row$gene_symbol
        h$trait <- row$trait
        heterogeneity[[key]] <- h
      }
    }

    if (nrow(dat) >= 3L) {
      s <- tryCatch(
        TwoSampleMR::mr(
          dat,
          method_list = c(
            "mr_ivw", "mr_egger_regression", "mr_weighted_median",
            "mr_weighted_mode"
          )
        ),
        error = function(e) NULL
      )
      if (!is.null(s)) {
        s$gene_ensembl <- row$gene_ensembl
        s$gene_symbol <- row$gene_symbol
        s$trait <- row$trait
        sensitivity[[key]] <- s
      }

      p <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL)
      if (!is.null(p)) {
        p$gene_ensembl <- row$gene_ensembl
        p$gene_symbol <- row$gene_symbol
        p$trait <- row$trait
        pleiotropy[[key]] <- p
      }
    }
  }

  list(
    sensitivity = bind_nonempty(sensitivity),
    heterogeneity = bind_nonempty(heterogeneity),
    pleiotropy = bind_nonempty(pleiotropy),
    steiger = bind_nonempty(steiger)
  )
}

log_step("Running sensitivity analyses for FDR-significant cis-MR pairs")
mr_sensitivity <- run_mr_sensitivity(mr_main, cis_mr$harmonised)
write_result(mr_sensitivity$sensitivity, "03_cisMR_sensitivity_multiSNP.csv")
write_result(mr_sensitivity$heterogeneity, "03_cisMR_heterogeneity_same_exposure.csv")
write_result(mr_sensitivity$pleiotropy, "03_cisMR_Egger_intercept.csv")
write_result(mr_sensitivity$steiger, "03_cisMR_Steiger.csv")


## 08. Colocalisation and prior sensitivity ----------------------------------

harmonise_coloc_region <- function(eqtl, gwas) {
  common <- intersect(eqtl$SNP, gwas$SNP)
  if (length(common) < config$coloc$minimum_common_snps) return(NULL)

  e <- eqtl[match(common, SNP)]
  g <- gwas[match(common, SNP)]

  keep <- is.finite(e$beta) & is.finite(e$se) & e$se > 0 &
    is.finite(e$maf) & e$maf >= config$coloc$minimum_maf & e$maf < 0.5 &
    is.finite(g$beta) & is.finite(g$se) & g$se > 0 &
    is.finite(g$maf) & g$maf >= config$coloc$minimum_maf & g$maf < 0.5
  e <- e[keep]
  g <- g[keep]
  if (nrow(e) < config$coloc$minimum_common_snps) return(NULL)

  same <- e$effect_allele == g$effect_allele &
    e$other_allele == g$other_allele
  flipped <- e$effect_allele == g$other_allele &
    e$other_allele == g$effect_allele
  ambiguous <- is_palindromic(e$effect_allele, e$other_allele) &
    (e$maf > 0.42 | g$maf > 0.42)
  keep <- (same | flipped) & !ambiguous
  e <- e[keep]
  g <- g[keep]
  flipped <- flipped[keep]
  if (any(flipped)) {
    old_effect <- e$effect_allele[flipped]
    e[flipped, `:=`(
      beta = -beta,
      eaf = 1 - eaf,
      effect_allele = other_allele,
      other_allele = old_effect
    )]
  }

  if (nrow(e) < config$coloc$minimum_common_snps) return(NULL)
  list(eqtl = e, gwas = g)
}

run_coloc_pair <- function(
  eqtl, gwas, gene_symbol, trait,
  exposure_n = config$gtex$sample_size
) {
  aligned <- harmonise_coloc_region(eqtl, gwas)
  if (is.null(aligned)) return(empty_dt())
  e <- aligned$eqtl
  g <- aligned$gwas

  d_eqtl <- list(
    snp = e$SNP,
    beta = e$beta,
    varbeta = e$se^2,
    MAF = e$maf,
    N = exposure_n,
    type = "quant",
    sdY = 1
  )
  d_gwas <- list(
    snp = g$SNP,
    beta = g$beta,
    varbeta = g$se^2,
    MAF = g$maf,
    N = stats::median(g$N, na.rm = TRUE)
  )
  if (config$outcome_type[[trait]] == "cc") {
    d_gwas$type <- "cc"
    d_gwas$s <- config$outcome_cases[[trait]] / config$outcome_n[[trait]]
  } else {
    sdY_gwas <- sqrt(stats::median(
      2 * g$maf * (1 - g$maf) * g$N * g$se^2,
      na.rm = TRUE
    ))
    if (!is.finite(sdY_gwas) || sdY_gwas <= 0) return(empty_dt())
    d_gwas$type <- "quant"
    d_gwas$sdY <- sdY_gwas
  }

  rows <- lapply(config$coloc$p12_values, function(p12) {
    fit <- tryCatch(
      coloc::coloc.abf(
        dataset1 = d_eqtl,
        dataset2 = d_gwas,
        p1 = config$coloc$p1,
        p2 = config$coloc$p2,
        p12 = p12
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- fit$summary
    data.table(
      gene_symbol = gene_symbol,
      trait = trait,
      p1 = config$coloc$p1,
      p2 = config$coloc$p2,
      p12 = p12,
      n_snps = nrow(e),
      PP.H0 = unname(s["PP.H0.abf"]),
      PP.H1 = unname(s["PP.H1.abf"]),
      PP.H2 = unname(s["PP.H2.abf"]),
      PP.H3 = unname(s["PP.H3.abf"]),
      PP.H4 = unname(s["PP.H4.abf"])
    )
  })
  bind_nonempty(rows)
}

make_locus_plot <- function(eqtl, gwas, gene_symbol, trait) {
  common <- intersect(eqtl$SNP, gwas$SNP)
  if (length(common) == 0L) return(invisible(NULL))
  e <- eqtl[match(common, SNP)]
  g <- gwas[match(common, SNP)]
  plot_dt <- rbind(
    # [REVISED: Reviewer 6] Both panels use matched GRCh37 GWAS positions.
    data.table(pos = g$pos[match(e$SNP, g$SNP)], minus_log10_p = -log10(pmax(e$pval, 1e-300)),
               dataset = "Skeletal-muscle eQTL"),
    data.table(pos = g$pos, minus_log10_p = -log10(pmax(g$pval, 1e-300)),
               dataset = trait)
  )
  p <- ggplot(plot_dt, aes(pos / 1e6, minus_log10_p, colour = dataset)) +
    geom_point(alpha = 0.65, size = 1.3) +
    facet_wrap(~dataset, ncol = 1, scales = "free_y") +
    labs(
      title = paste(gene_symbol, trait, sep = " - "),
      x = "GRCh37 genomic position (Mb)", y = expression(-log[10](P)), colour = NULL
    ) +
    theme_manuscript()
  ggsave(
    plot_path(paste0("locus_", gene_symbol, "_", trait, ".png")),
    p, width = 8, height = 6, dpi = 300
  )
}

run_colocalisation <- function(mr_main, gene_positions, gwas_list) {
  targets <- mr_main[mr_associated == TRUE, .(
    gene_ensembl, gene_symbol, gene_layer, trait
  )]
  results <- list()

  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    gp <- gene_positions[
      gene_positions$gene_ensembl == target$gene_ensembl,
      ,
      drop = FALSE
    ]
    if (nrow(gp) == 0L) next
    eqtl <- load_gene_eqtl(
      gp$gene_ensembl, gp$chr, gp$start, gp$end
    )
    if (nrow(eqtl) < config$coloc$minimum_common_snps) next

    fit <- run_coloc_pair(
      eqtl, gwas_list[[target$trait]], target$gene_symbol, target$trait
    )
    if (nrow(fit) == 0L) next
    fit[, `:=`(
      gene_ensembl = target$gene_ensembl,
      gene_layer = target$gene_layer
    )]
    results[[paste(target$gene_ensembl, target$trait, sep = "__")]] <- fit

    primary <- fit[p12 == config$coloc$primary_p12]
    if (config$optional$make_locus_plots &&
        config$optional$coordinates_harmonised_for_plots &&
        nrow(primary) == 1L &&
        primary$PP.H4 >= config$coloc$tier_a_pp4) {
      make_locus_plot(eqtl, gwas_list[[target$trait]],
                      target$gene_symbol, target$trait)
    }
  }
  bind_nonempty(results)
}

if (config$coloc$run) {
  log_step("Running colocalisation and prior-sensitivity analyses")
  coloc_all <- run_colocalisation(mr_main, gene_positions, gwas_list)

  if (nrow(coloc_all) > 0L) {
    coloc_all[, evidence_tier := dplyr::case_when(
      p12 != config$coloc$primary_p12 ~ "prior_sensitivity",
      PP.H4 >= config$coloc$tier_a_pp4 ~ "Tier_A",
      PP.H4 >= config$coloc$tier_b_pp4 & PP.H4 > PP.H3 ~ "Tier_B",
      PP.H4 >= config$coloc$tier_c_pp4 & PP.H4 > PP.H3 ~ "Tier_C",
      TRUE ~ "None"
    )]
    coloc_primary <- coloc_all[p12 == config$coloc$primary_p12]
  } else {
    coloc_primary <- empty_dt()
  }
} else {
  coloc_all <- empty_dt()
  coloc_primary <- empty_dt()
}

write_result(coloc_all, "04_coloc_all_priors.csv")
write_result(coloc_primary, "04_coloc_primary_prior.csv")
if (nrow(coloc_primary) > 0L) {
  write_result(
    coloc_primary[evidence_tier == "Tier_A"],
    "04_coloc_TierA_prioritised_candidates.csv"
  )
}


## 08b. Auditable multi-signal colocalisation with SuSiE ---------------------

initial_multisignal_qc <- function(target, gp, aligned = NULL) {
  data.table(
    gene_symbol = target$gene_symbol,
    gene_ensembl = target$gene_ensembl,
    gene_layer = target$gene_layer,
    trait = target$trait,
    chromosome = as.character(gp$chr),
    region_start = as.integer(gp$start),
    region_end = as.integer(gp$end),
    coloc_abf_n_snps = as.integer(target$n_snps),
    coloc_abf_primary_PP.H3 = as.numeric(target$PP.H3),
    coloc_abf_primary_PP.H4 = as.numeric(target$PP.H4),
    ld_reference_panel = paste0("1000G_", config$instruments$ancestry),
    ld_reference_prefix = config$instruments$ld_reference,
    plink_binary = config$instruments$plink_binary,
    plink_ld_command = "--r square gz --keep-allele-order",
    n_harmonised_snps = if (is.null(aligned)) NA_integer_ else nrow(aligned$eqtl),
    n_ld_reference_snps = NA_integer_,
    n_ld_allele_matched_snps = NA_integer_,
    n_susie_input_snps = NA_integer_,
    susie_input_limit = as.integer(config$coloc$susie_max_snps),
    variant_selection_status = NA_character_,
    variant_selection_rule = NA_character_,
    susie_status = NA_character_,
    susie_message = NA_character_,
    n_eqtl_credible_sets = NA_integer_,
    n_outcome_credible_sets = NA_integer_,
    n_primary_signal_pairs = NA_integer_,
    max_primary_PP.H3.abf = NA_real_,
    max_primary_PP.H4.abf = NA_real_,
    top_primary_signal_pair = NA_character_,
    locus_plot = NA_character_
  )
}

# [REVISED: Reviewer 5] Memory guard is not a significance-based SNP filter.
select_susie_input_snps <- function(e, g, bim, max_snps) {
  if (nrow(bim) > max_snps) stop(
    "Complete LD-matched region exceeds the SuSiE memory guard (", nrow(bim),
    " > ", max_snps, "). Increase susie_max_snps with sufficient memory; ",
    "no variants have been selected by association significance.")
  list(snps = bim$SNP, status = "complete_ld_region",
       rule = "All LD-reference allele-matched variants were used.")
}

prepare_local_ld_audited <- function(aligned, chromosome, target, gp) {
  qc <- initial_multisignal_qc(target, gp, aligned)
  chromosome <- as.character(chromosome)
  if (!chromosome %in% as.character(1:22)) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "Local LD reference covers autosomes 1-22 only."
    )]
    return(list(qc = qc))
  }

  bfile <- paste0(config$instruments$ld_reference, ".", chromosome)
  if (!file.exists(config$instruments$plink_binary) ||
      !file.exists(paste0(bfile, ".bed")) ||
      !file.exists(paste0(bfile, ".bim")) ||
      !file.exists(paste0(bfile, ".fam"))) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "PLINK binary or chromosome-specific LD reference files are missing."
    )]
    return(list(qc = qc))
  }

  e <- data.table::copy(aligned$eqtl)
  g <- data.table::copy(aligned$gwas)
  bim <- tryCatch(
    data.table::fread(
      paste0(bfile, ".bim"),
      header = FALSE,
      select = c(2L, 4L, 5L, 6L),
      col.names = c("SNP", "bp", "A1", "A2"),
      showProgress = FALSE
    ),
    error = function(err) NULL
  )
  if (is.null(bim) || nrow(bim) == 0L) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "The chromosome-specific BIM file could not be read."
    )]
    return(list(qc = qc))
  }

  bim[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  bim <- unique(bim[SNP %in% intersect(e$SNP, g$SNP)], by = "SNP")
  qc[, n_ld_reference_snps := nrow(bim)]
  if (nrow(bim) < config$coloc$minimum_common_snps) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "Fewer than the required common SNPs were present in the LD reference."
    )]
    return(list(qc = qc))
  }

  e_bim <- e[match(bim$SNP, SNP)]
  valid <- (
    e_bim$effect_allele == bim$A1 & e_bim$other_allele == bim$A2
  ) | (
    e_bim$effect_allele == bim$A2 & e_bim$other_allele == bim$A1
  )
  bim <- bim[valid]
  qc[, n_ld_allele_matched_snps := nrow(bim)]
  if (nrow(bim) < config$coloc$minimum_common_snps) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "Fewer than the required SNPs remained after LD allele matching."
    )]
    return(list(qc = qc))
  }

  selection <- select_susie_input_snps(
    e, g, bim, as.integer(config$coloc$susie_max_snps)
  )
  bim <- bim[match(selection$snps, SNP)]
  qc[, `:=`(
    n_susie_input_snps = nrow(bim),
    variant_selection_status = selection$status,
    variant_selection_rule = selection$rule
  )]

  tmp_prefix <- tempfile(pattern = "susie_ld_")
  on.exit(
    unlink(list.files(
      dirname(tmp_prefix),
      pattern = paste0("^", basename(tmp_prefix)),
      full.names = TRUE
    ), force = TRUE),
    add = TRUE
  )
  snp_file <- paste0(tmp_prefix, ".extract")
  data.table::fwrite(bim[, .(SNP)], snp_file, col.names = FALSE)
  output <- tryCatch(
    system2(
      config$instruments$plink_binary,
      args = c(
        "--bfile", shQuote(bfile),
        "--extract", shQuote(snp_file),
        "--keep-allele-order",
        "--r", "square", "gz",
        "--write-snplist",
        "--out", shQuote(tmp_prefix)
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(err) conditionMessage(err)
  )
  if (!file.exists(paste0(tmp_prefix, ".ld.gz")) ||
      !file.exists(paste0(tmp_prefix, ".snplist"))) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = paste(
        "PLINK LD calculation failed:",
        paste(tail(output, 5L), collapse = " ")
      )
    )]
    return(list(qc = qc))
  }

  snp_order <- readLines(paste0(tmp_prefix, ".snplist"), warn = FALSE)
  ld <- as.matrix(data.table::fread(
    paste0(tmp_prefix, ".ld.gz"), header = FALSE, showProgress = FALSE
  ))
  if (length(snp_order) != nrow(ld) || nrow(ld) != ncol(ld)) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "The PLINK LD matrix dimensions did not match the SNP list."
    )]
    return(list(qc = qc))
  }

  dimnames(ld) <- list(snp_order, snp_order)
  ld[!is.finite(ld)] <- 0
  ld <- (ld + t(ld)) / 2
  diag(ld) <- 1

  bim <- bim[match(snp_order, SNP)]
  e <- e[match(snp_order, SNP)]
  g <- g[match(snp_order, SNP)]
  flip <- e$effect_allele == bim$A2 & e$other_allele == bim$A1
  if (any(flip)) {
    old_e_effect <- e$effect_allele[flip]
    old_g_effect <- g$effect_allele[flip]
    e[flip, `:=`(
      beta = -beta,
      eaf = 1 - eaf,
      effect_allele = other_allele,
      other_allele = old_e_effect
    )]
    g[flip, `:=`(
      beta = -beta,
      eaf = 1 - eaf,
      effect_allele = other_allele,
      other_allele = old_g_effect
    )]
  }

  qc[, susie_status := "ld_prepared"]
  list(qc = qc, eqtl = e, gwas = g, LD = ld)
}

make_multisignal_locus_plot <- function(aligned, ld_data, target) {
  selected_snps <- ld_data$eqtl$SNP
  e <- data.table::copy(aligned$eqtl)
  g <- data.table::copy(aligned$gwas)
  plot_dt <- rbind(
    data.table(
      # [REVISED: Reviewer 6] Common plotting build; no effect estimates changed.
      SNP = e$SNP, pos = g$pos[match(e$SNP, g$SNP)], pval = e$pval,
      dataset = "Skeletal-muscle eQTL",
      retained_for_susie = e$SNP %in% selected_snps
    ),
    data.table(
      SNP = g$SNP, pos = g$pos, pval = g$pval,
      dataset = target$trait,
      retained_for_susie = g$SNP %in% selected_snps
    )
  )
  plot_dt[, minus_log10_p := -log10(pmax(pval, 1e-300))]
  lead_labels <- plot_dt[
    is.finite(minus_log10_p),
    .SD[order(-minus_log10_p)][seq_len(min(.N, 2L))],
    by = dataset
  ]
  plot_dt[, label := fifelse(SNP %in% lead_labels$SNP, SNP, NA_character_)]

  p <- ggplot(
    plot_dt,
    aes(pos / 1e6, minus_log10_p, colour = retained_for_susie)
  ) +
    geom_point(alpha = 0.65, size = 1.1) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 2.8,
      min.segment.length = 0,
      max.overlaps = Inf,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    facet_wrap(~dataset, ncol = 1, scales = "free_y") +
    scale_colour_manual(
      values = c("FALSE" = "grey78", "TRUE" = "#3366AA"),
      labels = c("FALSE" = "Not used in SuSiE", "TRUE" = "Used in SuSiE")
    ) +
    labs(
      title = paste(target$gene_symbol, target$trait, sep = " - "),
      subtitle = paste0(
        "Harmonised regional variants with the LD-feasible SuSiE set highlighted"
      ),
      x = "GRCh37 genomic position (Mb)",
      y = expression(-log[10](P)),
      colour = NULL
    ) +
    theme_manuscript()

  out_file <- plot_path(paste0(
    "coloc_multisignal_locus_", target$gene_symbol, "_", target$trait, ".png"
  ))
  png_device <- function(filename, width, height, ...) {
    grDevices::png(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = 300
    )
  }
  tryCatch(
    ggsave(
      out_file, p, width = 8, height = 6, dpi = 300,
      device = png_device
    ),
    warning = function(warn) warning(conditionMessage(warn)),
    error = function(err) warning("Locus plot was not written: ",
                                  conditionMessage(err))
  )
  pdf_file <- sub("\\.png$", ".pdf", out_file)
  tryCatch(
    ggsave(
      pdf_file, p, width = 8, height = 6,
      device = grDevices::cairo_pdf
    ),
    warning = function(warn) warning(conditionMessage(warn)),
    error = function(err) warning("Locus PDF was not written: ",
                                  conditionMessage(err))
  )
  out_file
}

make_susie_dataset_pair <- function(ld_data, trait) {
  e <- ld_data$eqtl
  g <- ld_data$gwas
  d_eqtl <- list(
    snp = e$SNP,
    beta = e$beta,
    varbeta = e$se^2,
    MAF = e$maf,
    N = config$gtex$sample_size,
    type = "quant",
    sdY = 1,
    LD = ld_data$LD
  )
  d_gwas <- list(
    snp = g$SNP,
    beta = g$beta,
    varbeta = g$se^2,
    MAF = g$maf,
    N = stats::median(g$N, na.rm = TRUE),
    LD = ld_data$LD
  )
  if (config$outcome_type[[trait]] == "cc") {
    d_gwas$type <- "cc"
    d_gwas$s <- config$outcome_cases[[trait]] / config$outcome_n[[trait]]
  } else {
    sdY_gwas <- sqrt(stats::median(
      2 * g$maf * (1 - g$maf) * g$N * g$se^2,
      na.rm = TRUE
    ))
    if (!is.finite(sdY_gwas) || sdY_gwas <= 0) return(NULL)
    d_gwas$type <- "quant"
    d_gwas$sdY <- sdY_gwas
  }
  list(eqtl = d_eqtl, gwas = d_gwas)
}

get_susie_set_count <- function(sfit) {
  if (is.null(sfit) || is.null(sfit$sets) || is.null(sfit$sets$cs)) return(0L)
  length(sfit$sets$cs)
}

get_susie_purity <- function(sfit, set_number, pattern) {
  if (is.null(sfit$sets$purity)) return(NA_real_)
  purity <- as.data.frame(sfit$sets$purity)
  if (nrow(purity) < set_number) return(NA_real_)
  col <- grep(pattern, names(purity), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(col)) return(NA_real_)
  as.numeric(purity[[col]][set_number])
}

extract_susie_credible_sets <- function(sfit, data_dt, dataset, target) {
  if (get_susie_set_count(sfit) == 0L) return(empty_dt())
  cs_list <- sfit$sets$cs
  cs_index <- sfit$sets$cs_index
  if (is.null(cs_index)) cs_index <- seq_along(cs_list)

  rows <- lapply(seq_along(cs_list), function(i) {
    members <- names(cs_list[[i]])
    if (is.null(members)) members <- as.character(cs_list[[i]])
    pips <- sfit$pip[members]
    pip_order <- order(
      ifelse(is.finite(pips), -pips, Inf),
      members
    )
    ordered_members <- members[pip_order]
    component <- as.integer(cs_index[i])
    alpha_cols <- intersect(members, colnames(sfit$alpha))
    coverage <- if (length(alpha_cols) > 0L &&
                    is.finite(component) &&
                    component <= nrow(sfit$alpha)) {
      sum(sfit$alpha[component, alpha_cols], na.rm = TRUE)
    } else {
      NA_real_
    }
    positions <- data_dt$pos[match(members, data_dt$SNP)]
    finite_positions <- positions[is.finite(positions)]
    lead_snp <- ordered_members[1]
    data.table(
      gene_symbol = target$gene_symbol,
      gene_ensembl = target$gene_ensembl,
      trait = target$trait,
      dataset = dataset,
      # [NEW: Reviewer 6] These are source coordinates, not a shared build.
      position_build = if (dataset == "skeletal_muscle_eQTL") "GRCh38" else "GRCh37",
      credible_set = paste0(dataset, "_CS", i),
      susie_component = component,
      cs_size = length(members),
      lead_snp = lead_snp,
      lead_pip = as.numeric(sfit$pip[lead_snp]),
      coverage = coverage,
      purity_min_abs_corr = get_susie_purity(sfit, i, "min"),
      purity_mean_abs_corr = get_susie_purity(sfit, i, "mean"),
      min_pos = if (length(finite_positions) > 0L) {
        min(finite_positions)
      } else {
        NA_integer_
      },
      max_pos = if (length(finite_positions) > 0L) {
        max(finite_positions)
      } else {
        NA_integer_
      },
      member_snps = paste(head(ordered_members, 100L), collapse = ";"),
      member_snps_truncated = length(ordered_members) > 100L
    )
  })
  bind_nonempty(rows)
}

extract_susie_credible_set_members <- function(sfit, data_dt, dataset, target) {
  if (get_susie_set_count(sfit) == 0L) return(empty_dt())
  cs_list <- sfit$sets$cs
  cs_index <- sfit$sets$cs_index
  if (is.null(cs_index)) cs_index <- seq_along(cs_list)

  rows <- lapply(seq_along(cs_list), function(i) {
    members <- names(cs_list[[i]])
    if (is.null(members)) members <- as.character(cs_list[[i]])
    source <- data_dt[match(members, SNP)]
    component <- as.integer(cs_index[i])
    data.table(
      gene_symbol = target$gene_symbol,
      gene_ensembl = target$gene_ensembl,
      trait = target$trait,
      dataset = dataset,
      # [NEW: Reviewer 6] These are source coordinates, not a shared build.
      position_build = if (dataset == "skeletal_muscle_eQTL") "GRCh38" else "GRCh37",
      credible_set = paste0(dataset, "_CS", i),
      susie_component = component,
      SNP = members,
      pos = source$pos,
      beta = source$beta,
      se = source$se,
      pval = source$pval,
      pip = as.numeric(sfit$pip[members])
    )[order(-pip, pval, SNP)]
  })
  bind_nonempty(rows)
}

run_coloc_susie_priors <- function(eqtl_fit, gwas_fit, target, n_ld_snps) {
  rows <- lapply(config$coloc$p12_values, function(p12) {
    fit <- tryCatch(
      coloc::coloc.susie(
        dataset1 = eqtl_fit,
        dataset2 = gwas_fit,
        p1 = config$coloc$p1,
        p2 = config$coloc$p2,
        p12 = p12
      ),
      error = function(err) NULL
    )
    if (is.null(fit) || is.null(fit$summary) || nrow(fit$summary) == 0L) {
      return(NULL)
    }
    out <- data.table::as.data.table(fit$summary, keep.rownames = "signal_pair")
    if (!"PP.H3.abf" %in% names(out)) out[, PP.H3.abf := NA_real_]
    if (!"PP.H4.abf" %in% names(out)) out[, PP.H4.abf := NA_real_]
    if (!"hit1" %in% names(out)) out[, hit1 := NA_character_]
    if (!"hit2" %in% names(out)) out[, hit2 := NA_character_]
    if (all(c("hit1", "hit2") %in% names(out))) {
      out[, signal_pair := paste(hit1, hit2, sep = "__")]
    }
    out[, `:=`(
      gene_symbol = target$gene_symbol,
      gene_ensembl = target$gene_ensembl,
      gene_layer = target$gene_layer,
      trait = target$trait,
      p1 = config$coloc$p1,
      p2 = config$coloc$p2,
      p12 = p12,
      n_ld_snps = n_ld_snps,
      PP.H3 = PP.H3.abf,
      PP.H4 = PP.H4.abf
    )]
    out
  })
  bind_nonempty(rows)
}

run_multisignal_coloc_pair <- function(target, gp, gwas_list) {
  eqtl <- load_gene_eqtl(gp$gene_ensembl, gp$chr, gp$start, gp$end)
  aligned <- harmonise_coloc_region(eqtl, gwas_list[[target$trait]])
  if (is.null(aligned)) {
    qc <- initial_multisignal_qc(target, gp)
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "Fewer than the required harmonised SNPs were available."
    )]
    return(list(qc = qc))
  }

  ld_data <- prepare_local_ld_audited(aligned, gp$chr, target, gp)
  qc <- ld_data$qc
  if (is.null(ld_data$LD)) return(list(qc = qc))

  plot_file <- make_multisignal_locus_plot(aligned, ld_data, target)
  qc[, locus_plot := plot_file]

  datasets <- make_susie_dataset_pair(ld_data, target$trait)
  if (is.null(datasets)) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = "Outcome sdY could not be estimated for SuSiE."
    )]
    return(list(qc = qc))
  }

  eqtl_fit <- tryCatch(
    coloc::runsusie(datasets$eqtl, suffix = 1, maxit = 200),
    error = function(err) err
  )
  gwas_fit <- tryCatch(
    coloc::runsusie(datasets$gwas, suffix = 2, maxit = 200),
    error = function(err) err
  )
  if (inherits(eqtl_fit, "error") || inherits(gwas_fit, "error")) {
    qc[, `:=`(
      susie_status = "not_run",
      susie_message = paste(
        "SuSiE failed:",
        if (inherits(eqtl_fit, "error")) conditionMessage(eqtl_fit) else "eQTL ok",
        if (inherits(gwas_fit, "error")) conditionMessage(gwas_fit) else "outcome ok"
      )
    )]
    return(list(qc = qc))
  }

  summary_all <- run_coloc_susie_priors(
    eqtl_fit, gwas_fit, target, nrow(ld_data$LD)
  )
  credible_sets <- bind_nonempty(list(
    extract_susie_credible_sets(
      eqtl_fit, ld_data$eqtl, "skeletal_muscle_eQTL", target
    ),
    extract_susie_credible_sets(
      gwas_fit, ld_data$gwas, paste0("outcome_", target$trait), target
    )
  ))
  credible_set_members <- bind_nonempty(list(
    extract_susie_credible_set_members(
      eqtl_fit, ld_data$eqtl, "skeletal_muscle_eQTL", target
    ),
    extract_susie_credible_set_members(
      gwas_fit, ld_data$gwas, paste0("outcome_", target$trait), target
    )
  ))

  primary <- if (nrow(summary_all) > 0L && "p12" %in% names(summary_all)) {
    summary_all[p12 == config$coloc$primary_p12]
  } else {
    empty_dt()
  }
  if (nrow(primary) > 0L) {
    top <- primary[order(-PP.H4.abf, -PP.H3.abf)][1L]
    qc[, `:=`(
      susie_status = "completed",
      susie_message = "completed",
      n_eqtl_credible_sets = get_susie_set_count(eqtl_fit),
      n_outcome_credible_sets = get_susie_set_count(gwas_fit),
      n_primary_signal_pairs = nrow(primary),
      max_primary_PP.H3.abf = max(primary$PP.H3.abf, na.rm = TRUE),
      max_primary_PP.H4.abf = max(primary$PP.H4.abf, na.rm = TRUE),
      top_primary_signal_pair = top$signal_pair
    )]
  } else {
    qc[, `:=`(
      susie_status = "completed_no_colocalising_cs_pair",
      susie_message = "SuSiE produced no tested credible-set pair."
    )]
  }

  list(
    qc = qc,
    summary_all = summary_all,
    credible_sets = credible_sets,
    credible_set_members = credible_set_members
  )
}

summarise_tier_a_priors <- function(coloc_all, susie_all) {
  tier_targets <- unique(
    coloc_primary[evidence_tier == "Tier_A",
                  .(gene_symbol, gene_ensembl, trait)]
  )
  abf <- merge(
    tier_targets,
    coloc_all[, .(
      gene_symbol, gene_ensembl, trait, p1, p2, p12,
      PP.H3, PP.H4, n_snps
    )],
    by = c("gene_symbol", "gene_ensembl", "trait"),
    all.x = FALSE
  )
  if (nrow(abf) > 0L) {
    abf[, `:=`(
      method = "coloc.abf_single_signal",
      n_signal_pairs = NA_integer_,
      top_signal_pair = NA_character_,
      n_variants_used = n_snps
    )]
  }

  susie <- empty_dt()
  if (nrow(susie_all) > 0L) {
    susie_counts <- susie_all[
      ,
      .(n_signal_pairs = .N),
      by = .(gene_symbol, gene_ensembl, trait, p1, p2, p12)
    ]
    susie_top <- susie_all[
      ,
      .SD[order(-PP.H4.abf, -PP.H3.abf)][1L],
      by = .(gene_symbol, gene_ensembl, trait, p1, p2, p12)
    ]
    susie <- merge(
      susie_top,
      susie_counts,
      by = c("gene_symbol", "gene_ensembl", "trait", "p1", "p2", "p12"),
      all.x = TRUE,
      sort = FALSE
    )[
      ,
      .(
        gene_symbol, gene_ensembl, trait, p1, p2, p12,
        PP.H3 = PP.H3.abf,
        PP.H4 = PP.H4.abf,
        n_signal_pairs,
        top_signal_pair = signal_pair,
        n_variants_used = n_ld_snps,
        method = "coloc.susie_top_signal_pair"
      )
    ]
  }

  out <- bind_nonempty(list(
    abf[, .(
      gene_symbol, gene_ensembl, trait, method, p1, p2, p12,
      n_variants_used, n_signal_pairs, top_signal_pair, PP.H3, PP.H4
    )],
    susie[, .(
      gene_symbol, gene_ensembl, trait, method, p1, p2, p12,
      n_variants_used, n_signal_pairs, top_signal_pair, PP.H3, PP.H4
    )]
  ))
  if (nrow(out) > 0L) {
    out[order(gene_symbol, trait, method, p12)]
  } else {
    out
  }
}

susie_coloc <- empty_dt()
susie_coloc_all_priors <- empty_dt()
susie_multisignal_qc <- empty_dt()
susie_credible_sets <- empty_dt()
susie_credible_set_members <- empty_dt()
tier_a_prior_summary <- empty_dt()
abcc8_pp4_stability <- empty_dt()

if (config$coloc$run_susie && nrow(coloc_primary) > 0L) {
  susie_targets <- coloc_primary[evidence_tier == "Tier_A"]
  susie_runs <- vector("list", nrow(susie_targets))
  for (i in seq_len(nrow(susie_targets))) {
    target <- susie_targets[i]
    gp <- gene_positions[
      gene_positions$gene_ensembl == target$gene_ensembl,
      ,
      drop = FALSE
    ]
    if (nrow(gp) == 0L) next
    gp <- gp[1L, , drop = FALSE]
    log_step(
      "Running audited coloc-SuSiE: ",
      target$gene_symbol, " - ", target$trait
    )
    susie_runs[[i]] <- run_multisignal_coloc_pair(target, gp, gwas_list)
  }
  susie_runs <- Filter(Negate(is.null), susie_runs)
  susie_multisignal_qc <- bind_nonempty(lapply(susie_runs, `[[`, "qc"))
  susie_coloc_all_priors <- bind_nonempty(lapply(susie_runs, `[[`, "summary_all"))
  susie_coloc <- if (
    nrow(susie_coloc_all_priors) > 0L &&
      "p12" %in% names(susie_coloc_all_priors)
  ) {
    susie_coloc_all_priors[p12 == config$coloc$primary_p12]
  } else {
    empty_dt()
  }
  susie_credible_sets <- bind_nonempty(lapply(susie_runs, `[[`, "credible_sets"))
  susie_credible_set_members <- bind_nonempty(
    lapply(susie_runs, `[[`, "credible_set_members")
  )
  tier_a_prior_summary <- summarise_tier_a_priors(
    coloc_all, susie_coloc_all_priors
  )
  if (nrow(tier_a_prior_summary) > 0L) {
    abcc8_pp4_stability <- tier_a_prior_summary[gene_symbol == "ABCC8"]
  }
}

write_result(susie_coloc, "04b_coloc_SuSiE_TierA.csv")
write_result(susie_multisignal_qc, "04c_coloc_SuSiE_TierA_QC.csv")
write_result(susie_coloc_all_priors, "04d_coloc_SuSiE_all_priors.csv")
write_result(susie_credible_sets, "04e_coloc_SuSiE_credible_sets.csv")
write_result(susie_credible_set_members, "04f_coloc_SuSiE_CS_members.csv")
write_result(tier_a_prior_summary, "04g_TierA_coloc_prior_summary.csv")
write_result(abcc8_pp4_stability, "04h_ABCC8_PP4_stability.csv")

## Reviewer 1 revision: explicit multi-signal audit -------------------------
## This append-only audit does not change the coloc classification. It records
## whether each Tier A locus used the complete LD-matched region and whether
## an interpretable paired credible-set result was available. PP3/PP4 values
## are intentionally summarized only for loci with a tested signal pair.
if (nrow(susie_multisignal_qc) > 0L) {
  susie_complete_region_audit <- susie_multisignal_qc[, .(
    gene_symbol,
    trait,
    n_coloc_abf_snps = coloc_abf_n_snps,
    n_ld_reference_snps,
    n_ld_allele_matched_snps,
    n_susie_input_snps,
    susie_input_limit,
    variant_selection_status,
    susie_status,
    n_eqtl_credible_sets,
    n_outcome_credible_sets,
    n_primary_signal_pairs,
    interpretable_paired_susie = susie_status == "completed" &
      n_primary_signal_pairs > 0L,
    locus_plot
  )]
} else {
  susie_complete_region_audit <- empty_dt()
}

if (nrow(susie_coloc_all_priors) > 0L) {
  susie_prior_applicability <- susie_coloc_all_priors[, .(
    n_signal_pairs = uniqueN(signal_pair),
    n_rows = .N,
    p12_values = paste(sort(unique(p12)), collapse = ";"),
    min_PP3 = min(PP.H3, na.rm = TRUE),
    max_PP3 = max(PP.H3, na.rm = TRUE),
    min_PP4 = min(PP.H4, na.rm = TRUE),
    max_PP4 = max(PP.H4, na.rm = TRUE),
    prior_results_available = TRUE
  ), by = .(gene_symbol, gene_ensembl, trait)]
} else {
  susie_prior_applicability <- empty_dt()
}

write_result(susie_complete_region_audit,
             "04i_coloc_SuSiE_complete_region_audit.csv")
write_result(susie_prior_applicability,
             "04j_coloc_SuSiE_prior_applicability.csv")


## 09. Canonical versus expanded-gene sensitivity ---------------------------

layer_summary <- mr_main[, .(
  n_tests = .N,
  n_genes = uniqueN(gene_symbol),
  n_mr_associated = sum(mr_associated, na.rm = TRUE),
  proportion_mr_associated = mean(mr_associated, na.rm = TRUE),
  n_layer_FDR_significant = sum(mr_associated_layered, na.rm = TRUE),
  proportion_layer_FDR_significant = mean(mr_associated_layered, na.rm = TRUE),
  n_positive = sum(beta > 0, na.rm = TRUE),
  n_negative = sum(beta < 0, na.rm = TRUE)
), by = .(trait, gene_layer)]

layer_tests <- lapply(unique(mr_main$trait), function(trait_name) {
  x <- mr_main[trait == trait_name]
  tab <- table(x$gene_layer, x$mr_associated)
  fisher_p <- if (all(dim(tab) == c(2L, 2L))) fisher.test(tab)$p.value else NA_real_
  wilcox_p <- tryCatch(
    wilcox.test(beta ~ gene_layer, data = x, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  data.table(trait = trait_name, fisher_p = fisher_p, wilcoxon_p = wilcox_p)
}) |>
  bind_nonempty()

if (nrow(layer_tests) > 0L) {
  layer_tests[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
  layer_tests[, wilcoxon_fdr := p.adjust(wilcoxon_p, method = "BH")]
}

if (nrow(coloc_primary) > 0L) {
  coloc_layer_summary <- coloc_primary[, .(
    n_coloc_tested = .N,
    n_tier_A = sum(evidence_tier == "Tier_A", na.rm = TRUE),
    n_tier_B = sum(evidence_tier == "Tier_B", na.rm = TRUE)
  ), by = .(trait, gene_layer)]
} else {
  coloc_layer_summary <- empty_dt()
}

write_result(layer_summary, "05_canonical_expanded_summary.csv")
write_result(layer_tests, "05_canonical_expanded_tests.csv")
write_result(coloc_layer_summary, "05_canonical_expanded_coloc_summary.csv")
write_result(
  mr_main[gene_layer == "canonical"],
  "05_canonical_only_cisMR_results.csv"
)


## 10. Exploratory GO enrichment ---------------------------------------------

if (config$optional$run_go_enrichment) {
  associated_genes <- unique(mr_main[mr_associated == TRUE, gene_symbol])
  tested_genes <- unique(mr_main$gene_symbol)

  if (length(associated_genes) >= 10L) {
    go_fit <- tryCatch(
      clusterProfiler::enrichGO(
        gene = associated_genes,
        universe = tested_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500
      ),
      error = function(e) NULL
    )

    if (!is.null(go_fit)) {
      go_results <- as.data.frame(go_fit)
      write_result(go_results, "06_GO_BP_exploratory.csv")
      if (nrow(go_results) > 0L) {
        p_go <- enrichplot::dotplot(go_fit, showCategory = 20) +
          labs(
            title = "Exploratory GO:BP enrichment",
            subtitle = "Universe: IR genes entering cis-MR"
          ) +
          theme_manuscript()
        ggsave(plot_path("GO_BP_exploratory.png"), p_go,
               width = 9, height = 8, dpi = 300)
      }
    }
  }
}


## 11. Main tables and figures ------------------------------------------------

mr_overview <- mr_main[, .(
  n_tested_genes = uniqueN(gene_symbol),
  n_MR_associated_pairs = sum(mr_associated, na.rm = TRUE),
  n_MR_associated_genes = uniqueN(gene_symbol[mr_associated == TRUE]),
  median_F = median(mean_F, na.rm = TRUE),
  minimum_F = min(min_F, na.rm = TRUE),
  n_single_SNP = sum(nsnp == 1L),
  n_two_SNP = sum(nsnp == 2L),
  n_three_or_more_SNP = sum(nsnp >= 3L)
), by = trait]

main_gene_table <- mr_main[mr_associated == TRUE, .(
  gene_symbol, gene_ensembl, gene_layer, trait, method, nsnp,
  beta, se, ci_lower, ci_upper, pval, fdr_within_trait, fdr_global,
  fdr_within_trait_layer, min_F, mean_F
)]

if (nrow(coloc_primary) > 0L) {
  evidence_table <- merge(
    main_gene_table,
    coloc_primary[, .(
      gene_symbol, trait, n_coloc_snps = n_snps,
      PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, evidence_tier
    )],
    by = c("gene_symbol", "trait"),
    all.x = TRUE
  )
} else {
  evidence_table <- copy(main_gene_table)
  evidence_table[, `:=`(
    n_coloc_snps = NA_integer_, PP.H0 = NA_real_, PP.H1 = NA_real_,
    PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
    evidence_tier = "coloc_not_available"
  )]
}

write_result(mr_overview, "07_Manuscript_Table1_cisMR_overview.csv")
write_result(evidence_table, "07_Manuscript_Table2_evidence_hierarchy.csv")

forest_data <- evidence_table |>
  dplyr::mutate(
    label = ifelse(evidence_tier == "Tier_A", paste0(gene_symbol, " *"), gene_symbol),
    trait = factor(trait, levels = names(config$files$outcomes))
  )

if (nrow(forest_data) > 0L) {
  p_forest <- ggplot(
    forest_data,
    aes(x = beta, y = reorder(label, beta), xmin = ci_lower, xmax = ci_upper,
        colour = trait)
  ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(height = 0.15) +
    geom_point(size = 2) +
    facet_wrap(~trait, scales = "free_y") +
    labs(
      title = "FDR-significant skeletal-muscle cis-MR associations",
      subtitle = "* Tier A colocalisation support",
      x = "MR effect estimate (95% CI)", y = NULL, colour = "Trait"
    ) +
    theme_manuscript()
  ggsave(plot_path("Figure_cisMR_forest.png"), p_forest,
         width = 11, height = 9, dpi = 300)
  ggsave(plot_path("Figure_cisMR_forest.pdf"), p_forest,
         width = 11, height = 9)
}

if (nrow(coloc_primary) > 0L) {
  p_coloc <- ggplot(
    coloc_primary,
    aes(x = PP.H3, y = PP.H4, colour = trait, label = gene_symbol)
  ) +
    geom_hline(yintercept = config$coloc$tier_a_pp4, linetype = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
    geom_point(size = 2.3) +
    ggrepel::geom_text_repel(
      data = coloc_primary[evidence_tier %in% c("Tier_A", "Tier_B")],
      size = 3, max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Colocalisation evidence",
      subtitle = "PP.H4 versus PP.H3 under the primary prior",
      x = "PP.H3: distinct causal variants",
      y = "PP.H4: shared causal variant",
      colour = "Trait"
    ) +
    theme_manuscript()
  ggsave(plot_path("Figure_coloc_PP3_PP4.png"), p_coloc,
         width = 8, height = 6, dpi = 300)
  ggsave(plot_path("Figure_coloc_PP3_PP4.pdf"), p_coloc,
         width = 8, height = 6)
}


## 12. Independent FUSION skeletal-muscle replication ------------------------

read_fusion_targets <- function(path, target_gene_ids, chunk_lines = 250000L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  selected <- list()
  chunk_number <- 0L

  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (length(lines) == 0L) break
    chunk_number <- chunk_number + 1L
    phenotype_id <- data.table::tstrsplit(
      lines, "\t", fixed = TRUE, keep = 7L
    )[[1]]
    keep <- strip_ensembl_version(phenotype_id) %in% target_gene_ids
    if (any(keep)) selected[[length(selected) + 1L]] <- lines[keep]
    if (chunk_number %% 20L == 0L) {
      log_step("FUSION rows scanned: ", format(chunk_number * chunk_lines, big.mark = ","))
    }
  }

  if (length(selected) == 0L) return(empty_dt())
  dt <- data.table::fread(
    text = paste(unlist(selected, use.names = FALSE), collapse = "\n"),
    header = FALSE,
    sep = "\t",
    showProgress = FALSE
  )
  data.table::setnames(dt, header)
  dt
}

standardise_fusion_eqtl <- function(dt, gene_id) {
  dt <- data.table::as.data.table(dt)[
    strip_ensembl_version(pheno_id) == strip_ensembl_version(gene_id)
  ]
  if (nrow(dt) == 0L) return(empty_dt())
  dt <- dt[, .(
    SNP = as.character(var_id),
    beta = as.numeric(beta),
    se = as.numeric(ste),
    pval = as.numeric(pv),
    eaf = as.numeric(effect_allele_freq),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele)),
    pos = as.integer(var_pos),
    chr = as.character(var_chr),
    gene_ensembl = strip_ensembl_version(gene_id)
  )]
  standardise_eqtl(dt, gene_id)
}

run_fusion_replication <- function(targets, fusion_dt, gwas_list) {
  rows <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    eqtl <- standardise_fusion_eqtl(fusion_dt, target$gene_ensembl)
    exposure <- format_exposure(
      eqtl, target$gene_symbol, exposure_n = config$fusion$sample_size
    )
    n_strong <- if (is.null(exposure)) 0L else nrow(exposure)
    exposure <- clump_exposure(exposure)
    n_clumped <- if (is.null(exposure)) 0L else nrow(exposure)

    mr_beta <- mr_se <- mr_pval <- NA_real_
    mr_nsnp <- 0L
    if (!is.null(exposure) && nrow(exposure) > 0L) {
      outcome <- format_outcome(
        gwas_list[[target$trait]], exposure$SNP, target$trait
      )
      if (!is.null(outcome)) {
        dat <- tryCatch(
          TwoSampleMR::harmonise_data(exposure, outcome, action = 3),
          error = function(e) NULL
        )
        if (!is.null(dat)) {
          dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
          if (nrow(dat) > 0L) {
            method <- if (nrow(dat) == 1L) "mr_wald_ratio" else "mr_ivw"
            estimate <- tryCatch(
              TwoSampleMR::mr(dat, method_list = method),
              error = function(e) NULL
            )
            if (!is.null(estimate) && nrow(estimate) > 0L) {
              mr_beta <- estimate$b[1]
              mr_se <- estimate$se[1]
              mr_pval <- estimate$pval[1]
              mr_nsnp <- estimate$nsnp[1]
            }
          }
        }
      }
    }

    coloc_fit <- run_coloc_pair(
      eqtl, gwas_list[[target$trait]], target$gene_symbol, target$trait,
      exposure_n = config$fusion$sample_size
    )
    coloc_fit <- coloc_fit[p12 == config$coloc$primary_p12]
    rows[[i]] <- data.table(
      gene_ensembl = target$gene_ensembl,
      gene_symbol = target$gene_symbol,
      trait = target$trait,
      discovery_beta = target$discovery_beta,
      discovery_PP.H4 = target$discovery_PP.H4,
      n_fusion_region_snps = nrow(eqtl),
      n_fusion_strong_instruments = n_strong,
      n_fusion_clumped_instruments = n_clumped,
      fusion_mr_nsnp = mr_nsnp,
      fusion_mr_beta = mr_beta,
      fusion_mr_se = mr_se,
      fusion_mr_pval = mr_pval,
      direction_concordant = ifelse(
        is.finite(mr_beta), sign(mr_beta) == sign(target$discovery_beta), NA
      ),
      fusion_PP.H3 = if (nrow(coloc_fit) == 1L) coloc_fit$PP.H3 else NA_real_,
      fusion_PP.H4 = if (nrow(coloc_fit) == 1L) coloc_fit$PP.H4 else NA_real_
    )
  }
  result <- bind_nonempty(rows)
  if (nrow(result) > 0L) {
    result[, fusion_mr_fdr := p.adjust(fusion_mr_pval, method = "BH")]
    result[, replicated :=
      direction_concordant %in% TRUE & fusion_mr_fdr < 0.05 & fusion_PP.H4 >= 0.80
    ]
  }
  result
}

fusion_replication <- empty_dt()
if (config$optional$run_fusion_replication && nrow(coloc_primary) > 0L) {
  fusion_targets <- merge(
    coloc_primary[evidence_tier == "Tier_A", .(
      gene_ensembl, gene_symbol, trait, discovery_PP.H4 = PP.H4
    )],
    mr_main[, .(gene_ensembl, trait, discovery_beta = beta)],
    by = c("gene_ensembl", "trait")
  )
  if (nrow(fusion_targets) > 0L) {
    assert_files_exist(config$files$fusion_eqtl, "FUSION eQTL")
    log_step("Reading FUSION records for Tier A discovery genes")
    fusion_dt <- read_fusion_targets(
      config$files$fusion_eqtl, unique(fusion_targets$gene_ensembl)
    )
    fusion_replication <- run_fusion_replication(
      fusion_targets, fusion_dt, gwas_list
    )
  }
}
write_result(fusion_replication, "08_FUSION_independent_replication.csv")


## 13. Optional eQTLGen cross-tissue comparison ------------------------------

if (config$optional$run_eqtlgen_comparison) {
  if (!file.exists(config$files$eqtlgen)) {
    warning("eQTLGen comparison requested but file not found: ",
            config$files$eqtlgen)
  } else {
    log_step("Reading eQTLGen for exploratory cross-tissue comparison")
    eqtlgen <- data.table::fread(config$files$eqtlgen)
    required <- c(
      "SNP", "GeneSymbol", "Pvalue", "AssessedAllele", "OtherAllele",
      "Zscore", "NrSamples"
    )
    if (!all(required %in% names(eqtlgen))) {
      warning("eQTLGen file does not contain all required columns; comparison skipped.")
    } else {
      eqtlgen <- eqtlgen[
        GeneSymbol %in% unique(mr_main[mr_associated == TRUE, gene_symbol])
      ]
      eqtlgen[, `:=`(
        beta_blood = Zscore / sqrt(NrSamples),
        se_blood = 1 / sqrt(NrSamples)
      )]
      data.table::setorder(eqtlgen, GeneSymbol, Pvalue)
      eqtlgen_summary <- eqtlgen[, .(
        n_significant_blood_eqtl = .N,
        top_SNP = first(SNP),
        top_pvalue = first(Pvalue),
        top_beta_blood = first(beta_blood),
        top_se_blood = first(se_blood),
        top_sample_size = first(NrSamples)
      ), by = GeneSymbol]
      write_result(
        eqtlgen_summary,
        "08_eQTLGen_cross_tissue_candidates.csv"
      )
    }
  }
}


## 13b. [NEW: Reviewer 7] GEFOS 2017 ALM external support -----------------------
##
## This is a targeted external outcome-support analysis, not a new discovery
## screen or definitive replication.
## The four ALM Tier A genes are fixed in config$external_alm$target_genes.
## The original GTEx v8 skeletal-muscle eQTL eligibility and LD-clumping rules
## are reused. No proxy SNP substitution is performed. External-outcome
## colocalisation is performed separately below as single-signal coloc.abf
## sensitivity evidence, with all three P12 settings reported.

normalise_external_colname <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

resolve_external_column <- function(
  dt, aliases, label, required = TRUE
) {
  source_names <- names(dt)
  source_keys <- normalise_external_colname(source_names)
  alias_keys <- normalise_external_colname(aliases)

  for (alias_key in alias_keys) {
    hit <- which(source_keys == alias_key)
    if (length(hit) > 0L) return(source_names[hit[1L]])
  }

  if (required) {
    stop(
      "GEFOS ALM file is missing the ", label, " column. Accepted names: ",
      paste(aliases, collapse = ", "),
      ". Observed columns: ", paste(source_names, collapse = ", ")
    )
  }
  NA_character_
}

read_gefos_alm <- function(
  path,
  fixed_n = config$external_alm$sample_size
) {
  log_step("Reading independent ALM outcome: ", path)
  raw <- data.table::fread(
    path,
    na.strings = c("", ".", "NA", "NaN", "-9"),
    check.names = FALSE,
    showProgress = TRUE
  )
  if (nrow(raw) == 0L) stop("The GEFOS ALM file contains no data rows.")

  column_map <- data.table::data.table(
    standard_column = c(
      "SNP", "effect_allele", "other_allele", "eaf",
      "beta", "se", "pval", "N"
    ),
    source_column = c(
      resolve_external_column(
        raw, c("MarkerName", "SNP", "rsid", "variant"), "SNP"
      ),
      resolve_external_column(
        raw, c("Allele1", "A1", "effect_allele", "effectallele"),
        "effect allele"
      ),
      resolve_external_column(
        raw, c("Allele2", "A2", "other_allele", "otherallele"),
        "other allele"
      ),
      resolve_external_column(
        raw, c("Freq1", "EAF", "effect_allele_frequency", "AF"),
        "effect-allele frequency"
      ),
      resolve_external_column(
        raw, c("Effect", "BETA", "beta", "b"), "effect estimate"
      ),
      resolve_external_column(
        raw, c("StdErr", "SE", "standard_error", "stderr"),
        "standard error"
      ),
      resolve_external_column(
        raw, c("P-value", "P", "PVAL", "pvalue"), "P value"
      ),
      resolve_external_column(
        raw, c("TotalSampleSize", "N", "samplesize", "sample_size"),
        "sample size", required = FALSE
      )
    )
  )

  source_for <- function(standard_name) {
    column_map[
      standard_column == standard_name,
      source_column
    ][1L]
  }

  n_source <- source_for("N")
  n_values <- if (is.na(n_source)) {
    rep(as.numeric(fixed_n), nrow(raw))
  } else {
    suppressWarnings(as.numeric(raw[[n_source]]))
  }
  n_values[!is.finite(n_values) | n_values <= 0] <- as.numeric(fixed_n)

  dt <- data.table::data.table(
    SNP = trimws(as.character(raw[[source_for("SNP")]])),
    effect_allele = toupper(trimws(
      as.character(raw[[source_for("effect_allele")]])
    )),
    other_allele = toupper(trimws(
      as.character(raw[[source_for("other_allele")]])
    )),
    eaf = suppressWarnings(as.numeric(raw[[source_for("eaf")]])),
    beta = suppressWarnings(as.numeric(raw[[source_for("beta")]])),
    se = suppressWarnings(as.numeric(raw[[source_for("se")]])),
    pval = suppressWarnings(as.numeric(raw[[source_for("pval")]])),
    N = n_values
  )

  dt[
    (!is.finite(pval) | is.na(pval)) &
      is.finite(beta) & is.finite(se) & se > 0,
    pval := normal_pvalue(beta, se)
  ]
  dt[, maf := pmin(eaf, 1 - eaf)]

  valid_rsid <- grepl("^rs[0-9]+$", dt$SNP, ignore.case = TRUE)
  valid_alleles <- dt$effect_allele %in% c("A", "C", "G", "T") &
    dt$other_allele %in% c("A", "C", "G", "T") &
    dt$effect_allele != dt$other_allele
  valid_numeric <- is.finite(dt$beta) & is.finite(dt$se) & dt$se > 0 &
    is.finite(dt$pval) & dt$pval >= 0 & dt$pval <= 1 &
    is.finite(dt$eaf) & dt$eaf > 0 & dt$eaf < 1

  n_valid_before_dedup <- sum(
    valid_rsid & valid_alleles & valid_numeric,
    na.rm = TRUE
  )
  filtered <- dt[valid_rsid & valid_alleles & valid_numeric]
  data.table::setorder(filtered, SNP, pval, -N)
  filtered <- unique(filtered, by = "SNP")

  safe_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    fun(x)
  }

  input_qc <- data.table::data.table(
    metric = c(
      "dataset",
      "input_file",
      "raw_rows",
      "rows_with_valid_rsid",
      "rows_with_valid_biallelic_snp_alleles",
      "rows_with_valid_effect_se_p_and_eaf",
      "valid_rows_before_rsid_deduplication",
      "unique_valid_rsid_rows",
      "duplicate_rsid_rows_removed",
      "median_per_variant_sample_size",
      "minimum_per_variant_sample_size",
      "maximum_per_variant_sample_size",
      "genome_build",
      "effect_allele_definition"
    ),
    value = as.character(c(
      config$external_alm$dataset_label,
      normalizePath(path, winslash = "/", mustWork = TRUE),
      nrow(raw),
      sum(valid_rsid, na.rm = TRUE),
      sum(valid_alleles, na.rm = TRUE),
      sum(valid_numeric, na.rm = TRUE),
      n_valid_before_dedup,
      nrow(filtered),
      n_valid_before_dedup - nrow(filtered),
      safe_summary(filtered$N, stats::median),
      safe_summary(filtered$N, min),
      safe_summary(filtered$N, max),
      "GRCh37 / hg19",
      "Allele1; Effect and Freq1 refer to Allele1"
    ))
  )

  list(gwas = filtered, input_qc = input_qc, column_map = column_map)
}

format_external_alm_outcome <- function(gwas, snps) {
  sub <- gwas[SNP %in% snps]
  if (nrow(sub) == 0L) return(NULL)

  out <- tryCatch(
    TwoSampleMR::format_data(
      as.data.frame(sub),
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval"
    ),
    error = function(e) {
      warning("GEFOS ALM outcome formatting failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(out) || nrow(out) == 0L) return(NULL)

  out$outcome <- config$external_alm$dataset_label
  out$id.outcome <- config$external_alm$dataset_label
  out$samplesize.outcome <- sub$N[match(out$SNP, sub$SNP)]
  out$samplesize.outcome[
    !is.finite(out$samplesize.outcome) | out$samplesize.outcome <= 0
  ] <- config$external_alm$sample_size
  out
}

prepare_external_alm_targets <- function() {
  target_genes <- config$external_alm$target_genes
  if (anyDuplicated(target_genes)) {
    stop("config$external_alm$target_genes contains duplicated symbols.")
  }
  if (nrow(coloc_primary) == 0L) {
    stop(
      "GEFOS ALM replication requires the completed primary ",
      "colocalisation results."
    )
  }

  alm_tier_a <- coloc_primary[
    trait == "ALM" & evidence_tier == "Tier_A"
  ]
  missing_tier_a <- setdiff(target_genes, alm_tier_a$gene_symbol)
  if (config$external_alm$require_primary_tier_a &&
      length(missing_tier_a) > 0L) {
    stop(
      "The fixed external-ALM target list is inconsistent with the primary ",
      "ALM Tier A results. Missing Tier A genes: ",
      paste(missing_tier_a, collapse = ", "),
      ". Do not silently change the target list after seeing external results."
    )
  }

  target_coloc <- unique(
    alm_tier_a[
      gene_symbol %in% target_genes,
      .(
        gene_ensembl, gene_symbol,
        discovery_PP.H3 = PP.H3,
        discovery_PP.H4 = PP.H4,
        discovery_evidence_tier = evidence_tier
      )
    ],
    by = c("gene_ensembl", "gene_symbol")
  )
  target_mr <- unique(
    mr_main[
      trait == "ALM" & gene_symbol %in% target_genes,
      .(
        gene_ensembl, gene_symbol,
        discovery_method = method,
        discovery_nsnp = nsnp,
        discovery_beta = beta,
        discovery_se = se,
        discovery_pval = pval,
        discovery_fdr = fdr_within_trait
      )
    ],
    by = c("gene_ensembl", "gene_symbol")
  )
  target_positions <- data.table::as.data.table(gene_positions)[
    gene_symbol %in% target_genes,
    .(gene_ensembl, gene_symbol, chr, start, end)
  ]

  targets <- Reduce(
    function(x, y) {
      merge(
        x, y,
        by = c("gene_ensembl", "gene_symbol"),
        all = FALSE,
        sort = FALSE
      )
    },
    list(target_coloc, target_mr, target_positions)
  )
  targets <- unique(targets, by = "gene_symbol")
  targets[, target_order := match(gene_symbol, target_genes)]
  data.table::setorder(targets, target_order)

  missing_targets <- setdiff(target_genes, targets$gene_symbol)
  if (length(missing_targets) > 0L) {
    stop(
      "Could not fully define all four GEFOS ALM replication targets: ",
      paste(missing_targets, collapse = ", ")
    )
  }
  targets
}

run_external_alm_replication <- function(targets, external_gwas) {
  result_rows <- vector("list", nrow(targets))
  harmonised_rows <- vector("list", nrow(targets))

  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    log_step(
      "Running targeted GEFOS ALM replication: ", target$gene_symbol
    )

    result_row <- data.table::data.table(
      target_order = target$target_order,
      gene_ensembl = target$gene_ensembl,
      gene_symbol = target$gene_symbol,
      discovery_outcome = "ALM",
      external_outcome = config$external_alm$dataset_label,
      discovery_method = target$discovery_method,
      discovery_nsnp = target$discovery_nsnp,
      discovery_beta = target$discovery_beta,
      discovery_se = target$discovery_se,
      discovery_pval = target$discovery_pval,
      discovery_fdr = target$discovery_fdr,
      discovery_PP.H3 = target$discovery_PP.H3,
      discovery_PP.H4 = target$discovery_PP.H4,
      n_eqtl_region = 0L,
      n_strong_instruments = 0L,
      n_clumped_instruments = 0L,
      n_external_outcome_snps = 0L,
      n_harmonised_instruments = 0L,
      external_method = NA_character_,
      external_nsnp = 0L,
      external_beta = NA_real_,
      external_se = NA_real_,
      external_ci_lower = NA_real_,
      external_ci_upper = NA_real_,
      external_pval = NA_real_,
      minimum_F = NA_real_,
      mean_F = NA_real_,
      test_status = "not_tested",
      exclusion_reason = NA_character_
    )

    eqtl <- load_gene_eqtl(
      target$gene_ensembl, target$chr, target$start, target$end
    )
    result_row[, n_eqtl_region := nrow(eqtl)]
    exposure <- format_exposure(eqtl, target$gene_symbol)
    result_row[, n_strong_instruments :=
      if (is.null(exposure)) 0L else nrow(exposure)]
    exposure <- clump_exposure(exposure)
    result_row[, n_clumped_instruments :=
      if (is.null(exposure)) 0L else nrow(exposure)]

    if (is.null(exposure) || nrow(exposure) == 0L) {
      result_row[, exclusion_reason :=
        "no strong independent skeletal-muscle cis-eQTL instrument"]
      result_rows[[i]] <- result_row
      next
    }

    result_row[, n_external_outcome_snps :=
      sum(exposure$SNP %in% external_gwas$SNP)]
    outcome <- format_external_alm_outcome(external_gwas, exposure$SNP)
    if (is.null(outcome) || nrow(outcome) == 0L) {
      result_row[, exclusion_reason :=
        "none of the prespecified instruments was available in GEFOS"]
      result_rows[[i]] <- result_row
      next
    }

    harmonisation_error <- NA_character_
    dat_all <- tryCatch(
      TwoSampleMR::harmonise_data(exposure, outcome, action = 3),
      error = function(e) {
        harmonisation_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(dat_all) || nrow(dat_all) == 0L) {
      result_row[, exclusion_reason := paste0(
        "allele harmonisation failed",
        ifelse(
          is.na(harmonisation_error),
          "",
          paste0(": ", harmonisation_error)
        )
      )]
      result_rows[[i]] <- result_row
      next
    }

    dat_audit <- data.table::as.data.table(dat_all)
    dat_audit[, `:=`(
      gene_ensembl = target$gene_ensembl,
      gene_symbol = target$gene_symbol,
      discovery_outcome = "ALM",
      external_outcome = config$external_alm$dataset_label,
      included_in_external_mr = mr_keep %in% TRUE
    )]
    harmonised_rows[[i]] <- dat_audit

    dat <- dat_all[dat_all$mr_keep %in% TRUE, , drop = FALSE]
    result_row[, n_harmonised_instruments := nrow(dat)]
    if (nrow(dat) == 0L) {
      result_row[, exclusion_reason :=
        "all overlapping instruments were removed during harmonisation"]
      result_rows[[i]] <- result_row
      next
    }

    method <- if (nrow(dat) == 1L) "mr_wald_ratio" else "mr_ivw"
    estimate <- tryCatch(
      TwoSampleMR::mr(dat, method_list = method),
      error = function(e) {
        warning(
          "GEFOS ALM MR failed for ", target$gene_symbol, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(estimate) || nrow(estimate) == 0L) {
      result_row[, exclusion_reason := "MR estimation failed"]
      result_rows[[i]] <- result_row
      next
    }

    result_row[, `:=`(
      external_method = as.character(estimate$method[1L]),
      external_nsnp = as.integer(estimate$nsnp[1L]),
      external_beta = as.numeric(estimate$b[1L]),
      external_se = as.numeric(estimate$se[1L]),
      external_ci_lower =
        as.numeric(estimate$b[1L]) - stats::qnorm(0.975) *
          as.numeric(estimate$se[1L]),
      external_ci_upper =
        as.numeric(estimate$b[1L]) + stats::qnorm(0.975) *
          as.numeric(estimate$se[1L]),
      external_pval = as.numeric(estimate$pval[1L]),
      minimum_F = min(
        (dat$beta.exposure / dat$se.exposure)^2,
        na.rm = TRUE
      ),
      mean_F = mean(
        (dat$beta.exposure / dat$se.exposure)^2,
        na.rm = TRUE
      ),
      test_status = "tested",
      exclusion_reason = NA_character_
    )]
    result_rows[[i]] <- result_row
  }

  results <- bind_nonempty(result_rows)
  instruments <- bind_nonempty(harmonised_rows)
  n_prespecified <- length(config$external_alm$target_genes)

  if (nrow(results) > 0L) {
    results[, external_p_bonferroni := stats::p.adjust(
      external_pval, method = "bonferroni", n = n_prespecified
    )]
    results[, external_fdr_4gene := stats::p.adjust(
      external_pval, method = "BH", n = n_prespecified
    )]
    results[, direction_concordant :=
      data.table::fifelse(
        test_status == "tested" &
          is.finite(discovery_beta) & is.finite(external_beta),
        sign(discovery_beta) == sign(external_beta),
        NA
      )]
    results[, nominal_directional_support :=
      direction_concordant %in% TRUE &
        external_pval < config$external_alm$nominal_alpha]
    results[, bonferroni_directional_support :=
      direction_concordant %in% TRUE &
        external_pval < config$external_alm$bonferroni_alpha]
    results[, replication_status := data.table::fcase(
      test_status != "tested",
      "not evaluable",
      bonferroni_directional_support,
      "directionally concordant and Bonferroni-supported",
      nominal_directional_support,
      "directionally concordant and nominally supported",
      direction_concordant %in% TRUE,
      "directionally concordant without nominal support",
      is.finite(external_pval) &
        external_pval < config$external_alm$nominal_alpha,
      "nominal association in the opposite direction",
      default = "opposite direction without nominal support"
    )]
    data.table::setorder(results, target_order)
  }

  list(results = results, instruments = instruments)
}

external_alm_validation <- list(
  results = empty_dt(),
  instruments = empty_dt(),
  targets = empty_dt(),
  input_qc = empty_dt(),
  column_map = empty_dt()
)

resolve_external_alm_file <- function(configured_path) {
  file_name <- basename(configured_path)
  candidates <- unique(c(
    configured_path,
    file.path(config$project_dir, file_name),
    file.path(config$project_dir, "outcome", file_name),
    file.path(config$project_dir, "data", "raw", "outcome", file_name)
  ))
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) return(configured_path)

  hits <- unique(normalizePath(hits, winslash = "/", mustWork = TRUE))
  if (length(hits) > 1L) {
    stop(
      "Multiple copies of the GEFOS ALM file were found. Set ",
      "config$files$external_alm to one explicit path:\n",
      paste(hits, collapse = "\n")
    )
  }
  hits[1L]
}

if (config$external_alm$run) {
  config$files$external_alm <- resolve_external_alm_file(
    config$files$external_alm
  )
  assert_files_exist(config$files$external_alm, "GEFOS 2017 ALM GWAS")
  external_alm_input <- read_gefos_alm(config$files$external_alm)
  external_alm_targets <- prepare_external_alm_targets()
  external_alm_fit <- run_external_alm_replication(
    external_alm_targets,
    external_alm_input$gwas
  )
  external_alm_validation <- c(
    external_alm_fit,
    list(
      targets = external_alm_targets,
      input_qc = external_alm_input$input_qc,
      column_map = external_alm_input$column_map
    )
  )
}

external_alm_readme <- data.table::data.table(
  item = c(
    "Purpose",
    "External outcome",
    "Input file",
    "Fixed discovery-derived targets",
    "Instrument definition",
    "Outcome harmonisation",
    "Primary supportive criterion",
    "Strict supportive criterion",
    "Multiple-testing family",
    "Analyses intentionally not run",
    "Interpretation boundary",
    "Citation"
  ),
  value = c(
    "Limited external outcome support for four discovery-derived ALM targets using the original GTEx instrument rules.",
    "GEFOS 2017 appendicular lean mass discovery meta-analysis; European ancestry; N up to 28,330; independent of UK Biobank.",
    config$files$external_alm,
    paste(config$external_alm$target_genes, collapse = ", "),
    "The same GTEx v8 skeletal-muscle cis-eQTL P < 5e-8, F >= 10 and EUR LD-clumping rules used in the primary analysis.",
    "rsID and alleles; harmonise_data action = 3; no proxy replacement.",
    "MR effect direction concordant with discovery and two-sided P < 0.05.",
    paste0(
      "MR effect direction concordant with discovery and Bonferroni P < ",
      format(config$external_alm$bonferroni_alpha, scientific = FALSE),
      " across four fixed genes."
    ),
    "Exactly four prespecified ALM genes, including non-evaluable targets.",
    "No re-screening of all IR genes, no proxy-SNP search, and no meta-analysis with the discovery outcome. External colocalisation is single-signal coloc.abf sensitivity evidence and is not definitive replication.",
    "External effect magnitudes need not be directly comparable because phenotype measurement and transformation differ; direction and prespecified P-value support are emphasized. A non-significant result is inconclusive when power or instrument coverage is limited.",
    "Zillikens et al. Nature Communications. 2017;8:80. doi:10.1038/s41467-017-00031-7; GWAS Catalog GCST005037."
  )
)

write_result(external_alm_readme, "09_GEFOS_ALM_README.csv")
write_result(
  external_alm_validation$column_map,
  "09_GEFOS_ALM_column_map.csv"
)
write_result(
  external_alm_validation$input_qc,
  "09_GEFOS_ALM_input_QC.csv"
)
write_result(
  external_alm_validation$targets,
  "09_GEFOS_ALM_targets.csv"
)
write_result(
  external_alm_validation$results,
  "09_GEFOS_ALM_replication.csv"
)
write_result(
  external_alm_validation$instruments,
  "09_GEFOS_ALM_instruments.csv"
)


# [NEW: Reviewer 7] External outcome coloc, using the same regional allele checks
# and all three priors. This is single-signal sensitivity evidence, not SuSiE.
config$outcome_type <- c(config$outcome_type, GEFOS_2017_ALM = "quant")
external_coloc_rows <- lapply(seq_len(nrow(external_alm_validation$targets)), function(i) {
  target <- external_alm_validation$targets[i]
  eqtl <- load_gene_eqtl(target$gene_ensembl, target$chr, target$start, target$end)
  fit <- run_coloc_pair(eqtl, external_alm_input$gwas, target$gene_symbol,
                        "GEFOS_2017_ALM")
  if (!nrow(fit)) stop("External coloc produced no result for ", target$gene_symbol)
  fit
})
external_coloc_all <- bind_nonempty(external_coloc_rows)
stopifnot(nrow(external_coloc_all) == 12L)
external_evidence <- merge(external_alm_validation$results,
  external_coloc_all[p12 == config$coloc$primary_p12,
                    .(gene_symbol, external_PP3 = PP.H3, external_PP4 = PP.H4)],
  by = "gene_symbol", all.x = TRUE)
external_evidence[, evidence_boundary :=
  "Limited external outcome support; original GTEx instruments reused; external coloc assumes one signal."]
# [NEW] Explicit analysis provenance for the external coloc table.
external_coloc_all[, `:=`(
  exposure_sdY = 1,
  outcome_sdY_method = "sqrt(median(2 * MAF * (1-MAF) * per-variant N * SE^2))",
  allele_rule = "main harmonise_coloc_region; MAF >=0.01 and <0.5; ambiguous palindromes removed",
  method = "coloc.abf; same assumptions as primary pipeline"
)]
write_result(external_coloc_all, "09_GEFOS_coloc_all_priors.csv")
write_result(external_evidence, "09_GEFOS_external_evidence.csv")

## 14. Reproducibility record and final summary ------------------------------

capture.output(sessionInfo(), file = result_path("sessionInfo.txt"))
saveRDS(config, result_path("analysis_config.rds"))

summary_lines <- c(
  paste0("Candidate IR genes: ", length(candidate_genes)),
  paste0("Genes mapped to coordinates: ", nrow(gene_positions)),
  paste0("cis-MR gene-trait tests: ", nrow(mr_main)),
  paste0("FDR-significant cis-MR pairs: ", sum(mr_main$mr_associated)),
  paste0("FDR-significant genes: ", uniqueN(mr_main[mr_associated == TRUE, gene_symbol])),
  paste0(
    "Tier A colocalisation-supported genes: ",
    if (nrow(coloc_primary) > 0L)
      uniqueN(coloc_primary[evidence_tier == "Tier_A", gene_symbol]) else 0L
  ),
  "Interpretation: MR-associated genes are candidates; Tier A genes are",
  "colocalisation-supported prioritised candidates, not confirmed causal genes."
)
writeLines(summary_lines, result_path("RUN_SUMMARY.txt"))
openxlsx::write.xlsx(
  mget(ls(result_tables), envir = result_tables, inherits = FALSE),
  file = result_path("manuscript_analysis_tables.xlsx"),
  overwrite = TRUE,
  asTable = TRUE
)
log_step("Analysis complete. Outputs: ", config$paths$results)


## 15. Final workbook compatibility repair and validation -------------------

normalise_ooxml_path <- function(path) {
  path_parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
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

repair_missing_ooxml_relationships <- function(xlsx_file) {
  if (!requireNamespace("xml2", quietly = TRUE) ||
      !requireNamespace("zip", quietly = TRUE)) {
    stop("Final workbook repair requires the xml2 and zip packages.")
  }
  if (!file.exists(xlsx_file)) {
    stop("Workbook not found: ", xlsx_file)
  }

  repair_dir <- tempfile("xlsx_repair_")
  repaired_file <- tempfile(fileext = ".xlsx")
  dir.create(repair_dir, recursive = TRUE)
  on.exit(unlink(repair_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(repaired_file, force = TRUE), add = TRUE)

  utils::unzip(xlsx_file, exdir = repair_dir)
  archive_files <- list.files(
    repair_dir, recursive = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = FALSE
  )
  archive_files <- gsub("\\\\", "/", archive_files)

  relationship_files <- list.files(
    repair_dir, pattern = "\\.rels$", recursive = TRUE,
    full.names = TRUE, all.files = TRUE
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
      relationship_xml, "//*[local-name()='Relationship']"
    )
    for (relationship_node in relationship_nodes) {
      if (identical(xml2::xml_attr(relationship_node, "TargetMode"),
                    "External")) {
        next
      }
      relationship_target <- xml2::xml_attr(relationship_node, "Target")
      target_path <- if (startsWith(relationship_target, "/")) {
        sub("^/+", "", relationship_target)
      } else {
        normalise_ooxml_path(file.path(
          relationship_base, relationship_target
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
    content_types_xml, "//*[local-name()='Override']"
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

  top_level_entries <- list.files(
    repair_dir, all.files = TRUE, no.. = TRUE, full.names = FALSE
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
    stop("Failed to replace the workbook after compatibility repair.")
  }

  list(
     removed_relationships = removed_relationships,
    removed_content_types = removed_content_types
  )
}


final_workbook <- result_path("manuscript_analysis_tables.xlsx")
workbook_repair <- repair_missing_ooxml_relationships(final_workbook)
expected_sheets <- ls(result_tables)
actual_sheets <- openxlsx::getSheetNames(final_workbook)
if (!identical(actual_sheets, expected_sheets)) {
  stop("Final workbook validation failed: worksheet names do not match.")}


log_step(
  "Final workbook validated; removed ",
  workbook_repair$removed_relationships,
  " missing-part relationship(s) and ",
  workbook_repair$removed_content_types,
  " missing-part content-type declaration(s)."
)



## 16. MAGIC insulin-resistance mechanistic extension ------------------------
##
## This section is appended to the completed primary pipeline. It does not
## redefine the original candidate set or alter the primary cis-MR results.
## It tests a prespecified mechanistic chain among Tier A candidates only:
##
## skeletal-muscle gene expression -> insulin-resistance trait
## insulin-resistance trait -> sarcopenia-related trait
##
## The product-of-coefficients results are exploratory two-step MR estimates.
## They are not interpreted as proof of mediation.

log_step("Starting the appended MAGIC insulin-resistance extension")

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  stop("The appended MAGIC analysis requires the rtracklayer package.")
}

magic_config <- list(
  files = c(
    HOMA_IR = file.path(config$project_dir, "HOMA_IR.vcf.gz"),
    ISI_adjBMI = file.path(config$project_dir, "ISI_adjBMI_EUR.gz"),
    FASTING_INSULIN = file.path(config$project_dir, "FI_EUR.gz")
  ),
  grch38_to_grch37_chain = file.path(
    config$project_dir, "hg38ToHg19.over.chain.gz"
  ),
  build = "GRCh37",
  chunk_size = 100000L,
  instrument_p = 5e-8,
  minimum_F = 10,
  minimum_position_concordance = 0.90,
  output_workbook = result_path("MAGIC_IR_mechanistic_extension.xlsx")
)

assert_files_exist(
  c(unname(magic_config$files), magic_config$grch38_to_grch37_chain),
  "MAGIC extension input"
)

if (nrow(coloc_primary) == 0L ||
    !"evidence_tier" %in% names(coloc_primary)) {
  stop("No primary colocalisation results are available for MAGIC extension.")
}

magic_targets <- unique(
  coloc_primary[evidence_tier == "Tier_A", .(
    gene_ensembl, gene_symbol, gene_layer
  )]
)
if (nrow(magic_targets) == 0L) {
  stop("No Tier A candidate genes are available for MAGIC extension.")
}

import_compressed_chain <- function(path) {
  temporary_chain <- tempfile(fileext = ".chain")
  input_connection <- gzfile(path, open = "rt")
  output_connection <- file(temporary_chain, open = "wt")
  on.exit({
    try(close(input_connection), silent = TRUE)
    try(close(output_connection), silent = TRUE)
    unlink(temporary_chain, force = TRUE)
  }, add = TRUE)

  repeat {
    chain_lines <- readLines(input_connection, n = 100000L, warn = FALSE)
    if (length(chain_lines) == 0L) break
    writeLines(chain_lines, output_connection)
  }
  close(input_connection)
  close(output_connection)
  rtracklayer::import.chain(temporary_chain)
}

complement_allele <- function(x) {
  chartr("ACGT", "TGCA", toupper(as.character(x)))
}

variant_pair_key <- function(chr, pos, allele1, allele2) {
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  paste(
    as.character(chr), as.integer(pos),
    pmin(allele1, allele2), pmax(allele1, allele2), sep = ":"
  )
}

log_step("Lifting Tier A skeletal-muscle eQTL coordinates to GRCh37")
magic_chain <- import_compressed_chain(magic_config$grch38_to_grch37_chain)

magic_eqtl_by_gene <- setNames(vector("list", nrow(magic_targets)),
                               magic_targets$gene_ensembl)
magic_eqtl_map_rows <- vector("list", nrow(magic_targets))

for (i in seq_len(nrow(magic_targets))) {
  target <- magic_targets[i]
  gp <- gene_positions[
    gene_positions$gene_ensembl == target$gene_ensembl,
    , drop = FALSE
  ]
  if (nrow(gp) != 1L) next

  eqtl <- load_gene_eqtl(
    gp$gene_ensembl, gp$chr, gp$start, gp$end
  )
  if (nrow(eqtl) == 0L) next
  eqtl <- unique(eqtl, by = "SNP")
  magic_eqtl_by_gene[[target$gene_ensembl]] <- eqtl

  query <- GenomicRanges::GRanges(
    seqnames = paste0("chr", eqtl$chr),
    ranges = IRanges::IRanges(eqtl$pos, width = 1L)
  )
  lifted <- rtracklayer::liftOver(query, magic_chain)
  unique_mapping <- lengths(lifted) == 1L
  chr37 <- rep(NA_character_, nrow(eqtl))
  pos37 <- rep(NA_integer_, nrow(eqtl))
  chr37[unique_mapping] <- sub(
    "^chr", "",
    as.character(GenomicRanges::seqnames(unlist(lifted[unique_mapping])))
  )
  pos37[unique_mapping] <- GenomicRanges::start(
    unlist(lifted[unique_mapping])
  )

  magic_eqtl_map_rows[[i]] <- data.table(
    gene_ensembl = target$gene_ensembl,
    gene_symbol = target$gene_symbol,
    SNP = eqtl$SNP,
    chr38 = as.character(eqtl$chr),
    pos38 = as.integer(eqtl$pos),
    chr37 = chr37,
    pos37 = pos37,
    effect_allele_eqtl = eqtl$effect_allele,
    other_allele_eqtl = eqtl$other_allele,
    eaf_eqtl = eqtl$eaf
  )
}

magic_eqtl_map <- bind_nonempty(magic_eqtl_map_rows)
magic_eqtl_map <- magic_eqtl_map[
  !is.na(chr37) & !is.na(pos37) & chr37 %in% as.character(1:22)
]
if (nrow(magic_eqtl_map) == 0L) {
  stop("No Tier A eQTL variants could be lifted from GRCh38 to GRCh37.")
}

magic_regions <- magic_eqtl_map[, .(
  region_start = min(pos37),
  region_end = max(pos37)
), by = .(gene_ensembl, gene_symbol, chr37)]

inside_magic_regions <- function(chr, pos, regions = magic_regions) {
  keep <- rep(FALSE, length(pos))
  for (i in seq_len(nrow(regions))) {
    keep <- keep |
      (as.character(chr) == regions$chr37[i] &
         as.integer(pos) >= regions$region_start[i] &
         as.integer(pos) <= regions$region_end[i])
  }
  keep
}

standardise_magic_rows <- function(dt, trait) {
  if (is.null(dt) || nrow(dt) == 0L) return(empty_dt())
  dt <- data.table::as.data.table(dt)
  for (column in c("SNP", "eaf")) {
    if (!column %in% names(dt)) dt[, (column) := NA]
  }
  dt[, `:=`(
    trait = trait,
    SNP = as.character(SNP),
    chr = sub("^chr", "", as.character(chr), ignore.case = TRUE),
    pos = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele)),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pval = as.numeric(pval),
    eaf = as.numeric(eaf),
    N = as.numeric(N)
  )]
  dt[, maf := pmin(eaf, 1 - eaf)]
  dt[
    chr %in% as.character(1:22) & is.finite(pos) &
      is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(pval) & pval >= 0 & pval <= 1 &
      effect_allele %in% c("A", "C", "G", "T") &
      other_allele %in% c("A", "C", "G", "T")
  ]
}

read_magic_table_stream <- function(path, trait, column_spec,
                                    target_snps = character()) {
  log_step("Streaming MAGIC table: ", trait)
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1]]
  file_columns <- unname(column_spec)
  selected <- match(file_columns, header)
  if (anyNA(selected)) {
    stop(
      trait, " is missing required column(s): ",
      paste(file_columns[is.na(selected)], collapse = ", ")
    )
  }

  regional_rows <- list()
  instrument_rows <- list()
  total_rows <- 0L
  chunk_index <- 0L

  repeat {
    lines <- readLines(
      connection, n = magic_config$chunk_size, warn = FALSE
    )
    if (length(lines) == 0L) break
    chunk_index <- chunk_index + 1L
    total_rows <- total_rows + length(lines)

    chunk <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE, sep = "\t", select = selected,
      showProgress = FALSE
    )
    data.table::setnames(chunk, names(column_spec))
    chunk <- standardise_magic_rows(chunk, trait)
    if (nrow(chunk) == 0L) next

    in_region <- inside_magic_regions(chunk$chr, chunk$pos)
    if (length(target_snps) > 0L) {
      in_region <- in_region | (!is.na(chunk$SNP) &
                                  chunk$SNP %in% target_snps)
    }
    is_instrument <- chunk$pval < magic_config$instrument_p &
      (chunk$beta / chunk$se)^2 >= magic_config$minimum_F

    if (any(in_region)) {
      regional_rows[[length(regional_rows) + 1L]] <- chunk[in_region]
    }
    if (any(is_instrument)) {
      instrument_rows[[length(instrument_rows) + 1L]] <- chunk[is_instrument]
    }
    if (chunk_index %% 25L == 0L) {
      log_step(trait, ": scanned ", format(total_rows, big.mark = ","),
               " variants")
    }
  }

  list(
    regional = unique(bind_nonempty(regional_rows),
                      by = c("chr", "pos", "effect_allele", "other_allele")),
    instruments = unique(bind_nonempty(instrument_rows),
                         by = c("chr", "pos", "effect_allele", "other_allele")),
    total_rows = total_rows
  )
}

read_magic_homa_stream <- function(path, trait, target_snps = character()) {
  log_step("Streaming MAGIC VCF: ", trait)
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    header_line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header_line) == 0L) stop("VCF header not found: ", path)
    if (startsWith(header_line, "#CHROM")) break
  }

  regional_rows <- list()
  instrument_rows <- list()
  total_rows <- 0L
  chunk_index <- 0L

  repeat {
    lines <- readLines(
      connection, n = magic_config$chunk_size, warn = FALSE
    )
    if (length(lines) == 0L) break
    chunk_index <- chunk_index + 1L
    total_rows <- total_rows + length(lines)

    chunk <- data.table::fread(
      text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t",
      select = c(1L, 2L, 3L, 4L, 5L, 9L, 10L),
      col.names = c(
        "chr", "pos", "SNP", "other_allele", "effect_allele",
        "FORMAT", "sample"
      ),
      showProgress = FALSE
    )
    chunk[, `:=`(beta = NA_real_, se = NA_real_, pval = NA_real_,
                 N = NA_real_, eaf = NA_real_)]

    for (format_template in unique(chunk$FORMAT)) {
      row_index <- which(chunk$FORMAT == format_template)
      format_fields <- strsplit(format_template, ":", fixed = TRUE)[[1]]
      format_values <- data.table::tstrsplit(
        chunk$sample[row_index], ":", fixed = TRUE
      )
      get_value <- function(field) {
        field_index <- match(field, format_fields)
        if (is.na(field_index)) return(rep(NA_real_, length(row_index)))
        as.numeric(format_values[[field_index]])
      }
      chunk[row_index, `:=`(
        beta = get_value("ES"),
        se = get_value("SE"),
        pval = 10^(-get_value("LP")),
        N = get_value("SS"),
        eaf = get_value("AF")
      )]
    }
    chunk[, c("FORMAT", "sample") := NULL]
    chunk <- standardise_magic_rows(chunk, trait)
    if (nrow(chunk) == 0L) next

    in_region <- inside_magic_regions(chunk$chr, chunk$pos) |
      (!is.na(chunk$SNP) & chunk$SNP %in% target_snps)
    is_instrument <- chunk$pval < magic_config$instrument_p &
      (chunk$beta / chunk$se)^2 >= magic_config$minimum_F

    if (any(in_region)) {
      regional_rows[[length(regional_rows) + 1L]] <- chunk[in_region]
    }
    if (any(is_instrument)) {
      instrument_rows[[length(instrument_rows) + 1L]] <- chunk[is_instrument]
    }
    if (chunk_index %% 25L == 0L) {
      log_step(trait, ": scanned ", format(total_rows, big.mark = ","),
               " variants")
    }
  }

  list(
    regional = unique(bind_nonempty(regional_rows), by = "SNP"),
    instruments = unique(bind_nonempty(instrument_rows), by = "SNP"),
    total_rows = total_rows
  )
}

target_magic_snps <- unique(magic_eqtl_map$SNP)

magic_raw <- list(
  HOMA_IR = read_magic_homa_stream(
    magic_config$files[["HOMA_IR"]], "HOMA_IR", target_magic_snps
  ),
  ISI_adjBMI = read_magic_table_stream(
    magic_config$files[["ISI_adjBMI"]], "ISI_adjBMI",
    c(
      chr = "chromosome", pos = "base_pair_location",
      effect_allele = "effect_allele", other_allele = "other_allele",
      beta = "beta", se = "standard_error", eaf = "effect_allele_frequency",
      pval = "p_value", N = "n", SNP = "variant_id"
    )
  ),
  FASTING_INSULIN = read_magic_table_stream(
    magic_config$files[["FASTING_INSULIN"]], "FASTING_INSULIN",
    c(
      SNP = "variant", chr = "chromosome", pos = "base_pair_location",
      effect_allele = "effect_allele", other_allele = "other_allele",
      eaf = "effect_allele_frequency", beta = "beta",
      se = "standard_error", pval = "p_value", N = "sample_size"
    ),
    target_magic_snps
  )
)

map_regional_magic_variants <- function(ir, trait, eqtl_map) {
  if (nrow(ir) == 0L) return(empty_dt())

  if (trait != "ISI_adjBMI") {
    mapped <- merge(
      ir[!is.na(SNP) & grepl("^rs", SNP)],
      unique(eqtl_map[, .(SNP, pos37)]),
      by = "SNP", all = FALSE, allow.cartesian = FALSE
    )
    if (nrow(mapped) == 0L) return(empty_dt())
    mapped[, position_concordant := pos == pos37]
    mapped[, pos37 := NULL]
    return(unique(mapped, by = "SNP"))
  }

  direct <- unique(eqtl_map[, .(
    lookup_key = variant_pair_key(
      chr37, pos37, effect_allele_eqtl, other_allele_eqtl
    ),
    SNP,
    strand_match = "direct"
  )])
  complemented <- unique(eqtl_map[, .(
    lookup_key = variant_pair_key(
      chr37, pos37,
      complement_allele(effect_allele_eqtl),
      complement_allele(other_allele_eqtl)
    ),
    SNP,
    strand_match = "complement"
  )])
  aliases <- unique(rbind(direct, complemented), by = "lookup_key")
  ir <- copy(ir)
  ir[, lookup_key := variant_pair_key(
    chr, pos, effect_allele, other_allele
  )]
  ir[, SNP := NULL]
  mapped <- merge(ir, aliases, by = "lookup_key", all = FALSE)
  mapped[strand_match == "complement", `:=`(
    effect_allele = complement_allele(effect_allele),
    other_allele = complement_allele(other_allele)
  )]
  mapped[, c("lookup_key", "strand_match") := NULL]
  mapped[, position_concordant := TRUE]
  unique(mapped, by = "SNP")
}

magic_regional <- lapply(names(magic_raw), function(trait) {
  map_regional_magic_variants(
    magic_raw[[trait]]$regional, trait, magic_eqtl_map
  )
})
names(magic_regional) <- names(magic_raw)

position_qc <- bind_nonempty(lapply(names(magic_regional), function(trait) {
  x <- magic_regional[[trait]]
  data.table(
    trait = trait,
    n_mapped_eqtl_variants = nrow(x),
    n_position_checked = sum(!is.na(x$position_concordant)),
    position_concordance = if (nrow(x) > 0L)
      mean(x$position_concordant, na.rm = TRUE) else NA_real_
  )
}))

if (any(position_qc$n_mapped_eqtl_variants < config$coloc$minimum_common_snps)) {
  warning(
    "At least one MAGIC trait has fewer than ",
    config$coloc$minimum_common_snps,
    " mapped variants across all Tier A loci; some colocalisation tests may be unavailable."
  )
}
if (any(
  position_qc$n_position_checked >= 20L &
    position_qc$position_concordance <
      magic_config$minimum_position_concordance,
  na.rm = TRUE
)) {
  stop("MAGIC coordinate validation failed; expected GRCh37 coordinates.")
}

format_magic_outcome <- function(ir, eqtl, snps, trait) {
  sub <- copy(ir[SNP %in% snps])
  if (nrow(sub) == 0L) return(NULL)
  reference_eaf <- eqtl$eaf[match(sub$SNP, eqtl$SNP)]
  missing_eaf <- !is.finite(sub$eaf) | sub$eaf <= 0 | sub$eaf >= 1
  sub[missing_eaf, eaf := reference_eaf[missing_eaf]]
  sub <- sub[is.finite(eaf) & eaf > 0 & eaf < 1]
  if (nrow(sub) == 0L) return(NULL)

  out <- tryCatch(
    TwoSampleMR::format_data(
      as.data.frame(sub), type = "outcome",
      snp_col = "SNP", beta_col = "beta", se_col = "se",
      eaf_col = "eaf", effect_allele_col = "effect_allele",
      other_allele_col = "other_allele", pval_col = "pval",
      chr_col = "chr", pos_col = "pos", samplesize_col = "N"
    ),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out$outcome <- trait
  out$id.outcome <- trait
  out
}

run_gene_to_magic_mr <- function(targets, regional_data) {
  rows <- list()
  harmonised <- list()

  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    eqtl <- magic_eqtl_by_gene[[target$gene_ensembl]]
    if (is.null(eqtl) || nrow(eqtl) == 0L) next
    exposure <- format_exposure(eqtl, target$gene_symbol)
    exposure <- clump_exposure(exposure)
    if (is.null(exposure) || nrow(exposure) == 0L) next

    for (trait in names(regional_data)) {
      outcome <- format_magic_outcome(
        regional_data[[trait]], eqtl, exposure$SNP, trait
      )
      if (is.null(outcome)) next
      dat <- tryCatch(
        TwoSampleMR::harmonise_data(exposure, outcome, action = 3),
        error = function(e) NULL
      )
      if (is.null(dat)) next
      dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
      if (nrow(dat) == 0L) next

      method <- if (nrow(dat) == 1L) "mr_wald_ratio" else "mr_ivw"
      estimate <- tryCatch(
        TwoSampleMR::mr(dat, method_list = method),
        error = function(e) NULL
      )
      if (is.null(estimate) || nrow(estimate) == 0L) next

      key <- paste(target$gene_ensembl, trait, sep = "__")
      harmonised[[key]] <- dat
      rows[[key]] <- data.table(
        gene_ensembl = target$gene_ensembl,
        gene_symbol = target$gene_symbol,
        gene_layer = target$gene_layer,
        ir_trait = trait,
        method = estimate$method[1],
        nsnp = estimate$nsnp[1],
        beta = estimate$b[1],
        se = estimate$se[1],
        pval = estimate$pval[1],
        min_F = min((dat$beta.exposure / dat$se.exposure)^2),
        mean_F = mean((dat$beta.exposure / dat$se.exposure)^2)
      )
    }
  }
  list(mr = bind_nonempty(rows), harmonised = harmonised)
}

log_step("Running Tier A gene-expression to MAGIC IR cis-MR")
gene_magic <- run_gene_to_magic_mr(magic_targets, magic_regional)
gene_to_ir_mr <- gene_magic$mr
if (nrow(gene_to_ir_mr) > 0L) {
  gene_to_ir_mr[, fdr_within_ir_trait := p.adjust(pval, "BH"),
                by = ir_trait]
  gene_to_ir_mr[, fdr_global := p.adjust(pval, "BH")]
  gene_to_ir_mr[, `:=`(
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    mr_fdr_significant = fdr_within_ir_trait < config$multiple_testing$alpha,
    interpretation_scale = fifelse(
      ir_trait == "ISI_adjBMI",
      "higher value indicates greater insulin sensitivity",
      "higher value indicates greater insulin resistance or insulin level"
    )
  )]
}

run_magic_coloc_pair <- function(eqtl, ir, gene_symbol, trait) {
  ir <- copy(ir)
  reference_eaf <- eqtl$eaf[match(ir$SNP, eqtl$SNP)]
  missing_eaf <- !is.finite(ir$eaf) | ir$eaf <= 0 | ir$eaf >= 1
  ir[missing_eaf, eaf := reference_eaf[missing_eaf]]
  ir[, maf := pmin(eaf, 1 - eaf)]
  aligned <- harmonise_coloc_region(eqtl, ir)
  if (is.null(aligned)) return(empty_dt())
  e <- aligned$eqtl
  g <- aligned$gwas

  dataset_eqtl <- list(
    snp = e$SNP, beta = e$beta, varbeta = e$se^2,
    MAF = e$maf, N = config$gtex$sample_size,
    type = "quant", sdY = 1
  )
  sdY_ir <- sqrt(stats::median(
    2 * g$maf * (1 - g$maf) * g$N * g$se^2,
    na.rm = TRUE
  ))
  if (!is.finite(sdY_ir) || sdY_ir <= 0) sdY_ir <- 1
  dataset_ir <- list(
    snp = g$SNP, beta = g$beta, varbeta = g$se^2,
    MAF = g$maf, N = stats::median(g$N, na.rm = TRUE),
    type = "quant", sdY = sdY_ir
  )

  bind_nonempty(lapply(config$coloc$p12_values, function(p12) {
    fit <- tryCatch(
      coloc::coloc.abf(
        dataset1 = dataset_eqtl, dataset2 = dataset_ir,
        p1 = config$coloc$p1, p2 = config$coloc$p2, p12 = p12
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    summary <- fit$summary
    data.table(
      gene_symbol = gene_symbol,
      ir_trait = trait,
      p1 = config$coloc$p1,
      p2 = config$coloc$p2,
      p12 = p12,
      n_snps = nrow(e),
      PP.H0 = unname(summary["PP.H0.abf"]),
      PP.H1 = unname(summary["PP.H1.abf"]),
      PP.H2 = unname(summary["PP.H2.abf"]),
      PP.H3 = unname(summary["PP.H3.abf"]),
      PP.H4 = unname(summary["PP.H4.abf"])
    )
  }))
}

log_step("Running Tier A gene-expression to MAGIC IR colocalisation")
gene_to_ir_coloc_rows <- list()
for (i in seq_len(nrow(magic_targets))) {
  target <- magic_targets[i]
  eqtl <- magic_eqtl_by_gene[[target$gene_ensembl]]
  if (is.null(eqtl) || nrow(eqtl) == 0L) next
  for (trait in names(magic_regional)) {
    fit <- run_magic_coloc_pair(
      eqtl, magic_regional[[trait]], target$gene_symbol, trait
    )
    if (nrow(fit) == 0L) next
    fit[, `:=`(
      gene_ensembl = target$gene_ensembl,
      gene_layer = target$gene_layer
    )]
    gene_to_ir_coloc_rows[[paste(target$gene_ensembl, trait, sep = "__")]] <- fit
  }
}
gene_to_ir_coloc <- bind_nonempty(gene_to_ir_coloc_rows)
if (nrow(gene_to_ir_coloc) > 0L) {
  gene_to_ir_coloc[, evidence_tier := dplyr::case_when(
    p12 != config$coloc$primary_p12 ~ "prior_sensitivity",
    PP.H4 >= config$coloc$tier_a_pp4 ~ "Tier_A",
    PP.H4 >= config$coloc$tier_b_pp4 & PP.H4 > PP.H3 ~ "Tier_B",
    PP.H4 >= config$coloc$tier_c_pp4 & PP.H4 > PP.H3 ~ "Tier_C",
    TRUE ~ "None"
  )]
}

map_position_instruments_to_rsid <- function(ir, reference) {
  if (nrow(ir) == 0L) return(empty_dt())
  reference <- unique(reference[, .(
    chr = as.character(chr), pos = as.integer(pos), SNP,
    ref_effect = effect_allele, ref_other = other_allele
  )], by = c("chr", "pos", "SNP"))
  ir_without_id <- copy(ir)
  ir_without_id[, SNP := NULL]
  mapped <- merge(
    ir_without_id, reference,
    by = c("chr", "pos"), all = FALSE, allow.cartesian = TRUE
  )
  same <- mapped$effect_allele == mapped$ref_effect &
    mapped$other_allele == mapped$ref_other
  flipped <- mapped$effect_allele == mapped$ref_other &
    mapped$other_allele == mapped$ref_effect
  strand_same <- complement_allele(mapped$effect_allele) ==
    mapped$ref_effect &
    complement_allele(mapped$other_allele) == mapped$ref_other
  strand_flipped <- complement_allele(mapped$effect_allele) ==
    mapped$ref_other &
    complement_allele(mapped$other_allele) == mapped$ref_effect
  keep <- same | flipped | strand_same | strand_flipped
  complement_rows <- (strand_same | strand_flipped)[keep]
  mapped <- mapped[keep]
  mapped[complement_rows, `:=`(
    effect_allele = complement_allele(effect_allele),
    other_allele = complement_allele(other_allele)
  )]
  mapped[, c("ref_effect", "ref_other") := NULL]
  unique(mapped, by = "SNP")
}

format_magic_exposure <- function(ir, trait) {
  if (nrow(ir) == 0L) return(NULL)
  ir <- copy(ir)
  ir[, F_stat := (beta / se)^2]
  ir <- ir[pval < magic_config$instrument_p &
             F_stat >= magic_config$minimum_F & grepl("^rs", SNP)]
  if (nrow(ir) == 0L) return(NULL)
  out <- tryCatch(
    TwoSampleMR::format_data(
      as.data.frame(ir), type = "exposure",
      snp_col = "SNP", beta_col = "beta", se_col = "se",
      eaf_col = "eaf", effect_allele_col = "effect_allele",
      other_allele_col = "other_allele", pval_col = "pval",
      chr_col = "chr", pos_col = "pos", samplesize_col = "N"
    ),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out$exposure <- trait
  out$id.exposure <- trait
  out
}

clump_multichromosome_exposure <- function(exposure) {
  if (is.null(exposure) || nrow(exposure) == 0L) return(NULL)
  pieces <- lapply(split(exposure, exposure$chr.exposure), clump_exposure)
  pieces <- Filter(function(x) !is.null(x) && nrow(x) > 0L, pieces)
  if (length(pieces) == 0L) return(NULL)
  out <- do.call(rbind, pieces)
  out[!duplicated(out$SNP), , drop = FALSE]
}

magic_instruments <- lapply(names(magic_raw), function(trait) {
  instruments <- copy(magic_raw[[trait]]$instruments)
  if (trait == "ISI_adjBMI") {
    instruments <- map_position_instruments_to_rsid(
      instruments, gwas_list[["ALM"]]
    )
  }
  if (nrow(instruments) > 0L) {
    reference_eaf <- gwas_list[["ALM"]]$eaf[
      match(instruments$SNP, gwas_list[["ALM"]]$SNP)
    ]
    missing_eaf <- !is.finite(instruments$eaf) |
      instruments$eaf <= 0 | instruments$eaf >= 1
    instruments[missing_eaf, eaf := reference_eaf[missing_eaf]]
    instruments[, maf := pmin(eaf, 1 - eaf)]
  }
  instruments
})
names(magic_instruments) <- names(magic_raw)

run_magic_to_sarcopenia_mr <- function(instruments, outcomes) {
  rows <- list()
  harmonised <- list()

  for (ir_trait in names(instruments)) {
    exposure <- format_magic_exposure(instruments[[ir_trait]], ir_trait)
    exposure <- clump_multichromosome_exposure(exposure)
    if (is.null(exposure) || nrow(exposure) == 0L) next

    for (sarcopenia_trait in names(outcomes)) {
      outcome <- format_outcome(
        outcomes[[sarcopenia_trait]], exposure$SNP, sarcopenia_trait
      )
      if (is.null(outcome)) next
      missing_eaf <- !is.finite(exposure$eaf.exposure)
      exposure$eaf.exposure[missing_eaf] <- outcome$eaf.outcome[
        match(exposure$SNP[missing_eaf], outcome$SNP)
      ]
      dat <- tryCatch(
        TwoSampleMR::harmonise_data(exposure, outcome, action = 3),
        error = function(e) NULL
      )
      if (is.null(dat)) next
      dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
      if (nrow(dat) == 0L) next

      methods <- if (nrow(dat) == 1L) {
        "mr_wald_ratio"
      } else {
        c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
      }
      estimates <- tryCatch(
        TwoSampleMR::mr(dat, method_list = methods),
        error = function(e) NULL
      )
      if (is.null(estimates) || nrow(estimates) == 0L) next
      estimates <- estimates[
        estimates$method %in% c("Wald ratio", "Inverse variance weighted"),
        , drop = FALSE
      ]
      if (nrow(estimates) == 0L) next

      key <- paste(ir_trait, sarcopenia_trait, sep = "__")
      harmonised[[key]] <- dat
      rows[[key]] <- data.table(
        ir_trait = ir_trait,
        sarcopenia_trait = sarcopenia_trait,
        method = estimates$method[1],
        nsnp = estimates$nsnp[1],
        beta = estimates$b[1],
        se = estimates$se[1],
        pval = estimates$pval[1],
        min_F = min((dat$beta.exposure / dat$se.exposure)^2),
        mean_F = mean((dat$beta.exposure / dat$se.exposure)^2)
      )
    }
  }
  list(mr = bind_nonempty(rows), harmonised = harmonised)
}

log_step("Running MAGIC IR traits to sarcopenia-related outcomes MR")
ir_sarcopenia <- run_magic_to_sarcopenia_mr(magic_instruments, gwas_list)
ir_to_sarcopenia_mr <- ir_sarcopenia$mr
if (nrow(ir_to_sarcopenia_mr) > 0L) {
  ir_to_sarcopenia_mr[, fdr_within_sarcopenia_trait := p.adjust(pval, "BH"),
                      by = sarcopenia_trait]
  ir_to_sarcopenia_mr[, fdr_global := p.adjust(pval, "BH")]
  ir_to_sarcopenia_mr[, `:=`(
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    mr_fdr_significant = fdr_within_sarcopenia_trait <
      config$multiple_testing$alpha,
    interpretation_scale = fifelse(
      ir_trait == "ISI_adjBMI",
      "effect per genetically higher insulin sensitivity",
      "effect per genetically higher insulin resistance or insulin level"
    )
  )]
}

two_step_chain <- empty_dt()
if (nrow(gene_to_ir_mr) > 0L && nrow(ir_to_sarcopenia_mr) > 0L) {
  step1 <- gene_to_ir_mr[, .(
    gene_ensembl, gene_symbol, ir_trait,
    beta_gene_ir = beta, se_gene_ir = se,
    p_gene_ir = pval, fdr_gene_ir = fdr_within_ir_trait,
    gene_ir_mr_supported = mr_fdr_significant
  )]
  step2 <- ir_to_sarcopenia_mr[, .(
    ir_trait, sarcopenia_trait,
    beta_ir_sarcopenia = beta, se_ir_sarcopenia = se,
    p_ir_sarcopenia = pval,
    fdr_ir_sarcopenia = fdr_within_sarcopenia_trait,
    ir_sarcopenia_mr_supported = mr_fdr_significant
  )]
  two_step_chain <- merge(step1, step2, by = "ir_trait",
                          allow.cartesian = TRUE)

  direct <- mr_main[gene_symbol %in% magic_targets$gene_symbol, .(
    gene_symbol, sarcopenia_trait = trait,
    beta_direct = beta, se_direct = se, p_direct = pval,
    fdr_direct = fdr_within_trait,
    direct_mr_supported = mr_associated
  )]
  two_step_chain <- merge(
    two_step_chain, direct,
    by = c("gene_symbol", "sarcopenia_trait"), all.x = TRUE
  )

  primary_magic_coloc <- if (nrow(gene_to_ir_coloc) > 0L) {
    gene_to_ir_coloc[
      p12 == config$coloc$primary_p12,
      .(gene_symbol, ir_trait, gene_ir_PP4 = PP.H4,
        gene_ir_coloc_tier = evidence_tier)
    ]
  } else {
    data.table(
      gene_symbol = character(), ir_trait = character(),
      gene_ir_PP4 = numeric(), gene_ir_coloc_tier = character()
    )
  }
  two_step_chain <- merge(
    two_step_chain, primary_magic_coloc,
    by = c("gene_symbol", "ir_trait"), all.x = TRUE
  )

  two_step_chain[, `:=`(
    indirect_beta = beta_gene_ir * beta_ir_sarcopenia,
    indirect_se = sqrt(
      beta_ir_sarcopenia^2 * se_gene_ir^2 +
        beta_gene_ir^2 * se_ir_sarcopenia^2
    )
  )]
  two_step_chain[, indirect_pval := normal_pvalue(indirect_beta, indirect_se)]
  two_step_chain[, indirect_fdr := p.adjust(indirect_pval, "BH")]
  two_step_chain[, `:=`(
    indirect_direction_matches_direct =
      sign(indirect_beta) == sign(beta_direct),
    complete_chain_supported =
      gene_ir_mr_supported %in% TRUE &
      gene_ir_PP4 >= config$coloc$tier_a_pp4 &
      ir_sarcopenia_mr_supported %in% TRUE &
      direct_mr_supported %in% TRUE,
    interpretation =
      "Exploratory product-of-coefficients estimate; not proof of mediation"
  )]
}

magic_input_qc <- bind_nonempty(lapply(names(magic_raw), function(trait) {
  raw <- magic_raw[[trait]]
  mapped <- magic_regional[[trait]]
  instruments <- magic_instruments[[trait]]
  data.table(
    trait = trait,
    source_file = basename(magic_config$files[[trait]]),
    stated_build = magic_config$build,
    total_variants_scanned = raw$total_rows,
    regional_variants_read = nrow(raw$regional),
    regional_variants_mapped_to_gtex = nrow(mapped),
    regional_missing_eaf = if (nrow(mapped) > 0L)
      sum(!is.finite(mapped$eaf)) else 0L,
    genomewide_strong_variants_before_clumping = nrow(raw$instruments),
    strong_variants_with_rsid = if (nrow(instruments) > 0L)
      sum(grepl("^rs", instruments$SNP)) else 0L
  )
}))
magic_input_qc <- merge(
  magic_input_qc, position_qc, by = "trait", all.x = TRUE
)

magic_readme <- data.table(
  item = c(
    "analysis_scope", "target_selection", "gene_to_ir_testing",
    "ir_to_sarcopenia_testing", "coordinate_harmonisation",
    "HOMA_IR_frequency", "ISI_direction", "other_IR_direction",
    "two_step_interpretation", "independent_replication"
  ),
  description = c(
    "Appended mechanistic extension; original pipeline and results are unchanged.",
    "Only candidates already classified as Tier A in the primary colocalisation analysis.",
    "Skeletal-muscle cis-eQTL MR plus coloc.abf with prior sensitivity; BH FDR is reported.",
    "Genome-wide significant MAGIC instruments, local EUR LD clumping, and MR against the four existing outcomes.",
    "GTEx v8 GRCh38 eQTL coordinates were lifted to MAGIC GRCh37 with the local hg38ToHg19 chain.",
    "HOMA-IR VCF lacks allele frequency; matched GTEx skeletal-muscle eQTL frequency is used only where required for harmonisation/colocalisation.",
    "Higher ISI indicates greater insulin sensitivity, so its beta direction is opposite to an insulin-resistance scale.",
    "Higher HOMA-IR or fasting insulin is interpreted as greater insulin resistance or insulin level.",
    "Indirect effects use the delta-method product of coefficients and are exploratory, not proof of mediation.",
    "FinnGen and Zenodo outcomes were not used because the corresponding files were not supplied; existing LOW_GRIP remains the older-adult multi-cohort outcome."
  )
)

magic_output_tables <- list(
  `01_input_QC` = as.data.frame(magic_input_qc),
  `02_target_eQTL_liftover_QC` = as.data.frame(
    magic_eqtl_map[, .(
      gene_ensembl, gene_symbol, SNP, chr38, pos38, chr37, pos37
    )]
  ),
  `03_gene_to_IR_MR` = as.data.frame(
    if (nrow(gene_to_ir_mr) > 0L) gene_to_ir_mr else
      data.table(status = "no_analysable_results")
  ),
  `04_gene_to_IR_coloc_all` = as.data.frame(
    if (nrow(gene_to_ir_coloc) > 0L) gene_to_ir_coloc else
      data.table(status = "no_analysable_results")
  ),
  `05_gene_to_IR_coloc_primary` = as.data.frame(
    if (nrow(gene_to_ir_coloc) > 0L)
      gene_to_ir_coloc[p12 == config$coloc$primary_p12] else
        data.table(status = "no_analysable_results")
  ),
  `06_IR_to_sarcopenia_MR` = as.data.frame(
    if (nrow(ir_to_sarcopenia_mr) > 0L) ir_to_sarcopenia_mr else
      data.table(status = "no_analysable_results")
  ),
  `07_two_step_evidence_chain` = as.data.frame(
    if (nrow(two_step_chain) > 0L) two_step_chain else
      data.table(status = "no_analysable_results")
  ),
  `08_README` = as.data.frame(magic_readme)
)

openxlsx::write.xlsx(
  magic_output_tables,
  file = magic_config$output_workbook,
  overwrite = TRUE,
  asTable = TRUE
)
magic_workbook_repair <- repair_missing_ooxml_relationships(
  magic_config$output_workbook
)
magic_actual_sheets <- openxlsx::getSheetNames(magic_config$output_workbook)
if (!identical(magic_actual_sheets, names(magic_output_tables))) {
  stop("MAGIC extension workbook validation failed: worksheet names differ.")
}

log_step(
  "MAGIC extension complete: ", magic_config$output_workbook,
  "; Tier A genes = ", nrow(magic_targets),
  "; gene-to-IR MR tests = ", nrow(gene_to_ir_mr),
  "; IR-to-sarcopenia MR tests = ", nrow(ir_to_sarcopenia_mr),
  "; complete two-step chains = ",
  if (nrow(two_step_chain) > 0L)
    sum(two_step_chain$complete_chain_supported, na.rm = TRUE) else 0L
)


## 17. Human skeletal-muscle transcriptomic validation -----------------------
##
## This appended section validates the prespecified Tier A genes in two
## independent public human skeletal-muscle transcriptomic datasets. Existing
## analyses and outputs above are unchanged. The primary contrast in both
## datasets is sarcopenia versus age-matched older controls, with age included
## as a covariate. GSE167186 young and unclassified samples are retained in the
## metadata/QC record but excluded from the primary disease contrast.

log_step("Starting appended GEO human skeletal-muscle validation")

geo_required_packages <- c("DESeq2", "edgeR")
geo_missing_packages <- geo_required_packages[
  !vapply(geo_required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(geo_missing_packages) > 0L) {
  stop(
    "Please install the following packages for the GEO extension: ",
    paste(geo_missing_packages, collapse = ", ")
  )
}

geo_config <- list(
  files = list(
    GSE111016_counts = file.path(
      config$project_dir,
      "GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv.gz"
    ),
    GSE111016_metadata = file.path(
      config$project_dir, "GSE111016_series_matrix.txt.gz"
    ),
    GSE167186_counts = file.path(
      config$project_dir, "GSE167186_counts.csv.gz"
    ),
    GSE167186_metadata = file.path(
      config$project_dir, "GSE167186-GPL20301_series_matrix.txt.gz"
    )
  ),
  minimum_count = 10L,
  candidate_fdr = 0.05,
  nominal_p = 0.05,
  output_workbook = result_path("GEO_human_muscle_validation.xlsx"),
  output_forest_pdf = result_path("GEO_TierA_candidate_forest.pdf"),
  output_forest_png = result_path("GEO_TierA_candidate_forest.png")
)

assert_files_exist(unlist(geo_config$files), "GEO transcriptomic validation")

if (!exists("magic_targets", inherits = FALSE) || nrow(magic_targets) == 0L) {
  stop("Tier A candidate genes are unavailable for the GEO extension.")
}

geo_targets <- unique(
  magic_targets[, .(
    gene_ensembl = strip_ensembl_version(gene_ensembl),
    gene_symbol = as.character(gene_symbol),
    gene_layer = as.character(gene_layer)
  )]
)
data.table::setorder(geo_targets, gene_symbol)

read_geo_series_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  parse_sample_line <- function(prefix) {
    line <- lines[startsWith(lines, prefix)]
    if (length(line) == 0L) return(character())
    values <- strsplit(line[1], "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", values)
  }

  sample_title <- parse_sample_line("!Sample_title")
  geo_accession <- parse_sample_line("!Sample_geo_accession")
  if (length(sample_title) == 0L ||
      length(sample_title) != length(geo_accession)) {
    stop("Unable to recover aligned GEO sample titles/accessions from: ", path)
  }

  characteristic_lines <- lines[
    startsWith(lines, "!Sample_characteristics_ch1")
  ]
  long_rows <- lapply(seq_along(characteristic_lines), function(i) {
    values <- strsplit(characteristic_lines[i], "\t", fixed = TRUE)[[1]][-1]
    values <- gsub('^"|"$', "", values)
    length(values) <- length(sample_title)
    data.table(
      sample_index = seq_along(sample_title),
      characteristic_row = i,
      raw_value = trimws(values)
    )
  })
  characteristic_long <- bind_nonempty(long_rows)
  characteristic_long <- characteristic_long[
    !is.na(raw_value) & nzchar(raw_value) & grepl(":", raw_value, fixed = TRUE)
  ]
  characteristic_long[, `:=`(
    key = tolower(trimws(sub(":.*$", "", raw_value))),
    value = trimws(sub("^[^:]+:[[:space:]]*", "", raw_value))
  )]

  extract_characteristic <- function(pattern) {
    selected <- characteristic_long[grepl(pattern, key, perl = TRUE)]
    if (nrow(selected) == 0L) return(rep(NA_character_, length(sample_title)))
    selected <- selected[order(characteristic_row)]
    first_value <- selected[, .(value = value[1]), by = sample_index]
    answer <- rep(NA_character_, length(sample_title))
    answer[first_value$sample_index] <- first_value$value
    answer
  }

  metadata <- data.table(
    sample_index = seq_along(sample_title),
    sample_title = sample_title,
    geo_accession = geo_accession,
    sarcopenia_status = extract_characteristic("^sarcopenia status$"),
    group = extract_characteristic("^group$"),
    sex = extract_characteristic("^sex$"),
    age_text = extract_characteristic("^age( \\(yr\\))?$")
  )
  metadata[, age := suppressWarnings(as.numeric(age_text))]
  metadata
}

prepare_count_matrix <- function(count_table, gene_column, sample_columns) {
  gene_ids <- as.character(count_table[[gene_column]])
  count_matrix <- as.matrix(count_table[, ..sample_columns])
  storage.mode(count_matrix) <- "numeric"
  if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
    stop("Count matrix contains missing, non-finite, or negative values.")
  }
  count_matrix <- round(count_matrix)
  rownames(count_matrix) <- gene_ids
  if (anyDuplicated(gene_ids)) {
    count_matrix <- rowsum(count_matrix, group = gene_ids, reorder = FALSE)
  }
  storage.mode(count_matrix) <- "integer"
  count_matrix
}

run_geo_deseq <- function(
    count_matrix, metadata, dataset, gene_symbols,
    control_label = "Control", case_label = "Sarcopenia") {
  if (!identical(colnames(count_matrix), metadata$count_column)) {
    stop(dataset, ": count columns and metadata are not identically aligned.")
  }
  if (any(!is.finite(metadata$age))) {
    stop(dataset, ": non-finite age in primary analysis samples.")
  }

  metadata[, analysis_group := factor(
    analysis_group, levels = c(control_label, case_label)
  )]
  if (!identical(levels(metadata$analysis_group), c(control_label, case_label)) ||
      any(table(metadata$analysis_group) == 0L)) {
    stop(dataset, ": both case and control groups are required.")
  }
  metadata[, age_z := as.numeric(scale(age))]

  keep <- edgeR::filterByExpr(
    count_matrix, group = metadata$analysis_group,
    min.count = geo_config$minimum_count
  )
  if (sum(keep) < 100L) {
    stop(dataset, ": too few genes passed expression filtering.")
  }
  filtered_counts <- count_matrix[keep, , drop = FALSE]

  col_data <- as.data.frame(metadata)
  rownames(col_data) <- col_data$count_column
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = col_data,
    design = ~ age_z + analysis_group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  fit <- DESeq2::results(
    dds,
    contrast = c("analysis_group", case_label, control_label),
    independentFiltering = TRUE,
    alpha = geo_config$candidate_fdr
  )
  fit <- as.data.table(as.data.frame(fit), keep.rownames = "gene_id")
  data.table::setnames(
    fit,
    c("log2FoldChange", "pvalue", "padj"),
    c("log2FC_sarcopenia_vs_control", "p_value", "whole_transcriptome_FDR")
  )
  fit[, `:=`(
    dataset = dataset,
    gene_symbol = unname(gene_symbols[gene_id]),
    contrast = paste0(case_label, "_vs_", control_label),
    n_control = sum(metadata$analysis_group == control_label),
    n_sarcopenia = sum(metadata$analysis_group == case_label)
  )]
  data.table::setcolorder(
    fit,
    c(
      "dataset", "gene_id", "gene_symbol", "contrast",
      "n_control", "n_sarcopenia", "baseMean",
      "log2FC_sarcopenia_vs_control", "lfcSE", "stat", "p_value",
      "whole_transcriptome_FDR"
    )
  )
  data.table::setorder(fit, p_value)

  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  log2_normalized <- log2(normalized_counts + 1)
  variable_genes <- head(
    order(apply(log2_normalized, 1L, stats::var), decreasing = TRUE),
    min(500L, nrow(log2_normalized))
  )
  pca <- stats::prcomp(
    t(log2_normalized[variable_genes, , drop = FALSE]),
    center = TRUE, scale. = FALSE
  )
  pca_table <- cbind(
    metadata[, .(
      dataset, sample_title, geo_accession, count_column,
      analysis_group = as.character(analysis_group), age,
      library_size
    )],
    data.table(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2]
    )
  )

  list(
    full_results = fit,
    normalized_counts = normalized_counts,
    metadata = metadata,
    pca = pca_table,
    n_genes_input = nrow(count_matrix),
    n_genes_tested = nrow(filtered_counts)
  )
}

log_step("Reading and analysing GSE111016")
gse111_metadata <- read_geo_series_metadata(
  geo_config$files$GSE111016_metadata
)
gse111_counts_dt <- data.table::fread(geo_config$files$GSE111016_counts)
data.table::setnames(gse111_counts_dt, 1L, "gene_id")
gse111_counts_dt[, gene_id := strip_ensembl_version(gene_id)]
gse111_sample_columns <- setdiff(names(gse111_counts_dt), "gene_id")
gse111_metadata[, count_column := sub("[[:space:]]*\\[sss\\]$", "", sample_title)]
gse111_metadata <- gse111_metadata[match(gse111_sample_columns, count_column)]
if (anyNA(gse111_metadata$count_column) ||
    !identical(gse111_metadata$count_column, gse111_sample_columns)) {
  stop("GSE111016 sample titles do not match count-matrix columns.")
}
gse111_metadata[, analysis_group := fifelse(
  tolower(sarcopenia_status) == "yes", "Sarcopenia",
  fifelse(tolower(sarcopenia_status) == "no", "Control", NA_character_)
)]
if (anyNA(gse111_metadata$analysis_group)) {
  stop("GSE111016 sarcopenia status could not be assigned for every sample.")
}
gse111_count_matrix <- prepare_count_matrix(
  gse111_counts_dt, "gene_id", gse111_sample_columns
)
gse111_metadata[, `:=`(
  dataset = "GSE111016",
  library_size = colSums(gse111_count_matrix)
)]
gse111_gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(gse111_count_matrix),
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"
)
gse111_fit <- run_geo_deseq(
  gse111_count_matrix, gse111_metadata, "GSE111016",
  gse111_gene_symbols
)

log_step("Reading and analysing GSE167186")
gse167_metadata_all <- read_geo_series_metadata(
  geo_config$files$GSE167186_metadata
)
gse167_counts_dt <- data.table::fread(geo_config$files$GSE167186_counts)
data.table::setnames(gse167_counts_dt, 1L, "gene_symbol")
gse167_sample_columns <- setdiff(names(gse167_counts_dt), "gene_symbol")
gse167_metadata_all[, count_column := sample_title]
gse167_metadata_all <- gse167_metadata_all[match(
  gse167_sample_columns, count_column
)]
if (anyNA(gse167_metadata_all$count_column) ||
    !identical(gse167_metadata_all$count_column, gse167_sample_columns)) {
  stop("GSE167186 sample titles do not match count-matrix columns.")
}
gse167_metadata_all[, analysis_group := fifelse(
  tolower(group) == "sarcopenia", "Sarcopenia",
  fifelse(tolower(group) == "old healthy", "Control", NA_character_)
)]
gse167_count_matrix_all <- prepare_count_matrix(
  gse167_counts_dt, "gene_symbol", gse167_sample_columns
)
gse167_metadata_all[, `:=`(
  dataset = "GSE167186",
  library_size = colSums(gse167_count_matrix_all),
  included_primary_contrast = !is.na(analysis_group)
)]
gse167_primary_metadata <- gse167_metadata_all[included_primary_contrast == TRUE]
gse167_primary_columns <- gse167_primary_metadata$count_column
gse167_count_matrix <- gse167_count_matrix_all[
  , gse167_primary_columns, drop = FALSE
]
gse167_gene_symbols <- setNames(
  rownames(gse167_count_matrix), rownames(gse167_count_matrix)
)
gse167_fit <- run_geo_deseq(
  gse167_count_matrix, gse167_primary_metadata, "GSE167186",
  gse167_gene_symbols
)

build_candidate_results <- function(fit_object, targets, id_type) {
  fit <- copy(fit_object$full_results)
  if (id_type == "ensembl") {
    fit[, gene_symbol := NULL]
    candidate <- merge(
      targets,
      fit,
      by.x = "gene_ensembl", by.y = "gene_id",
      all.x = TRUE, sort = FALSE
    )
  } else {
    candidate <- merge(
      targets,
      fit,
      by = "gene_symbol",
      all.x = TRUE, sort = FALSE
    )
  }
  candidate[, candidate_set_FDR := p.adjust(
    p_value, method = "BH", n = nrow(targets)
  )]
  candidate[, `:=`(
    detected_after_filtering = is.finite(baseMean),
    nominally_associated = is.finite(p_value) & p_value < geo_config$nominal_p,
    candidate_FDR_significant =
      is.finite(candidate_set_FDR) & candidate_set_FDR < geo_config$candidate_fdr
  )]
  data.table::setorder(candidate, gene_symbol)
  candidate
}

gse111_candidates <- build_candidate_results(
  gse111_fit, geo_targets, "ensembl"
)
gse167_candidates <- build_candidate_results(
  gse167_fit, geo_targets, "symbol"
)

candidate_meta <- merge(
  gse111_candidates[, .(
    gene_symbol, gene_ensembl, gene_layer,
    beta_GSE111016 = log2FC_sarcopenia_vs_control,
    se_GSE111016 = lfcSE,
    p_GSE111016 = p_value,
    fdr6_GSE111016 = candidate_set_FDR
  )],
  gse167_candidates[, .(
    gene_symbol,
    beta_GSE167186 = log2FC_sarcopenia_vs_control,
    se_GSE167186 = lfcSE,
    p_GSE167186 = p_value,
    fdr6_GSE167186 = candidate_set_FDR
  )],
  by = "gene_symbol", all = TRUE
)
candidate_meta[, `:=`(
  direction_concordant =
    is.finite(beta_GSE111016) & is.finite(beta_GSE167186) &
    sign(beta_GSE111016) == sign(beta_GSE167186),
  meta_beta = {
    w1 <- 1 / se_GSE111016^2
    w2 <- 1 / se_GSE167186^2
    (w1 * beta_GSE111016 + w2 * beta_GSE167186) / (w1 + w2)
  },
  meta_se = sqrt(1 / (1 / se_GSE111016^2 + 1 / se_GSE167186^2))
)]
candidate_meta[
  !is.finite(se_GSE111016) | se_GSE111016 <= 0 |
    !is.finite(se_GSE167186) | se_GSE167186 <= 0,
  `:=`(meta_beta = NA_real_, meta_se = NA_real_)
]
candidate_meta[, `:=`(
  meta_z = meta_beta / meta_se,
  meta_p = 2 * stats::pnorm(-abs(meta_beta / meta_se))
)]
candidate_meta[, meta_FDR6 := p.adjust(
  meta_p, method = "BH", n = nrow(geo_targets)
)]
candidate_meta[, support_tier := dplyr::case_when(
  direction_concordant %in% TRUE &
    p_GSE111016 < geo_config$nominal_p &
    p_GSE167186 < geo_config$nominal_p ~ "Replicated_in_both_cohorts",
  direction_concordant %in% TRUE & meta_FDR6 < geo_config$candidate_fdr &
    (p_GSE111016 < geo_config$nominal_p |
       p_GSE167186 < geo_config$nominal_p) ~ "Meta_supported",
  direction_concordant %in% TRUE &
    (p_GSE111016 < geo_config$nominal_p |
       p_GSE167186 < geo_config$nominal_p) ~ "Single_cohort_signal",
  TRUE ~ "No_cross_cohort_support"
)]
data.table::setorder(candidate_meta, meta_p)

geo_discovery_context <- merge(
  geo_targets,
  mr_main[gene_symbol %in% geo_targets$gene_symbol, .(
    gene_symbol, trait, discovery_MR_beta = beta,
    discovery_MR_se = se, discovery_MR_p = pval,
    discovery_MR_FDR = fdr_within_trait
  )],
  by = "gene_symbol", all.x = TRUE
)
geo_discovery_context <- merge(
  geo_discovery_context,
  coloc_primary[
    evidence_tier == "Tier_A" & gene_symbol %in% geo_targets$gene_symbol,
    .(gene_symbol, trait, discovery_coloc_PP4 = PP.H4)
  ],
  by = c("gene_symbol", "trait"), all.x = TRUE
)
if (exists("fusion_replication", inherits = FALSE) &&
    nrow(fusion_replication) > 0L) {
  geo_discovery_context <- merge(
    geo_discovery_context,
    fusion_replication[, .(
      gene_symbol, trait, fusion_MR_beta = fusion_mr_beta,
      fusion_MR_FDR = fusion_mr_fdr,
      fusion_coloc_PP4 = fusion_PP.H4,
      fusion_replicated = replicated
    )],
    by = c("gene_symbol", "trait"), all.x = TRUE
  )
}

geo_sample_qc <- rbindlist(
  list(
    gse111_fit$pca[, .(
      dataset, sample_title, geo_accession, count_column,
      original_group = gse111_metadata$analysis_group,
      primary_analysis_group = analysis_group,
      age, library_size, included_primary_contrast = TRUE, PC1, PC2
    )],
    merge(
      gse167_metadata_all[, .(
        dataset, sample_title, geo_accession, count_column,
        original_group = group, age, library_size,
        included_primary_contrast
      )],
      gse167_fit$pca[, .(
        count_column, primary_analysis_group = analysis_group, PC1, PC2
      )],
      by = "count_column", all.x = TRUE, sort = FALSE
    )
  ),
  use.names = TRUE, fill = TRUE
)

geo_dataset_qc <- data.table(
  dataset = c("GSE111016", "GSE167186"),
  n_samples_downloaded = c(nrow(gse111_metadata), nrow(gse167_metadata_all)),
  n_samples_primary_contrast = c(
    nrow(gse111_fit$metadata), nrow(gse167_fit$metadata)
  ),
  n_control = c(
    sum(gse111_fit$metadata$analysis_group == "Control"),
    sum(gse167_fit$metadata$analysis_group == "Control")
  ),
  n_sarcopenia = c(
    sum(gse111_fit$metadata$analysis_group == "Sarcopenia"),
    sum(gse167_fit$metadata$analysis_group == "Sarcopenia")
  ),
  n_genes_input = c(gse111_fit$n_genes_input, gse167_fit$n_genes_input),
  n_genes_tested = c(gse111_fit$n_genes_tested, gse167_fit$n_genes_tested),
  design = "~ age_z + analysis_group",
  contrast = "Sarcopenia_vs_older_control"
)

forest_rows <- rbindlist(list(
  gse111_candidates[, .(
    gene_symbol, dataset = "GSE111016",
    beta = log2FC_sarcopenia_vs_control, se = lfcSE,
    p_value, candidate_FDR = candidate_set_FDR
  )],
  gse167_candidates[, .(
    gene_symbol, dataset = "GSE167186",
    beta = log2FC_sarcopenia_vs_control, se = lfcSE,
    p_value, candidate_FDR = candidate_set_FDR
  )],
  candidate_meta[, .(
    gene_symbol, dataset = "Fixed-effect meta-analysis",
    beta = meta_beta, se = meta_se,
    p_value = meta_p, candidate_FDR = meta_FDR6
  )]
), use.names = TRUE, fill = TRUE)
forest_rows[, `:=`(
  lower_95CI = beta - 1.96 * se,
  upper_95CI = beta + 1.96 * se,
  dataset = factor(
    dataset,
    levels = c("GSE111016", "GSE167186", "Fixed-effect meta-analysis")
  ),
  gene_symbol = factor(gene_symbol, levels = rev(sort(unique(gene_symbol))))
)]

geo_forest <- ggplot(
  forest_rows[is.finite(beta) & is.finite(se)],
  aes(x = beta, y = gene_symbol, colour = dataset, shape = dataset)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(
    aes(xmin = lower_95CI, xmax = upper_95CI),
    width = 0.18, orientation = "y",
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(
    "GSE111016" = "#0072B2",
    "GSE167186" = "#D55E00",
    "Fixed-effect meta-analysis" = "#009E73"
  )) +
  labs(
    x = "Log2 fold change: sarcopenia versus older control (95% CI)",
    y = NULL,
    colour = NULL,
    shape = NULL,
    title = "Human skeletal-muscle expression of Tier A candidates",
    subtitle = "Age-adjusted DESeq2 estimates; meta-analysis is exploratory"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  geo_config$output_forest_pdf, geo_forest,
  width = 8.2, height = 5.2, units = "in"
)
ggsave(
  geo_config$output_forest_png, geo_forest,
  width = 8.2, height = 5.2, units = "in", dpi = 320
)

geo_readme <- data.table(
  item = c(
    "scope", "prespecified_targets", "primary_contrast",
    "GSE111016", "GSE167186", "model", "multiple_testing",
    "meta_analysis", "support_tier", "interpretation"
  ),
  description = c(
    "Appended public human skeletal-muscle validation; all preceding code and outputs are unchanged.",
    "Only the six candidates classified as Tier A in the primary colocalisation analysis are used for candidate-set inference.",
    "Sarcopenia versus older non-sarcopenic control; positive log2FC means higher expression in sarcopenic muscle.",
    "Forty male participants of Chinese descent (20 sarcopenia, 20 control); vastus lateralis; age adjusted.",
    "Primary analysis excludes Young Healthy and UNCLASSIFIED samples and compares Sarcopenia with Old Healthy; age adjusted.",
    "DESeq2 negative-binomial model after edgeR filterByExpr; design is standardized age plus disease group.",
    "Both whole-transcriptome BH FDR and prespecified six-gene BH FDR are reported. Nominal P values are not treated as definitive validation.",
    "Inverse-variance fixed-effect synthesis of the two DESeq2 log2 fold changes; exploratory because only two heterogeneous cohorts are available.",
    "Replicated requires same direction and nominal P<0.05 in both cohorts; Meta_supported requires concordant direction, meta FDR<0.05, and nominal support in at least one cohort.",
    "Expression concordance supports human disease-tissue relevance but does not by itself establish causality, drug efficacy, or mediation by insulin resistance."
  )
)

geo_output_tables <- list(
  `01_dataset_QC` = as.data.frame(geo_dataset_qc),
  `02_sample_QC` = as.data.frame(geo_sample_qc),
  `03_GSE111016_all_DE` = as.data.frame(gse111_fit$full_results),
  `04_GSE111016_TierA` = as.data.frame(gse111_candidates),
  `05_GSE167186_all_DE` = as.data.frame(gse167_fit$full_results),
  `06_GSE167186_TierA` = as.data.frame(gse167_candidates),
  `07_cross_cohort_TierA` = as.data.frame(candidate_meta),
  `08_discovery_context` = as.data.frame(geo_discovery_context),
  `09_forest_plot_data` = as.data.frame(forest_rows),
  `10_README` = as.data.frame(geo_readme)
)

openxlsx::write.xlsx(
  geo_output_tables,
  file = geo_config$output_workbook,
  overwrite = TRUE,
  asTable = TRUE
)
geo_workbook <- openxlsx::loadWorkbook(geo_config$output_workbook)
for (sheet_name in names(geo_output_tables)) {
  sheet_data <- geo_output_tables[[sheet_name]]
  openxlsx::freezePane(geo_workbook, sheet = sheet_name, firstRow = TRUE)
  openxlsx::showGridLines(geo_workbook, sheet = sheet_name, showGridLines = FALSE)
  if (ncol(sheet_data) > 0L) {
    column_widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- as.character(sheet_data[[column_index]])
      values <- values[!is.na(values)]
      values <- head(values, 500L)
      observed_width <- if (length(values) > 0L) max(nchar(values)) else 0L
      min(
        28,
        max(10, nchar(names(sheet_data)[column_index]) + 2L, observed_width + 2L)
      )
    }, numeric(1))
    if (identical(sheet_name, "02_cluster_annotation")) {
      widths <- c(18, 22, 32, 70)
      openxlsx::addStyle(
        sn_workbook,
        sheet = sheet_name,
        style = openxlsx::createStyle(wrapText = TRUE, valign = "top"),
        rows = 2:(nrow(sheet_data) + 1L),
        cols = 4L,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
    if (identical(sheet_name, "10_README")) {
      column_widths <- c(24, 80)
      readme_wrap_style <- openxlsx::createStyle(
        wrapText = TRUE, valign = "top"
      )
      openxlsx::addStyle(
        geo_workbook, sheet = sheet_name, style = readme_wrap_style,
        rows = 2:(nrow(sheet_data) + 1L), cols = 2L,
        gridExpand = TRUE, stack = TRUE
      )
    }
    openxlsx::setColWidths(
      geo_workbook, sheet = sheet_name,
      cols = seq_len(ncol(sheet_data)), widths = column_widths
    )
  }
}
openxlsx::saveWorkbook(
  geo_workbook, geo_config$output_workbook, overwrite = TRUE
)
geo_workbook_repair <- repair_missing_ooxml_relationships(
  geo_config$output_workbook
)
geo_actual_sheets <- openxlsx::getSheetNames(geo_config$output_workbook)
if (!identical(geo_actual_sheets, names(geo_output_tables))) {
  stop("GEO extension workbook validation failed: worksheet names differ.")
}

log_step(
  "GEO validation complete: ", geo_config$output_workbook,
  "; Tier A genes = ", nrow(geo_targets),
  "; replicated in both cohorts = ",
  sum(candidate_meta$support_tier == "Replicated_in_both_cohorts", na.rm = TRUE),
  "; meta-supported = ",
  sum(candidate_meta$support_tier == "Meta_supported", na.rm = TRUE)
)




# 18. GSE167186 single-nucleus cell-type localization and donor-level analysis------


sn_required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "edgeR",
  "ggplot2", "openxlsx", "patchwork"
)
sn_missing_packages <- sn_required_packages[
  !vapply(sn_required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(sn_missing_packages) > 0L) {
  stop(
    "Section 18 requires installed packages: ",
    paste(sn_missing_packages, collapse = ", ")
  )
}

sn_config <- list(
  raw_tar = "GSE167186_RAW.tar",
  metadata_xlsx = "GSE167186_SimplifiedMetadataSheet.xlsx",
  output_dir = "results_public_release",
  output_workbook = file.path(
    "results_public_release", "GSE167186_snRNA_celltype_validation.xlsx"
  ),
  output_umap_pdf = file.path(
    "results_public_release", "GSE167186_snRNA_celltype_UMAP.pdf"
  ),
  output_umap_png = file.path(
    "results_public_release", "GSE167186_snRNA_celltype_UMAP.png"
  ),
  output_dotplot_pdf = file.path(
    "results_public_release", "GSE167186_TierA_celltype_expression.pdf"
  ),
  output_dotplot_png = file.path(
    "results_public_release", "GSE167186_TierA_celltype_expression.png"
  ),
  min_features = 200L,
  variable_features = 2000L,
  dimensions = 1:10,
  cluster_resolution = 0.2,
  minimum_nuclei_per_pseudobulk = 20L,
  seed = 20260726L
)

if (!file.exists(sn_config$raw_tar)) {
  stop("Missing GSE167186 single-nucleus archive: ", sn_config$raw_tar)
}
if (!file.exists(sn_config$metadata_xlsx)) {
  stop("Missing GSE167186 metadata workbook: ", sn_config$metadata_xlsx)
}
dir.create(sn_config$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(sn_config$seed)

sn_metadata_raw <- openxlsx::read.xlsx(
  sn_config$metadata_xlsx,
  sheet = 1,
  colNames = FALSE,
  skipEmptyRows = FALSE,
  skipEmptyCols = FALSE
)
sn_header_row <- which(
  grepl("Sample name", as.character(sn_metadata_raw[[1L]]), fixed = TRUE)
)[1L]
if (is.na(sn_header_row)) {
  stop("Could not identify the sample-name row in GSE167186 metadata.")
}

sn_metadata <- data.table::data.table(
  sample_name = trimws(as.character(
    sn_metadata_raw[(sn_header_row + 1L):nrow(sn_metadata_raw), 1L]
  )),
  geo_accession = trimws(as.character(
    sn_metadata_raw[(sn_header_row + 1L):nrow(sn_metadata_raw), 2L]
  )),
  raw_file = trimws(as.character(
    sn_metadata_raw[(sn_header_row + 1L):nrow(sn_metadata_raw), 5L]
  ))
)
sn_metadata <- sn_metadata[
  !is.na(sample_name) & !is.na(geo_accession) & !is.na(raw_file)
]
sn_metadata[, sample_name := gsub("\\u00a0", "", sample_name)]
sn_metadata[, geo_accession := gsub("\\u00a0", "", geo_accession)]
sn_metadata[, raw_file := gsub("\\u00a0", "", raw_file)]
sn_metadata[, group := ifelse(
  grepl("^Old", sample_name), "Old",
  ifelse(grepl("^You", sample_name), "Young", NA_character_)
)]
sn_metadata[, donor_id := sub("^(Old|You)", "", sample_name)]
sn_metadata[, archive_file := paste0(geo_accession, "_", raw_file)]
sn_metadata[, csv_archive_file := sub(
  "_filtered_feature_bc_matrix\\.h5$", ".csv.gz", archive_file
)]

if (nrow(sn_metadata) != 17L || anyNA(sn_metadata$group)) {
  stop("Expected 17 classified single-nucleus donors (11 Old and 6 Young).")
}
if (!identical(
  as.integer(table(sn_metadata$group)[c("Old", "Young")]),
  c(11L, 6L)
)) {
  stop("GSE167186 donor groups do not match the expected 11 Old and 6 Young.")
}

sn_archive_members <- utils::untar(sn_config$raw_tar, list = TRUE)
sn_missing_members <- setdiff(sn_metadata$csv_archive_file, sn_archive_members)
if (length(sn_missing_members) > 0L) {
  stop(
    "The GSE167186 archive is missing processed count matrices: ",
    paste(sn_missing_members, collapse = ", ")
  )
}

sn_extract_dir <- tempfile("GSE167186_snRNA_")
dir.create(sn_extract_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(sn_extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(
  sn_config$raw_tar,
  files = sn_metadata$csv_archive_file,
  exdir = sn_extract_dir
)

read_gse167186_processed_counts <- function(file_path, sample_label) {
  counts_table <- data.table::fread(
    file_path,
    check.names = FALSE,
    showProgress = FALSE
  )
  gene_symbols <- make.unique(as.character(counts_table[[1L]]))
  counts_table[[1L]] <- NULL
  counts_dense <- as.matrix(counts_table)
  rm(counts_table)
  counts_sparse <- methods::as(counts_dense, "dgCMatrix")
  rm(counts_dense)
  rownames(counts_sparse) <- gene_symbols
  colnames(counts_sparse) <- paste0(sample_label, "_", colnames(counts_sparse))
  nuclei_input <- ncol(counts_sparse)

  object <- Seurat::CreateSeuratObject(
    counts = counts_sparse,
    project = "GSE167186",
    min.features = sn_config$min_features
  )
  object$sample_name <- sample_label
  object$group <- sn_metadata[sample_name == sample_label, group][1L]
  object$geo_accession <- sn_metadata[
    sample_name == sample_label, geo_accession
  ][1L]
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object, pattern = "^MT-"
  )
  attr(object, "nuclei_input") <- nuclei_input
  object
}

sn_objects <- vector("list", nrow(sn_metadata))
names(sn_objects) <- sn_metadata$sample_name
sn_input_qc <- vector("list", nrow(sn_metadata))

for (sample_index in seq_len(nrow(sn_metadata))) {
  sample_row <- sn_metadata[sample_index]
  sample_file <- file.path(sn_extract_dir, sample_row$csv_archive_file)
  sample_object <- read_gse167186_processed_counts(
    sample_file, sample_row$sample_name
  )
  sn_objects[[sample_row$sample_name]] <- sample_object
  sn_input_qc[[sample_index]] <- data.table::data.table(
    sample_name = sample_row$sample_name,
    geo_accession = sample_row$geo_accession,
    group = sample_row$group,
    nuclei_input = attr(sample_object, "nuclei_input"),
    nuclei_retained = ncol(sample_object),
    median_nFeature_RNA = stats::median(sample_object$nFeature_RNA),
    median_nCount_RNA = stats::median(sample_object$nCount_RNA),
    median_percent_mt = stats::median(sample_object$percent.mt)
  )
  rm(sample_object)
  invisible(gc())
}
sn_input_qc <- data.table::rbindlist(sn_input_qc)

sn_object <- merge(
  x = sn_objects[[1L]],
  y = sn_objects[-1L],
  merge.data = FALSE
)
rm(sn_objects)
invisible(gc())

sn_object <- SeuratObject::JoinLayers(sn_object)
sn_object$group <- factor(sn_object$group, levels = c("Young", "Old"))
sn_object <- Seurat::NormalizeData(sn_object, verbose = FALSE)
sn_object <- Seurat::FindVariableFeatures(
  sn_object,
  selection.method = "vst",
  nfeatures = sn_config$variable_features,
  verbose = FALSE
)
sn_object <- Seurat::ScaleData(
  sn_object,
  features = Seurat::VariableFeatures(sn_object),
  verbose = FALSE
)
sn_object <- Seurat::RunPCA(
  sn_object,
  features = Seurat::VariableFeatures(sn_object),
  npcs = max(sn_config$dimensions),
  verbose = FALSE
)
sn_object <- Seurat::RunUMAP(
  sn_object,
  dims = sn_config$dimensions,
  seed.use = sn_config$seed,
  verbose = FALSE
)
sn_object <- Seurat::FindNeighbors(
  sn_object,
  dims = sn_config$dimensions,
  verbose = FALSE
)
sn_object <- Seurat::FindClusters(
  sn_object,
  resolution = sn_config$cluster_resolution,
  random.seed = sn_config$seed,
  verbose = FALSE
)

sn_marker_sets <- list(
  `Fast skeletal fiber` = c("MYH1", "MYH2", "TNNT3", "TNNI2"),
  `Slow skeletal fiber` = c("MYH7", "TNNT1", "TNNI1", "MYL3"),
  `Fibro-adipogenic progenitor` = c(
    "PDGFRA", "COL1A1", "COL1A2", "DCN", "LUM"
  ),
  `Satellite cell` = c("PAX7", "NCAM1", "VCAM1", "MYF5"),
  `Smooth muscle/pericyte` = c("RGS5", "PDGFRB", "CSPG4", "ACTA2", "MCAM"),
  `Endothelial cell` = c("PECAM1", "VWF", "EMCN", "KDR", "EGFL7"),
  `Immune cell` = c("PTPRC", "TYROBP", "LST1", "CD3D", "NKG7")
)
sn_marker_sets <- lapply(
  sn_marker_sets,
  function(markers) intersect(markers, rownames(sn_object))
)
if (any(lengths(sn_marker_sets) < 2L)) {
  stop("Insufficient canonical markers for one or more cell types.")
}

sn_object <- Seurat::AddModuleScore(
  sn_object,
  features = unname(sn_marker_sets),
  name = "celltype_marker_score_",
  seed = sn_config$seed,
  search = FALSE
)
sn_score_columns <- paste0(
  "celltype_marker_score_", seq_along(sn_marker_sets)
)
sn_cluster_scores <- data.table::as.data.table(
  sn_object@meta.data[, c("seurat_clusters", sn_score_columns), drop = FALSE]
)[, lapply(.SD, mean), by = seurat_clusters, .SDcols = sn_score_columns]
sn_cluster_scores_long <- data.table::melt(
  sn_cluster_scores,
  id.vars = "seurat_clusters",
  variable.name = "score_column",
  value.name = "mean_marker_score"
)
sn_cluster_scores_long[, cell_type := names(sn_marker_sets)[
  match(score_column, sn_score_columns)
]]
sn_support_cell_types <- c(
  "Fibro-adipogenic progenitor", "Satellite cell",
  "Smooth muscle/pericyte", "Endothelial cell", "Immune cell"
)
sn_cluster_annotation <- sn_cluster_scores_long[
  , {
    support_scores <- .SD[cell_type %in% sn_support_cell_types]
    myofiber_scores <- .SD[cell_type %in% c(
      "Fast skeletal fiber", "Slow skeletal fiber"
    )]
    best_support <- support_scores[which.max(mean_marker_score)]
    best_myofiber <- myofiber_scores[which.max(mean_marker_score)]
    if (best_support$mean_marker_score >= 0.25) {
      best_support
    } else {
      best_myofiber
    }
  },
  by = seurat_clusters
]
sn_cluster_annotation[, score_column := NULL]
sn_cluster_annotation[, annotation_method := paste0(
  "Canonical-marker module score with support-cell override >=0.25; ",
  "otherwise fast/slow myofiber maximum; resolution ",
  sn_config$cluster_resolution,
  ". This prevents ambient myofiber RNA from masking support-cell identity."
)]

if (!setequal(sn_cluster_annotation$cell_type, names(sn_marker_sets))) {
  stop(
    "Cell annotation did not recover all seven prespecified muscle cell classes."
  )
}

sn_object$cell_type <- sn_cluster_annotation$cell_type[
  match(
    as.character(sn_object$seurat_clusters),
    sn_cluster_annotation$seurat_clusters
  )
]
sn_object$cell_type <- factor(
  sn_object$cell_type,
  levels = names(sn_marker_sets)
)

sn_targets <- c("ABCC8", "SMAD3", "MAPK1", "YWHAZ", "RXRA", "ZBTB7B")
sn_targets_present <- intersect(sn_targets, rownames(sn_object))
if (length(sn_targets_present) == 0L) {
  stop("None of the prespecified Tier A genes are present in GSE167186.")
}

sn_counts <- Seurat::GetAssayData(sn_object, assay = "RNA", layer = "counts")
sn_log_data <- Seurat::GetAssayData(sn_object, assay = "RNA", layer = "data")
sn_cell_groups <- data.table::data.table(
  cell_index = seq_len(ncol(sn_object)),
  sample_name = sn_object$sample_name,
  group = as.character(sn_object$group),
  cell_type = as.character(sn_object$cell_type)
)

sn_expression_summary <- data.table::rbindlist(lapply(
  split(
    sn_cell_groups$cell_index,
    interaction(
      sn_cell_groups$group,
      sn_cell_groups$cell_type,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(cell_indices) {
    data.table::data.table(
      gene = sn_targets_present,
      group = sn_cell_groups$group[cell_indices[1L]],
      cell_type = sn_cell_groups$cell_type[cell_indices[1L]],
      nuclei = length(cell_indices),
      mean_log_normalized_expression = Matrix::rowMeans(
        sn_log_data[sn_targets_present, cell_indices, drop = FALSE]
      ),
      percent_detected = 100 * Matrix::rowMeans(
        sn_counts[sn_targets_present, cell_indices, drop = FALSE] > 0
      )
    )
  }
))

sn_candidate_localization <- sn_expression_summary[
  , .(
    nuclei = sum(nuclei),
    mean_log_normalized_expression = stats::weighted.mean(
      mean_log_normalized_expression, nuclei
    ),
    percent_detected = stats::weighted.mean(percent_detected, nuclei)
  ),
  by = .(gene, cell_type)
]
sn_candidate_localization[
  , localization_rank := data.table::frank(
    -mean_log_normalized_expression, ties.method = "min"
  ),
  by = gene
]

sn_cell_groups[, pseudobulk_id := paste(sample_name, cell_type, sep = "||")]
sn_pseudobulk_levels <- unique(sn_cell_groups$pseudobulk_id)
sn_membership <- Matrix::sparseMatrix(
  i = sn_cell_groups$cell_index,
  j = match(sn_cell_groups$pseudobulk_id, sn_pseudobulk_levels),
  x = 1,
  dims = c(ncol(sn_object), length(sn_pseudobulk_levels)),
  dimnames = list(colnames(sn_object), sn_pseudobulk_levels)
)
sn_pseudobulk_counts <- sn_counts %*% sn_membership
sn_split_ids <- data.table::tstrsplit(
  colnames(sn_pseudobulk_counts), "\\|\\|", fixed = FALSE
)
sn_pseudobulk_meta <- data.table::data.table(
  pseudobulk_id = colnames(sn_pseudobulk_counts),
  sample_name = sn_split_ids[[1L]],
  cell_type = sn_split_ids[[2L]],
  nuclei = as.numeric(Matrix::colSums(sn_membership))
)
sn_pseudobulk_meta <- merge(
  sn_pseudobulk_meta,
  sn_metadata[, .(sample_name, geo_accession, group)],
  by = "sample_name",
  all.x = TRUE,
  sort = FALSE
)

run_candidate_pseudobulk <- function(cell_type_value) {
  cell_meta <- sn_pseudobulk_meta[
    cell_type == cell_type_value &
      nuclei >= sn_config$minimum_nuclei_per_pseudobulk
  ]
  group_counts <- table(cell_meta$group)
  young_donors <- if ("Young" %in% names(group_counts)) {
    unname(group_counts["Young"])
  } else 0L
  old_donors <- if ("Old" %in% names(group_counts)) {
    unname(group_counts["Old"])
  } else 0L

  if (young_donors < 3L || old_donors < 3L) {
    return(data.table::data.table(
      gene = sn_targets,
      cell_type = cell_type_value,
      young_donors = young_donors,
      old_donors = old_donors,
      logFC_Old_vs_Young = NA_real_,
      PValue = NA_real_,
      test_status = "Not tested: fewer than 3 donors in one group"
    ))
  }

  matrix_columns <- match(
    cell_meta$pseudobulk_id, colnames(sn_pseudobulk_counts)
  )
  group_factor <- factor(cell_meta$group, levels = c("Young", "Old"))
  expression_object <- edgeR::DGEList(
    counts = sn_pseudobulk_counts[, matrix_columns, drop = FALSE],
    group = group_factor
  )
  expressed_genes <- edgeR::filterByExpr(
    expression_object, group = group_factor
  )
  expression_object <- expression_object[
    expressed_genes, , keep.lib.sizes = FALSE
  ]
  expression_object <- edgeR::calcNormFactors(expression_object)
  design_matrix <- stats::model.matrix(~ group_factor)
  expression_object <- edgeR::estimateDisp(
    expression_object, design_matrix, robust = TRUE
  )
  fitted_model <- edgeR::glmQLFit(
    expression_object, design_matrix, robust = TRUE
  )
  test_result <- edgeR::glmQLFTest(fitted_model, coef = "group_factorOld")
  full_result <- data.table::as.data.table(
    edgeR::topTags(test_result, n = Inf, sort.by = "none")$table,
    keep.rownames = "gene"
  )

  candidate_result <- merge(
    data.table::data.table(gene = sn_targets),
    full_result[, .(gene, logFC_Old_vs_Young = logFC, PValue)],
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  candidate_result[, `:=`(
    cell_type = cell_type_value,
    young_donors = young_donors,
    old_donors = old_donors,
    test_status = ifelse(
      is.na(PValue), "Not tested: low pseudobulk expression", "Tested"
    )
  )]
  candidate_result[]
}

sn_candidate_de <- data.table::rbindlist(
  lapply(names(sn_marker_sets), run_candidate_pseudobulk),
  fill = TRUE
)
sn_candidate_de[, candidate_FDR_global := stats::p.adjust(
  PValue, method = "BH"
)]
sn_candidate_de[, interpretation := data.table::fcase(
  is.na(PValue), test_status,
  candidate_FDR_global < 0.05 & logFC_Old_vs_Young > 0,
    "Higher in old muscle for this cell type",
  candidate_FDR_global < 0.05 & logFC_Old_vs_Young < 0,
    "Lower in old muscle for this cell type",
  default = "No FDR-significant age association"
)]
data.table::setcolorder(
  sn_candidate_de,
  c(
    "gene", "cell_type", "young_donors", "old_donors",
    "logFC_Old_vs_Young", "PValue", "candidate_FDR_global",
    "test_status", "interpretation"
  )
)

sn_composition_grid <- data.table::CJ(
  sample_name = sn_metadata$sample_name,
  cell_type = names(sn_marker_sets),
  unique = TRUE
)
sn_composition <- sn_cell_groups[, .(nuclei = .N), by = .(sample_name, cell_type)]
sn_composition <- merge(
  sn_composition_grid,
  sn_composition,
  by = c("sample_name", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
sn_composition[is.na(nuclei), nuclei := 0L]
sn_composition <- merge(
  sn_composition,
  sn_metadata[, .(sample_name, group)],
  by = "sample_name",
  all.x = TRUE,
  sort = FALSE
)
sn_composition[, proportion := nuclei / sum(nuclei), by = sample_name]
sn_composition_tests <- sn_composition[
  , {
    young_values <- proportion[group == "Young"]
    old_values <- proportion[group == "Old"]
    list(
      young_donors = length(young_values),
      old_donors = length(old_values),
      young_median_proportion = stats::median(young_values),
      old_median_proportion = stats::median(old_values),
      median_difference_Old_minus_Young =
        stats::median(old_values) - stats::median(young_values),
      PValue = stats::wilcox.test(
        old_values, young_values, exact = FALSE
      )$p.value
    )
  },
  by = cell_type
]
sn_composition_tests[, FDR := stats::p.adjust(PValue, method = "BH")]

sn_umap <- Seurat::DimPlot(
  sn_object,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "group",
  raster = TRUE,
  pt.size = 0.08
) +
  ggplot2::labs(
    title = "GSE167186 human skeletal-muscle single-nucleus atlas",
    subtitle = "Cell identities assigned using prespecified canonical markers"
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank()
  )

ggplot2::ggsave(
  sn_config$output_umap_pdf, sn_umap,
  width = 12, height = 6.5, units = "in"
)
ggplot2::ggsave(
  sn_config$output_umap_png, sn_umap,
  width = 12, height = 6.5, units = "in", dpi = 320
)

sn_dotplot <- ggplot2::ggplot(
  sn_expression_summary,
  ggplot2::aes(
    x = gene,
    y = cell_type,
    size = percent_detected,
    colour = mean_log_normalized_expression
  )
) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ group, ncol = 2) +
  ggplot2::scale_colour_viridis_c(option = "C") +
  ggplot2::scale_size(range = c(1.5, 8)) +
  ggplot2::labs(
    title = "Cell-type localization of prespecified Tier A genes",
    subtitle = "GSE167186; descriptive single-nucleus expression",
    x = NULL,
    y = NULL,
    colour = "Mean log-normalized\nexpression",
    size = "Nuclei detected (%)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggplot2::ggsave(
  sn_config$output_dotplot_pdf, sn_dotplot,
  width = 10, height = 6.5, units = "in"
)
ggplot2::ggsave(
  sn_config$output_dotplot_png, sn_dotplot,
  width = 10, height = 6.5, units = "in", dpi = 320
)

sn_readme <- data.table::data.table(
  item = c(
    "scope", "cohort", "input_matrix", "quality_control",
    "dimension_reduction", "cell_annotation", "candidate_genes",
    "candidate_expression", "differential_expression",
    "multiple_testing", "composition", "interpretation_limit"
  ),
  description = c(
    "Appended cell-type localization and aging analysis; all preceding code and outputs are unchanged.",
    "GSE167186 single-nucleus cohort: 11 old and 6 young human vastus-lateralis donors; this is not a sarcopenia-versus-age-matched-control contrast.",
    "Processed per-sample CSV count matrices from GSE167186_RAW.tar are used because they contain the authors' 143,051 post-QC nuclei; the accompanying H5 files do not reproduce that processed matrix.",
    "The authors' processed matrices are retained with the published minimum of 200 detected genes per nucleus.",
    "Log normalization, 2,000 variable genes, scaling, 10 PCs, UMAP, and Louvain clustering at resolution 0.2 reproduce the published framework; the donor-level age model has no batch term because usable batch covariates are unavailable. All 17 donors are male, so sex cannot be estimated. Residual confounding remains possible.",
    "Clusters are assigned to seven prespecified skeletal-muscle cell classes using canonical-marker module scores. A support-cell score of at least 0.25 overrides ambient myofiber signal; otherwise the higher fast/slow myofiber score is used. All scores are exported for audit.",
    "Only ABCC8, SMAD3, MAPK1, YWHAZ, RXRA, and ZBTB7B are used for candidate-set inference.",
    "Mean log-normalized expression and percentage of nuclei detected are descriptive localization summaries, not independent-sample significance tests.",
    "Raw counts are summed within donor and cell type. edgeR quasi-likelihood models compare Old versus Young using donors, not nuclei, as biological replicates; cell types require at least 20 nuclei per donor and at least 3 donors per group.",
    "BH FDR is calculated globally across all tested candidate-gene-by-cell-type comparisons.",
    "Cell-type proportions are calculated per donor and compared by Wilcoxon tests with BH FDR across seven cell types.",
    "Aging-associated cell localization can strengthen biological plausibility but cannot establish sarcopenia-specific replication, mediation, causality, or clinical efficacy."
  )
)

sn_output_tables <- list(
  `01_sample_QC` = as.data.frame(sn_input_qc),
  `02_cluster_annotation` = as.data.frame(sn_cluster_annotation),
  `03_all_marker_scores` = as.data.frame(sn_cluster_scores_long),
  `04_cell_composition` = as.data.frame(sn_composition),
  `05_composition_tests` = as.data.frame(sn_composition_tests),
  `06_TierA_expression` = as.data.frame(sn_expression_summary),
  `07_TierA_localization` = as.data.frame(sn_candidate_localization),
  `08_TierA_pseudobulk_DE` = as.data.frame(sn_candidate_de),
  `09_pseudobulk_QC` = as.data.frame(sn_pseudobulk_meta),
  `10_README` = as.data.frame(sn_readme)
)

openxlsx::write.xlsx(
  sn_output_tables,
  file = sn_config$output_workbook,
  overwrite = TRUE,
  asTable = TRUE
)
sn_workbook <- openxlsx::loadWorkbook(sn_config$output_workbook)
for (sheet_name in names(sn_output_tables)) {
  sheet_data <- sn_output_tables[[sheet_name]]
  openxlsx::freezePane(sn_workbook, sheet = sheet_name, firstRow = TRUE)
  openxlsx::showGridLines(
    sn_workbook, sheet = sheet_name, showGridLines = FALSE
  )
  if (ncol(sheet_data) > 0L) {
    widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- head(as.character(sheet_data[[column_index]]), 500L)
      observed <- if (length(values) > 0L) {
        max(nchar(values), na.rm = TRUE)
      } else 0L
      min(
        28,
        max(
          10,
          nchar(names(sheet_data)[column_index]) + 2L,
          observed + 2L
        )
      )
    }, numeric(1))
    if (identical(sheet_name, "10_README")) {
      widths <- c(24, 90)
      openxlsx::addStyle(
        sn_workbook,
        sheet = sheet_name,
        style = openxlsx::createStyle(wrapText = TRUE, valign = "top"),
        rows = 2:(nrow(sheet_data) + 1L),
        cols = 2L,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
    openxlsx::setColWidths(
      sn_workbook,
      sheet = sheet_name,
      cols = seq_len(ncol(sheet_data)),
      widths = widths
    )
  }
}
openxlsx::saveWorkbook(
  sn_workbook, sn_config$output_workbook, overwrite = TRUE
)
sn_workbook_repair <- repair_missing_ooxml_relationships(
  sn_config$output_workbook
)

sn_actual_sheets <- openxlsx::getSheetNames(sn_config$output_workbook)
if (!identical(sn_actual_sheets, names(sn_output_tables))) {
  stop("GSE167186 single-nucleus workbook validation failed.")
}

log_step(
  "GSE167186 single-nucleus analysis complete: ",
  sn_config$output_workbook,
  "; donors = ", nrow(sn_metadata),
  "; nuclei = ", ncol(sn_object),
  "; annotated cell types = ", length(unique(sn_object$cell_type)),
  "; FDR-significant candidate-cell-type associations = ",
  sum(sn_candidate_de$candidate_FDR_global < 0.05, na.rm = TRUE)
)



# 19. Public-database clinical translation assessment--------------------
# HPA + Open Targets + ChEMBL + ClinicalTrials.gov
# Appended analysis: the preceding code is intentionally unchanged.


clinical_translation_packages <- c(
  "httr2", "jsonlite", "xml2", "openxlsx", "data.table",
  "AnnotationDbi", "org.Hs.eg.db")

clinical_translation_missing <- clinical_translation_packages[
  !vapply(clinical_translation_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(clinical_translation_missing) > 0L) {
  stop(
    "Section 19 requires installed packages: ",
    paste(clinical_translation_missing, collapse = ", "))}

ct_config <- list(
  output_workbook = file.path(
    "results_public_release", "clinical_translation_public_databases.xlsx"
  ),
  discovery_workbook = file.path(
    "results_public_release", "GEO_human_muscle_validation.xlsx"
  ),
  candidate_genes = c("ABCC8", "MAPK1", "RXRA", "SMAD3", "YWHAZ", "ZBTB7B"),
  candidate_ensembl = c(
    "ENSG00000006071", "ENSG00000100030", "ENSG00000186350",
    "ENSG00000166949", "ENSG00000164924", "ENSG00000160685"
  ),
  clinical_context = paste(
    "sarcopenia OR muscle weakness OR muscle atrophy OR frailty OR",
    "insulin resistance OR diabetes"
  ),
  max_trial_drugs_per_gene = 5L,
  max_trials_per_query = 100L
)
dir.create(dirname(ct_config$output_workbook), recursive = TRUE, showWarnings = FALSE)

ct_query_log <- new.env(parent = emptyenv())
ct_query_log$rows <- list()
ct_log_query <- function(source, query_label, status, detail = "") {
  ct_query_log$rows[[length(ct_query_log$rows) + 1L]] <- data.frame(
    source = source,
    query_label = query_label,
    status = status,
    detail = substr(as.character(detail), 1L, 1000L),
    queried_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
}

ct_text <- function(x, collapse = "; ") {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  values <- unlist(x, recursive = TRUE, use.names = FALSE)
  values <- trimws(as.character(values))
  values <- unique(values[!is.na(values) & nzchar(values)])
  if (length(values) == 0L) NA_character_ else paste(values, collapse = collapse)
}

ct_num <- function(x) {
  value <- suppressWarnings(as.numeric(ct_text(x)))
  if (length(value) == 0L || is.na(value[1L])) NA_real_ else value[1L]
}

ct_empty_table <- function(columns) {
  as.data.frame(setNames(replicate(
    length(columns), character(0), simplify = FALSE
  ), columns), stringsAsFactors = FALSE)
}

ct_get_json <- function(url, query = list(), source, query_label) {
  tryCatch({
    parsed_url <- httr2::url_parse(url)
    parsed_url$query <- c(parsed_url$query, query)
    request <- httr2::request(httr2::url_build(parsed_url))
    request <- httr2::req_user_agent(
      request, "IR-sarcopenia-public-database-analysis/1.0"
    )
    request <- httr2::req_timeout(request, seconds = 60)
    request <- httr2::req_retry(request, max_tries = 3)
    response <- httr2::req_perform(request)
    httr2::resp_check_status(response)
    result <- httr2::resp_body_json(response, simplifyVector = FALSE)
    ct_log_query(source, query_label, "success")
    result
  }, error = function(error) {
    ct_log_query(source, query_label, "failed", conditionMessage(error))
    NULL
  })
}

ct_get_text <- function(url, source, query_label) {
  tryCatch({
    request <- httr2::request(url)
    request <- httr2::req_user_agent(
      request, "IR-sarcopenia-public-database-analysis/1.0"
    )
    request <- httr2::req_timeout(request, seconds = 60)
    request <- httr2::req_retry(request, max_tries = 3)
    response <- httr2::req_perform(request)
    httr2::resp_check_status(response)
    result <- httr2::resp_body_string(response)
    ct_log_query(source, query_label, "success")
    result
  }, error = function(error) {
    ct_log_query(source, query_label, "failed", conditionMessage(error))
    NULL
  })
}

ct_post_graphql <- function(query_string, variables, query_label) {
  tryCatch({
    request <- httr2::request(
      "https://api.platform.opentargets.org/api/v4/graphql"
    )
    request <- httr2::req_user_agent(
      request, "IR-sarcopenia-public-database-analysis/1.0"
    )
    request <- httr2::req_body_json(
      request, list(query = query_string, variables = variables)
    )
    request <- httr2::req_timeout(request, seconds = 60)
    request <- httr2::req_retry(request, max_tries = 3)
    response <- httr2::req_perform(request)
    httr2::resp_check_status(response)
    result <- httr2::resp_body_json(response, simplifyVector = FALSE)
    if (!is.null(result$errors)) {
      stop(ct_text(vapply(result$errors, function(x) x$message, character(1))))
    }
    ct_log_query("Open Targets", query_label, "success")
    result$data$target
  }, error = function(error) {
    ct_log_query(
      "Open Targets", query_label, "failed", conditionMessage(error)
    )
    NULL
  })
}

ct_candidates <- data.frame(
  gene_symbol = ct_config$candidate_genes,
  gene_ensembl = ct_config$candidate_ensembl,
  stringsAsFactors = FALSE
)

ct_annotation <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = ct_config$candidate_genes,
  keytype = "SYMBOL",
  columns = c("ENSEMBL", "ENTREZID", "UNIPROT")
)
ct_annotation <- as.data.frame(ct_annotation, stringsAsFactors = FALSE)
ct_annotation <- ct_annotation[!duplicated(ct_annotation$SYMBOL), , drop = FALSE]
ct_candidates$entrez_id <- ct_annotation$ENTREZID[
  match(ct_candidates$gene_symbol, ct_annotation$SYMBOL)
]
ct_candidates$uniprot_orgdb <- ct_annotation$UNIPROT[
  match(ct_candidates$gene_symbol, ct_annotation$SYMBOL)
]

# Recover the direction supported by the existing MR/colocalisation results.
ct_context <- openxlsx::read.xlsx(
  ct_config$discovery_workbook, sheet = "08_discovery_context"
)
ct_context <- as.data.frame(ct_context, stringsAsFactors = FALSE)
ct_context <- ct_context[
  ct_context$gene_symbol %in% ct_config$candidate_genes &
    !is.na(ct_context$discovery_MR_FDR) &
    ct_context$discovery_MR_FDR < 0.05 &
    !is.na(ct_context$discovery_coloc_PP4),
  , drop = FALSE
]
ct_context <- ct_context[order(
  ct_context$gene_symbol,
  -ct_context$discovery_coloc_PP4,
  ct_context$discovery_MR_FDR
), , drop = FALSE]
ct_context <- ct_context[!duplicated(ct_context$gene_symbol), , drop = FALSE]
ct_context$phenotype_interpretation <- ifelse(
  ct_context$trait %in% c("ALM", "GRIP", "WALK"),
  "Higher phenotype value is treated as favourable",
  "Trait direction requires phenotype-specific interpretation"
)
ct_context$genetically_supported_expression_direction <- ifelse(
  ct_context$trait %in% c("ALM", "GRIP", "WALK") &
    ct_context$discovery_MR_beta > 0,
  "increase expression",
  ifelse(
    ct_context$trait %in% c("ALM", "GRIP", "WALK") &
      ct_context$discovery_MR_beta < 0,
    "decrease expression",
    "uncertain"
  )
)
ct_context$direction_caveat <- paste(
  "Expression-direction inference is not equivalent to acute protein",
  "agonism, antagonism, inhibition, activation, dose response, or clinical efficacy."
)

# Human Protein Atlas: summary annotation and skeletal-muscle IHC evidence.
ct_hpa_summary <- list()
ct_hpa_muscle <- list()
for (candidate_index in seq_len(nrow(ct_candidates))) {
  gene_symbol <- ct_candidates$gene_symbol[candidate_index]
  ensembl_id <- ct_candidates$gene_ensembl[candidate_index]
  hpa_json_url <- paste0("https://www.proteinatlas.org/", ensembl_id, ".json")
  hpa_xml_url <- paste0("https://www.proteinatlas.org/", ensembl_id, ".xml")
  hpa_json <- ct_get_json(
    hpa_json_url, source = "Human Protein Atlas",
    query_label = paste(gene_symbol, "JSON")
  )
  if (is.null(hpa_json)) hpa_json <- list()
  ct_hpa_summary[[length(ct_hpa_summary) + 1L]] <- data.frame(
    gene_symbol = gene_symbol,
    gene_ensembl = ensembl_id,
    hpa_gene = ct_text(hpa_json[["Gene"]]),
    uniprot = ct_text(hpa_json[["Uniprot"]]),
    gene_description = ct_text(hpa_json[["Gene description"]]),
    protein_class = ct_text(hpa_json[["Protein class"]]),
    biological_process = ct_text(hpa_json[["Biological process"]]),
    molecular_function = ct_text(hpa_json[["Molecular function"]]),
    disease_involvement = ct_text(hpa_json[["Disease involvement"]]),
    evidence = ct_text(hpa_json[["Evidence"]]),
    rna_tissue_specificity = ct_text(hpa_json[["RNA tissue specificity"]]),
    rna_tissue_distribution = ct_text(hpa_json[["RNA tissue distribution"]]),
    rna_single_cell_specificity = ct_text(
      hpa_json[["RNA single cell type specificity"]]
    ),
    protein_cell_type_specificity = ct_text(
      hpa_json[["Protein cell type specificity"]]
    ),
    protein_tissue_specificity = ct_text(
      hpa_json[["Protein tissue specificity"]]
    ),
    subcellular_location = ct_text(hpa_json[["Subcellular location"]]),
    hpa_url = paste0("https://www.proteinatlas.org/", ensembl_id),
    stringsAsFactors = FALSE
  )

  hpa_xml_text <- ct_get_text(
    hpa_xml_url, source = "Human Protein Atlas",
    query_label = paste(gene_symbol, "XML")
  )
  muscle_rows_before <- length(ct_hpa_muscle)
  if (!is.null(hpa_xml_text)) {
    hpa_xml <- xml2::read_xml(hpa_xml_text)
    muscle_nodes <- xml2::xml_find_all(
      hpa_xml,
      paste0(
        "//*[local-name()='tissueExpression']",
        "[@technology='IHC' and @assayType='tissue']",
        "/*[local-name()='data']",
        "[*[local-name()='tissue'][normalize-space(text())='Skeletal muscle']]"
      )
    )
    for (muscle_node in muscle_nodes) {
      antibody_node <- xml2::xml_find_first(
        muscle_node, "ancestor::*[local-name()='antibody'][1]"
      )
      tissue_level <- xml2::xml_text(xml2::xml_find_first(
        muscle_node, "./*[local-name()='level'][@type='expression']"
      ))
      tissue_cells <- xml2::xml_find_all(
        muscle_node, "./*[local-name()='tissueCell']"
      )
      if (length(tissue_cells) == 0L) tissue_cells <- xml2::xml_find_all(
        muscle_node, "."
      )
      for (cell_index in seq_along(tissue_cells)) {
        tissue_cell <- tissue_cells[[cell_index]]
        cell_type <- xml2::xml_text(xml2::xml_find_first(
          tissue_cell, "./*[local-name()='cellType']"
        ))
        cell_level <- xml2::xml_text(xml2::xml_find_first(
          tissue_cell, "./*[local-name()='level'][@type='expression']"
        ))
        ct_hpa_muscle[[length(ct_hpa_muscle) + 1L]] <- data.frame(
          gene_symbol = gene_symbol,
          gene_ensembl = ensembl_id,
          tissue = "Skeletal muscle",
          cell_type = ifelse(nzchar(cell_type), cell_type, NA_character_),
          tissue_expression_level = ifelse(
            nzchar(tissue_level), tissue_level, NA_character_
          ),
          cell_expression_level = ifelse(
            nzchar(cell_level), cell_level, NA_character_
          ),
          antibody_id = xml2::xml_attr(antibody_node, "id"),
          antibody_reliability = xml2::xml_attr(antibody_node, "reliability"),
          evidence_type = "HPA immunohistochemistry",
          source_url = hpa_xml_url,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(ct_hpa_muscle) == muscle_rows_before) {
    ct_hpa_muscle[[length(ct_hpa_muscle) + 1L]] <- data.frame(
      gene_symbol = gene_symbol, gene_ensembl = ensembl_id,
      tissue = "Skeletal muscle", cell_type = NA_character_,
      tissue_expression_level = "not reported",
      cell_expression_level = "not reported", antibody_id = NA_character_,
      antibody_reliability = NA_character_,
      evidence_type = "HPA immunohistochemistry", source_url = hpa_xml_url,
      stringsAsFactors = FALSE
    )
  }
}
ct_hpa_summary <- data.table::rbindlist(ct_hpa_summary, fill = TRUE)
ct_hpa_muscle <- data.table::rbindlist(ct_hpa_muscle, fill = TRUE)
ct_candidates$uniprot <- ct_hpa_summary$uniprot[
  match(ct_candidates$gene_symbol, ct_hpa_summary$gene_symbol)
]
ct_candidates$uniprot[is.na(ct_candidates$uniprot)] <-
  ct_candidates$uniprot_orgdb[is.na(ct_candidates$uniprot)]

# Open Targets: relevant disease links, tractability and safety annotations.
ct_ot_query <- paste(
  "query targetClinicalTranslation($ensemblId: String!) {",
  "target(ensemblId: $ensemblId) {",
  "id approvedSymbol approvedName",
  "tractability { label modality value }",
  "safetyLiabilities { event datasource literature",
  "biosamples { tissueLabel cellLabel } }",
  "associatedDiseases(page: {index: 0, size: 500}) {",
  "count rows { score disease { id name } } }",
  "} }"
)
ct_relevant_pattern <- paste0(
  "sarcopen|frailty|muscle weakness|muscle atrophy|muscle mass|",
  "insulin resistan|diabetes|glucose|hypogly|obesity|body mass|",
  "walking|grip strength"
)
ct_ot_associations <- list()
ct_ot_tractability <- list()
ct_ot_safety <- list()
for (candidate_index in seq_len(nrow(ct_candidates))) {
  gene_symbol <- ct_candidates$gene_symbol[candidate_index]
  ensembl_id <- ct_candidates$gene_ensembl[candidate_index]
  target <- ct_post_graphql(
    ct_ot_query, list(ensemblId = ensembl_id), gene_symbol
  )
  if (is.null(target)) next
  association_rows <- target$associatedDiseases$rows
  if (length(association_rows) > 0L) {
    for (association in association_rows) {
      disease_name <- ct_text(association$disease$name)
      if (!is.na(disease_name) && grepl(
        ct_relevant_pattern, disease_name, ignore.case = TRUE
      )) {
        ct_ot_associations[[length(ct_ot_associations) + 1L]] <- data.frame(
          gene_symbol = gene_symbol,
          gene_ensembl = ensembl_id,
          disease_id = ct_text(association$disease$id),
          disease_name = disease_name,
          association_score = ct_num(association$score),
          relevance_rule = "Prespecified muscle/metabolic keyword filter",
          source_url = paste0(
            "https://platform.opentargets.org/target/", ensembl_id,
            "/associations"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(target$tractability) > 0L) {
    for (tractability in target$tractability) {
      ct_ot_tractability[[length(ct_ot_tractability) + 1L]] <- data.frame(
        gene_symbol = gene_symbol, gene_ensembl = ensembl_id,
        modality = ct_text(tractability$modality),
        assessment = ct_text(tractability$label),
        supported = isTRUE(tractability$value),
        source = "Open Targets tractability",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(target$safetyLiabilities) > 0L) {
    for (safety in target$safetyLiabilities) {
      ct_ot_safety[[length(ct_ot_safety) + 1L]] <- data.frame(
        gene_symbol = gene_symbol, gene_ensembl = ensembl_id,
        safety_event = ct_text(safety$event),
        datasource = ct_text(safety$datasource),
        tissue = ct_text(lapply(safety$biosamples, `[[`, "tissueLabel")),
        cell = ct_text(lapply(safety$biosamples, `[[`, "cellLabel")),
        literature = ct_text(safety$literature),
        interpretation = paste(
          "Safety annotation is a flag for review, not proof that a",
          "candidate-targeted intervention is unsafe."
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
ct_ot_associations <- data.table::rbindlist(ct_ot_associations, fill = TRUE)
ct_ot_tractability <- data.table::rbindlist(ct_ot_tractability, fill = TRUE)
ct_ot_safety <- data.table::rbindlist(ct_ot_safety, fill = TRUE)
if (ncol(ct_ot_associations) == 0L) ct_ot_associations <- ct_empty_table(c(
  "gene_symbol", "gene_ensembl", "disease_id", "disease_name",
  "association_score", "relevance_rule", "source_url"
))
if (ncol(ct_ot_tractability) == 0L) ct_ot_tractability <- ct_empty_table(c(
  "gene_symbol", "gene_ensembl", "modality", "assessment", "supported", "source"
))
ct_ot_tractability$supported <- ct_ot_tractability$supported %in% TRUE
if (ncol(ct_ot_safety) == 0L) ct_ot_safety <- ct_empty_table(c(
  "gene_symbol", "gene_ensembl", "safety_event", "datasource", "tissue",
  "cell", "literature", "interpretation"
))

# ChEMBL: all human targets containing the candidate protein, including complexes.
ct_chembl_targets <- list()
ct_chembl_mechanisms <- list()
ct_molecule_cache <- new.env(parent = emptyenv())
ct_get_molecule <- function(molecule_id) {
  if (exists(molecule_id, envir = ct_molecule_cache, inherits = FALSE)) {
    return(get(molecule_id, envir = ct_molecule_cache, inherits = FALSE))
  }
  molecule <- ct_get_json(
    paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", molecule_id, ".json"),
    source = "ChEMBL", query_label = paste("molecule", molecule_id)
  )
  if (is.null(molecule)) molecule <- list()
  assign(molecule_id, molecule, envir = ct_molecule_cache)
  molecule
}

for (candidate_index in seq_len(nrow(ct_candidates))) {
  gene_symbol <- ct_candidates$gene_symbol[candidate_index]
  accession <- ct_candidates$uniprot[candidate_index]
  if (is.na(accession) || !nzchar(accession)) {
    ct_log_query("ChEMBL", gene_symbol, "skipped", "No UniProt accession")
    next
  }
  target_result <- ct_get_json(
    "https://www.ebi.ac.uk/chembl/api/data/target.json",
    query = list(target_components__accession = accession, limit = 100L),
    source = "ChEMBL", query_label = paste(gene_symbol, "targets")
  )
  if (is.null(target_result) || length(target_result$targets) == 0L) next
  human_targets <- Filter(
    function(target) identical(ct_text(target$organism), "Homo sapiens"),
    target_result$targets
  )
  for (target in human_targets) {
    target_id <- ct_text(target$target_chembl_id)
    ct_chembl_targets[[length(ct_chembl_targets) + 1L]] <- data.frame(
      gene_symbol = gene_symbol, uniprot = accession,
      target_chembl_id = target_id,
      target_name = ct_text(target$pref_name),
      target_type = ct_text(target$target_type),
      complex_context = !identical(
        ct_text(target$target_type), "SINGLE PROTEIN"
      ),
      source_url = paste0(
        "https://www.ebi.ac.uk/chembl/explore/target/", target_id
      ),
      stringsAsFactors = FALSE
    )
    mechanism_result <- ct_get_json(
      "https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
      query = list(target_chembl_id = target_id, limit = 1000L),
      source = "ChEMBL",
      query_label = paste(gene_symbol, target_id, "mechanisms")
    )
    if (is.null(mechanism_result) || length(mechanism_result$mechanisms) == 0L) {
      next
    }
    for (mechanism in mechanism_result$mechanisms) {
      molecule_id <- ct_text(mechanism$molecule_chembl_id)
      molecule <- ct_get_molecule(molecule_id)
      molecule_name <- ct_text(molecule$pref_name)
      if (is.na(molecule_name)) molecule_name <- molecule_id
      ct_chembl_mechanisms[[length(ct_chembl_mechanisms) + 1L]] <- data.frame(
        gene_symbol = gene_symbol, uniprot = accession,
        target_chembl_id = target_id,
        target_name = ct_text(target$pref_name),
        target_type = ct_text(target$target_type),
        molecule_chembl_id = molecule_id,
        molecule_name = molecule_name,
        action_type = ct_text(mechanism$action_type),
        mechanism_of_action = ct_text(mechanism$mechanism_of_action),
        direct_interaction = isTRUE(as.logical(mechanism$direct_interaction)),
        disease_efficacy = isTRUE(as.logical(mechanism$disease_efficacy)),
        max_phase = ct_num(mechanism$max_phase),
        first_approval = ct_text(molecule$first_approval),
        molecule_type = ct_text(molecule$molecule_type),
        withdrawn_flag = isTRUE(as.logical(molecule$withdrawn_flag)),
        source_url = paste0(
          "https://www.ebi.ac.uk/chembl/explore/compound/", molecule_id
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
ct_chembl_targets <- unique(data.table::rbindlist(ct_chembl_targets, fill = TRUE))
ct_chembl_mechanisms <- unique(data.table::rbindlist(
  ct_chembl_mechanisms, fill = TRUE
))
if (ncol(ct_chembl_targets) == 0L) ct_chembl_targets <- ct_empty_table(c(
  "gene_symbol", "uniprot", "target_chembl_id", "target_name", "target_type",
  "complex_context", "source_url"
))
if (ncol(ct_chembl_mechanisms) == 0L) ct_chembl_mechanisms <- ct_empty_table(c(
  "gene_symbol", "uniprot", "target_chembl_id", "target_name", "target_type",
  "molecule_chembl_id", "molecule_name", "action_type", "mechanism_of_action",
  "direct_interaction", "disease_efficacy", "max_phase", "first_approval",
  "molecule_type", "withdrawn_flag", "source_url"
))

# Compare genetic expression direction with reported pharmacological action.
ct_direction_check <- merge(
  as.data.frame(ct_chembl_mechanisms),
  ct_context[, c(
    "gene_symbol", "trait", "discovery_MR_beta", "discovery_MR_FDR",
    "discovery_coloc_PP4", "genetically_supported_expression_direction"
  )],
  by = "gene_symbol", all.x = TRUE
)
if (nrow(ct_direction_check) > 0L) {
  activating_pattern <- "AGONIST|ACTIVATOR|STIMULATOR|POSITIVE MODULATOR"
  inhibiting_pattern <- "INHIBITOR|ANTAGONIST|BLOCKER|NEGATIVE MODULATOR"
  ct_direction_check$pharmacology_direction <- ifelse(
    grepl(activating_pattern, ct_direction_check$action_type, ignore.case = TRUE),
    "increase activity",
    ifelse(
      grepl(inhibiting_pattern, ct_direction_check$action_type, ignore.case = TRUE),
      "decrease activity", "uncertain"
    )
  )
  ct_direction_check$direction_concordance <- ifelse(
    ct_direction_check$genetically_supported_expression_direction ==
      "increase expression" &
      ct_direction_check$pharmacology_direction == "increase activity",
    "directionally concordant",
    ifelse(
      ct_direction_check$genetically_supported_expression_direction ==
        "decrease expression" &
        ct_direction_check$pharmacology_direction == "decrease activity",
      "directionally concordant",
      ifelse(
        ct_direction_check$pharmacology_direction == "uncertain" |
          ct_direction_check$genetically_supported_expression_direction == "uncertain",
        "uncertain", "directionally discordant"
      )
    )
  )
  ct_direction_check$interpretation_limit <- paste(
    "Concordance is a hypothesis-generating comparison only; expression MR",
    "does not establish the effect of acute or tissue-specific drug modulation."
  )
}
if (nrow(ct_direction_check) == 0L) ct_direction_check <- ct_empty_table(c(
  "gene_symbol", "molecule_name", "action_type", "mechanism_of_action",
  "genetically_supported_expression_direction", "pharmacology_direction",
  "direction_concordance", "interpretation_limit"
))

# ClinicalTrials.gov: gene-context search plus a bounded search for clinical drugs.
ct_trial_queries <- list()
for (candidate_gene in ct_config$candidate_genes) {
  ct_trial_queries[[length(ct_trial_queries) + 1L]] <- data.frame(
    gene_symbol = candidate_gene,
    search_term = candidate_gene,
    match_basis = "gene symbol",
    target_linkage = paste(
      "Gene mentioned in record; a target intervention is not established"
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(ct_chembl_mechanisms) > 0L) {
    drug_rows <- ct_chembl_mechanisms[
      ct_chembl_mechanisms$gene_symbol == candidate_gene &
        !is.na(ct_chembl_mechanisms$molecule_name) &
        nzchar(ct_chembl_mechanisms$molecule_name) &
        !is.na(ct_chembl_mechanisms$max_phase) &
        ct_chembl_mechanisms$max_phase >= 1,
      , drop = FALSE
    ]
    if (nrow(drug_rows) > 0L) {
      drug_rows <- drug_rows[order(-drug_rows$max_phase), , drop = FALSE]
      drug_rows <- drug_rows[!duplicated(drug_rows$molecule_name), , drop = FALSE]
      drug_rows <- head(drug_rows, ct_config$max_trial_drugs_per_gene)
      for (drug_index in seq_len(nrow(drug_rows))) {
        ct_trial_queries[[length(ct_trial_queries) + 1L]] <- data.frame(
          gene_symbol = candidate_gene,
          search_term = drug_rows$molecule_name[drug_index],
          match_basis = "ChEMBL clinical molecule",
          target_linkage = paste0(
            "ChEMBL mechanism links molecule to ",
            drug_rows$target_name[drug_index], " (",
            drug_rows$action_type[drug_index], ")"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
ct_trial_queries <- unique(data.table::rbindlist(ct_trial_queries, fill = TRUE))
ct_clinical_trials <- list()
for (trial_query_index in seq_len(nrow(ct_trial_queries))) {
  query_row <- ct_trial_queries[trial_query_index, ]
  query_expression <- paste0(
    '"', query_row$search_term, '" AND (', ct_config$clinical_context, ')'
  )
  trial_result <- ct_get_json(
    "https://clinicaltrials.gov/api/v2/studies",
    query = list(
      `query.term` = query_expression,
      pageSize = ct_config$max_trials_per_query,
      countTotal = "true",
      format = "json"
    ),
    source = "ClinicalTrials.gov",
    query_label = paste(query_row$gene_symbol, query_row$search_term)
  )
  if (is.null(trial_result) || length(trial_result$studies) == 0L) next
  for (study in trial_result$studies) {
    protocol <- study$protocolSection
    identification <- protocol$identificationModule
    status <- protocol$statusModule
    design <- protocol$designModule
    conditions <- protocol$conditionsModule$conditions
    interventions <- protocol$armsInterventionsModule$interventions
    outcomes <- protocol$outcomesModule$primaryOutcomes
    sponsors <- protocol$sponsorCollaboratorsModule
    ct_clinical_trials[[length(ct_clinical_trials) + 1L]] <- data.frame(
      gene_symbol = query_row$gene_symbol,
      matched_term = query_row$search_term,
      match_basis = query_row$match_basis,
      target_linkage = query_row$target_linkage,
      nct_id = ct_text(identification$nctId),
      brief_title = ct_text(identification$briefTitle),
      overall_status = ct_text(status$overallStatus),
      study_type = ct_text(design$studyType),
      phase = ct_text(design$phases),
      enrollment = ct_num(design$enrollmentInfo$count),
      start_date = ct_text(status$startDateStruct$date),
      completion_date = ct_text(status$completionDateStruct$date),
      conditions = ct_text(conditions),
      interventions = ct_text(lapply(interventions, function(x) {
        paste(ct_text(x$type), ct_text(x$name), sep = ": ")
      })),
      primary_outcomes = ct_text(lapply(outcomes, `[[`, "measure")),
      lead_sponsor = ct_text(sponsors$leadSponsor$name),
      has_results = isTRUE(study$hasResults),
      records_returned_for_query = length(trial_result$studies),
      total_records_for_query = ct_num(trial_result$totalCount),
      retrieval_cap = ct_config$max_trials_per_query,
      source_url = paste0(
        "https://clinicaltrials.gov/study/", ct_text(identification$nctId)
      ),
      stringsAsFactors = FALSE
    )
  }
}
ct_clinical_trials <- unique(data.table::rbindlist(
  ct_clinical_trials, fill = TRUE
))
if (ncol(ct_clinical_trials) == 0L) ct_clinical_trials <- ct_empty_table(c(
  "gene_symbol", "matched_term", "match_basis", "target_linkage", "nct_id",
  "brief_title", "overall_status", "study_type", "phase", "enrollment",
  "start_date", "completion_date", "conditions", "interventions",
  "primary_outcomes", "lead_sponsor", "has_results",
  "records_returned_for_query", "total_records_for_query", "retrieval_cap",
  "source_url"
))

# Gene-level evidence matrix; this is not a validated clinical prediction score.
ct_gene_priority <- ct_candidates[, c(
  "gene_symbol", "gene_ensembl", "entrez_id", "uniprot"
)]
ct_gene_priority$genetic_trait <- ct_context$trait[
  match(ct_gene_priority$gene_symbol, ct_context$gene_symbol)
]
ct_gene_priority$expression_direction <-
  ct_context$genetically_supported_expression_direction[
    match(ct_gene_priority$gene_symbol, ct_context$gene_symbol)
  ]
ct_gene_priority$hpa_skeletal_muscle_protein_reported <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) any(
    ct_hpa_muscle$gene_symbol == gene &
      !ct_hpa_muscle$tissue_expression_level %in% c("not reported", NA)
  ), logical(1)
)
ct_gene_priority$max_relevant_ot_association <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) {
    values <- ct_ot_associations$association_score[
      ct_ot_associations$gene_symbol == gene
    ]
    if (length(values) == 0L || all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
  }, numeric(1)
)
ct_gene_priority$open_targets_supported_tractability_items <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) sum(
    ct_ot_tractability$gene_symbol == gene & ct_ot_tractability$supported,
    na.rm = TRUE
  ), integer(1)
)
ct_gene_priority$chembl_mechanism_count <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) sum(ct_chembl_mechanisms$gene_symbol == gene), integer(1)
)
ct_gene_priority$highest_chembl_phase <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) {
    values <- ct_chembl_mechanisms$max_phase[
      ct_chembl_mechanisms$gene_symbol == gene
    ]
    if (length(values) == 0L || all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
  }, numeric(1)
)
ct_gene_priority$directionally_concordant_mechanisms <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) sum(
    ct_direction_check$gene_symbol == gene &
      ct_direction_check$direction_concordance == "directionally concordant",
    na.rm = TRUE
  ), integer(1)
)
ct_gene_priority$drug_linked_trial_records <- vapply(
  ct_gene_priority$gene_symbol,
  function(gene) length(unique(ct_clinical_trials$nct_id[
    ct_clinical_trials$gene_symbol == gene &
      ct_clinical_trials$match_basis == "ChEMBL clinical molecule" &
      !is.na(ct_clinical_trials$nct_id)
  ])), integer(1)
)
ct_gene_priority$translation_evidence_level <- ifelse(
  ct_gene_priority$directionally_concordant_mechanisms > 0 &
    ct_gene_priority$drug_linked_trial_records > 0,
  "Level 1: direction-concordant mechanism plus drug-linked trial record",
  ifelse(
    ct_gene_priority$chembl_mechanism_count > 0,
    "Level 2: ChEMBL mechanism; direction or trial linkage incomplete",
    ifelse(
      ct_gene_priority$open_targets_supported_tractability_items > 0,
      "Level 3: tractability annotation without identified ChEMBL mechanism",
      "Level 4: limited public-database translation evidence"
    )
  )
)
ct_gene_priority$interpretation <- paste(
  "Evidence level organizes follow-up only and is not a validated prediction",
  "of efficacy, safety, or benefit in sarcopenia."
)

ct_query_log_table <- data.table::rbindlist(ct_query_log$rows, fill = TRUE)
ct_sources <- data.frame(
  database = c(
    "Human Protein Atlas", "Open Targets Platform", "ChEMBL",
    "ClinicalTrials.gov"
  ),
  access_method = c(
    "Individual-entry JSON and XML", "GraphQL API v4",
    "ChEMBL REST API", "ClinicalTrials.gov API v2"
  ),
  official_url = c(
    "https://www.proteinatlas.org/about/help/dataaccess",
    "https://platform-docs.opentargets.org/data-access/graphql-api",
    "https://www.ebi.ac.uk/chembl/api/data/docs",
    "https://clinicaltrials.gov/data-api/about-api"
  ),
  role = c(
    "Human tissue, cell-type and protein-expression context",
    "Target-disease associations, tractability and safety flags",
    "Target-containing complexes, drug mechanisms and clinical phase",
    "Registered studies found by prespecified gene/drug-context queries"
  ),
  stringsAsFactors = FALSE
)
ct_readme <- data.frame(
  item = c(
    "Purpose", "Candidate set", "Direction rule", "Complex targets",
    "Disease filter", "Trial linkage", "Trial retrieval cap",
    "Evidence-level meaning", "Negative/empty findings", "Reproducibility"
  ),
  description = c(
    paste(
      "Public-database triangulation to prioritize experimental and clinical",
      "follow-up; it does not demonstrate treatment efficacy."
    ),
    paste(ct_config$candidate_genes, collapse = ", "),
    paste(
      "For ALM/GRIP/WALK, positive MR beta implies a hypothesized increase-expression",
      "direction and negative beta a decrease-expression direction. Protein drug",
      "actions are compared only as a hypothesis-generating check."
    ),
    paste(
      "ChEMBL targets containing the candidate UniProt accession are retained,",
      "including protein complexes, because excluding complexes can omit valid",
      "mechanisms such as K-ATP-channel pharmacology for ABCC8."
    ),
    paste(
      "Open Targets associations are limited by a prespecified muscle/metabolic",
      "keyword filter; absence from the sheet is not proof of no association."
    ),
    paste(
      "Only records found through a ChEMBL molecule name are called drug-linked.",
      "Gene-symbol hits do not establish that the intervention targets the gene."
    ),
    paste(
      "At most", ct_config$max_trials_per_query,
      "records are retrieved per query; total_records_for_query identifies truncation."
    ),
    paste(
      "Levels 1-4 are descriptive evidence categories created for this analysis,",
      "not a clinically validated score."
    ),
    paste(
      "Empty or discordant results are retained to prevent selective reporting;",
      "they should be discussed proportionately rather than interpreted as proof",
      "of absence."
    ),
    paste(
      "All queries, status messages, source URLs and retrieval times are retained",
      "in this workbook. Database contents can change after the query date."
    )
  ),
  stringsAsFactors = FALSE
)

if (ncol(ct_ot_associations) == 0L) {
  ct_ot_associations <- ct_empty_table(c(
    "gene_symbol", "gene_ensembl", "disease_id", "disease_name",
    "association_score", "relevance_rule", "source_url"
  ))
}
if (ncol(ct_ot_tractability) == 0L) {
  ct_ot_tractability <- ct_empty_table(c(
    "gene_symbol", "gene_ensembl", "modality", "assessment", "supported", "source"
  ))
}
if (ncol(ct_ot_safety) == 0L) {
  ct_ot_safety <- ct_empty_table(c(
    "gene_symbol", "gene_ensembl", "safety_event", "datasource", "tissue",
    "cell", "literature", "interpretation"
  ))
}
if (ncol(ct_chembl_targets) == 0L) {
  ct_chembl_targets <- ct_empty_table(c(
    "gene_symbol", "uniprot", "target_chembl_id", "target_name", "target_type",
    "complex_context", "source_url"
  ))
}
if (ncol(ct_chembl_mechanisms) == 0L) {
  ct_chembl_mechanisms <- ct_empty_table(c(
    "gene_symbol", "uniprot", "target_chembl_id", "target_name", "target_type",
    "molecule_chembl_id", "molecule_name", "action_type", "mechanism_of_action",
    "direct_interaction", "disease_efficacy", "max_phase", "first_approval",
    "molecule_type", "withdrawn_flag", "source_url"
  ))
}
if (ncol(ct_direction_check) == 0L) {
  ct_direction_check <- ct_empty_table(c(
    "gene_symbol", "molecule_name", "action_type", "mechanism_of_action",
    "genetically_supported_expression_direction", "pharmacology_direction",
    "direction_concordance", "interpretation_limit"
  ))
}
if (ncol(ct_clinical_trials) == 0L) {
  ct_clinical_trials <- ct_empty_table(c(
    "gene_symbol", "matched_term", "match_basis", "target_linkage", "nct_id",
    "brief_title", "overall_status", "study_type", "phase", "enrollment",
    "start_date", "completion_date", "conditions", "interventions",
    "primary_outcomes", "lead_sponsor", "has_results",
    "records_returned_for_query", "total_records_for_query", "retrieval_cap",
    "source_url"
  ))
}

ct_output_tables <- list(
  `00_README` = ct_readme,
  `01_gene_evidence_matrix` = as.data.frame(ct_gene_priority),
  `02_candidate_context` = ct_context,
  `03_HPA_summary` = as.data.frame(ct_hpa_summary),
  `04_HPA_muscle_IHC` = as.data.frame(ct_hpa_muscle),
  `05_OT_associations` = as.data.frame(ct_ot_associations),
  `06_OT_tractability` = as.data.frame(ct_ot_tractability),
  `07_OT_safety` = as.data.frame(ct_ot_safety),
  `08_ChEMBL_targets` = as.data.frame(ct_chembl_targets),
  `09_ChEMBL_mechanisms` = as.data.frame(ct_chembl_mechanisms),
  `10_direction_check` = as.data.frame(ct_direction_check),
  `11_clinical_trials` = as.data.frame(ct_clinical_trials),
  `12_query_log` = as.data.frame(ct_query_log_table),
  `13_sources` = ct_sources
)

ct_workbook <- openxlsx::createWorkbook()
ct_header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
ct_note_style <- openxlsx::createStyle(
  fgFill = "#D9EAF7", valign = "top", wrapText = TRUE
)
ct_link_style <- openxlsx::createStyle(
  fontColour = "#0563C1", textDecoration = "underline"
)
for (sheet_name in names(ct_output_tables)) {
  sheet_data <- ct_output_tables[[sheet_name]]
  openxlsx::addWorksheet(ct_workbook, sheet_name, gridLines = FALSE)
  openxlsx::writeDataTable(
    ct_workbook, sheet_name, sheet_data, tableStyle = "TableStyleMedium2"
  )
  if (ncol(sheet_data) > 0L) {
    openxlsx::addStyle(
      ct_workbook, sheet_name, ct_header_style,
      rows = 1L, cols = seq_len(ncol(sheet_data)), gridExpand = TRUE
    )
    openxlsx::freezePane(ct_workbook, sheet_name, firstRow = TRUE)
    widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- head(as.character(sheet_data[[column_index]]), 500L)
      observed_lengths <- nchar(values)
      observed <- if (length(observed_lengths) == 0L ||
          all(is.na(observed_lengths))) 0L else max(
        observed_lengths, na.rm = TRUE
      )
      min(45, max(11, nchar(names(sheet_data)[column_index]) + 2L, observed + 2L))
    }, numeric(1))
    openxlsx::setColWidths(
      ct_workbook, sheet_name, cols = seq_len(ncol(sheet_data)), widths = widths
    )
    url_columns <- grep("url$|official_url$", names(sheet_data), ignore.case = TRUE)
    if (length(url_columns) > 0L && nrow(sheet_data) > 0L) {
      openxlsx::addStyle(
        ct_workbook, sheet_name, ct_link_style,
        rows = 2:(nrow(sheet_data) + 1L), cols = url_columns,
        gridExpand = TRUE, stack = TRUE
      )
    }
  }
  if (sheet_name == "00_README" && nrow(sheet_data) > 0L) {
    openxlsx::addStyle(
      ct_workbook, sheet_name, ct_note_style,
      rows = 2:(nrow(sheet_data) + 1L), cols = 2L,
      gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setColWidths(ct_workbook, sheet_name, cols = 1:2, widths = c(26, 90))
  }
}
openxlsx::saveWorkbook(
  ct_workbook, ct_config$output_workbook, overwrite = TRUE
)
if (exists("repair_missing_ooxml_relationships", mode = "function")) {
  repair_missing_ooxml_relationships(ct_config$output_workbook)
}
ct_actual_sheets <- openxlsx::getSheetNames(ct_config$output_workbook)
if (!identical(ct_actual_sheets, names(ct_output_tables))) {
  stop("Clinical-translation workbook validation failed: sheet names differ.")
}
if (nrow(ct_hpa_summary) != length(ct_config$candidate_genes) ||
    nrow(ct_gene_priority) != length(ct_config$candidate_genes)) {
  stop("Clinical-translation workbook validation failed: candidate rows missing.")
}

log_step(
  "Public-database clinical translation assessment complete: ",
  ct_config$output_workbook,
  "; ChEMBL mechanisms = ", nrow(ct_chembl_mechanisms),
  "; trial records = ", nrow(ct_clinical_trials),
  "; failed API queries = ", sum(ct_query_log_table$status == "failed")
)


## 20. GTEx v11 skeletal-muscle sQTL extension ------------------------------
## This appended section asks whether the six prespecified Tier A genes also
## show evidence at the splicing layer. GTEx v11 sGenes provide one
## permutation-tested lead LeafCutter phenotype and lead variant per gene.
## Eligible lead sQTLs are tested only against each gene's prespecified primary
## sarcopenia trait. The significant-pairs archive is used for event inventory
## and input validation, not for strict colocalisation: it omits non-significant
## cis-region variants required by coloc.

sq_config <- list(
  sqtl_tar = file.path(config$project_dir, "GTEx_Analysis_v11_sQTL.tar"),
  groups_tar = file.path(
    config$project_dir, "GTEx_Analysis_v11_sQTL_groups.tar"
  ),
  sqtl_members = c(
    sgenes = paste0(
      "GTEx_Analysis_v11_sQTL/",
      "Muscle_Skeletal.v11.sGenes.txt.gz"
    ),
    pairs = paste0(
      "GTEx_Analysis_v11_sQTL/",
      "Muscle_Skeletal.v11.sQTLs.signif_pairs.parquet"
    )
  ),
  groups_member = paste0(
    "GTEx_Analysis_v11_sQTL_phenotype_matrices/",
    "Muscle_Skeletal.v11.sQTL_phenotype_groups.txt.gz"
  ),
  candidate_genes = c("ABCC8", "MAPK1", "RXRA", "SMAD3", "YWHAZ", "ZBTB7B"),
  target_trait = c(
    ABCC8 = "ALM", MAPK1 = "ALM", RXRA = "WALK",
    SMAD3 = "GRIP", YWHAZ = "ALM", ZBTB7B = "ALM"
  ),
  sgene_fdr = 0.05,
  minimum_F = 10,
  maximum_af_difference = 0.20,
  output_workbook = result_path("GTEx_v11_sQTL_TierA_analysis.xlsx")
)

assert_files_exist(
  c(sq_config$sqtl_tar, sq_config$groups_tar),
  "GTEx v11 skeletal-muscle sQTL archive"
)
if (!requireNamespace("arrow", quietly = TRUE)) {
  stop("The GTEx v11 sQTL extension requires the arrow package.")
}

sq_extract_member <- function(tar_file, member, destination) {
  utils::untar(tar_file, files = member, exdir = destination)
  extracted <- file.path(destination, member)
  if (!file.exists(extracted)) {
    stop("Archive member was not extracted: ", member)
  }
  extracted
}

sq_empty_table <- function(columns) {
  out <- as.data.frame(
    setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
  out
}

sq_harmonise_lead <- function(exposure_row, outcome_table, trait) {
  outcome_row <- data.table::as.data.table(outcome_table)[
    SNP == exposure_row$SNP
  ]
  if (nrow(outcome_row) == 0L) {
    return(data.table::data.table(
      gene_symbol = exposure_row$gene_symbol,
      phenotype_id = exposure_row$phenotype_id,
      trait = trait,
      SNP = exposure_row$SNP,
      harmonisation_status = "outcome_variant_unavailable"
    ))
  }
  outcome_row <- outcome_row[1L]

  exposure_effect <- toupper(exposure_row$effect_allele)
  exposure_other <- toupper(exposure_row$other_allele)
  outcome_effect <- toupper(outcome_row$effect_allele)
  outcome_other <- toupper(outcome_row$other_allele)
  exact_match <- exposure_effect == outcome_effect &&
    exposure_other == outcome_other
  swapped_match <- exposure_effect == outcome_other &&
    exposure_other == outcome_effect

  if (!exact_match && !swapped_match) {
    return(data.table::data.table(
      gene_symbol = exposure_row$gene_symbol,
      phenotype_id = exposure_row$phenotype_id,
      trait = trait,
      SNP = exposure_row$SNP,
      exposure_alleles = paste0(exposure_other, "/", exposure_effect),
      outcome_alleles = paste0(outcome_other, "/", outcome_effect),
      harmonisation_status = "allele_pair_mismatch"
    ))
  }

  aligned_beta <- if (exact_match) outcome_row$beta else -outcome_row$beta
  aligned_eaf <- if (exact_match) outcome_row$eaf else 1 - outcome_row$eaf
  af_difference <- abs(exposure_row$eaf - aligned_eaf)
  status <- if (
    is.finite(af_difference) &&
      af_difference <= sq_config$maximum_af_difference
  ) "included" else "allele_frequency_mismatch"

  data.table::data.table(
    gene_symbol = exposure_row$gene_symbol,
    gene_ensembl = exposure_row$gene_ensembl,
    phenotype_id = exposure_row$phenotype_id,
    trait = trait,
    SNP = exposure_row$SNP,
    effect_allele = exposure_effect,
    other_allele = exposure_other,
    exposure_beta = exposure_row$beta,
    exposure_se = exposure_row$se,
    exposure_eaf = exposure_row$eaf,
    outcome_beta_aligned = aligned_beta,
    outcome_se = outcome_row$se,
    outcome_eaf_aligned = aligned_eaf,
    af_difference = af_difference,
    harmonisation_action = if (exact_match) {
      "exact genomic REF/ALT allele match"
    } else {
      "outcome effect reversed to the sQTL ALT allele"
    },
    harmonisation_status = status
  )
}

run_sqtl_extension <- function() {
  extraction_dir <- tempfile("gtex_v11_sqtl_")
  dir.create(extraction_dir, recursive = TRUE)
  on.exit(unlink(extraction_dir, recursive = TRUE, force = TRUE), add = TRUE)

  sgenes_file <- sq_extract_member(
    sq_config$sqtl_tar, sq_config$sqtl_members[["sgenes"]], extraction_dir
  )
  pairs_file <- sq_extract_member(
    sq_config$sqtl_tar, sq_config$sqtl_members[["pairs"]], extraction_dir
  )
  groups_file <- sq_extract_member(
    sq_config$groups_tar, sq_config$groups_member, extraction_dir
  )

  sgenes_all <- data.table::fread(sgenes_file, showProgress = FALSE)
  required_sgene_columns <- c(
    "phenotype_id", "gene_id", "gene_name", "variant_id",
    "rs_id_dbSNP157_GRCh38p14", "ref", "alt", "af",
    "pval_nominal", "slope", "slope_se", "qval",
    "pval_nominal_threshold"
  )
  missing_sgene_columns <- setdiff(required_sgene_columns, names(sgenes_all))
  if (length(missing_sgene_columns) > 0L) {
    stop(
      "GTEx v11 sGenes columns missing: ",
      paste(missing_sgene_columns, collapse = ", ")
    )
  }

  candidate_sgenes <- sgenes_all[gene_name %in% sq_config$candidate_genes]
  if (data.table::uniqueN(candidate_sgenes$gene_name) !=
      length(sq_config$candidate_genes)) {
    stop("Not all six Tier A genes were found in the Muscle_Skeletal sGenes file.")
  }
  candidate_sgenes[, `:=`(
    gene_symbol = gene_name,
    gene_ensembl = strip_ensembl_version(gene_id),
    SNP = as.character(rs_id_dbSNP157_GRCh38p14),
    effect_allele = toupper(as.character(alt)),
    other_allele = toupper(as.character(ref)),
    beta = as.numeric(slope),
    se = as.numeric(slope_se),
    eaf = as.numeric(af),
    F_stat = (as.numeric(slope) / as.numeric(slope_se))^2,
    target_trait = unname(sq_config$target_trait[gene_name]),
    sgene_fdr_significant = as.numeric(qval) <= sq_config$sgene_fdr,
    lead_pair_empirically_significant =
      as.numeric(pval_nominal) <= as.numeric(pval_nominal_threshold)
  )]

  candidate_gene_ids <- unique(candidate_sgenes$gene_id)
  pairs_dataset <- arrow::open_dataset(pairs_file, format = "parquet")
  candidate_pairs <- pairs_dataset |>
    dplyr::filter(group_id %in% candidate_gene_ids) |>
    dplyr::collect() |>
    data.table::as.data.table()
  candidate_pairs <- merge(
    candidate_pairs,
    unique(candidate_sgenes[, .(group_id = gene_id, gene_symbol)]),
    by = "group_id", all.x = TRUE, sort = FALSE
  )
  candidate_pairs[, `:=`(
    F_stat = (as.numeric(slope) / as.numeric(slope_se))^2,
    conventional_p_lt_5e_8 = as.numeric(pval_nominal) < 5e-8,
    empirical_pair_significant =
      as.numeric(pval_nominal) <= as.numeric(pval_nominal_threshold)
  )]

  phenotype_groups <- data.table::fread(
    groups_file, header = FALSE,
    col.names = c("phenotype_id", "group_gene_id"),
    showProgress = FALSE
  )
  candidate_group_map <- phenotype_groups[
    group_gene_id %in% candidate_gene_ids
  ]
  pair_map_check <- merge(
    unique(candidate_pairs[, .(phenotype_id, group_id)]),
    candidate_group_map,
    by = "phenotype_id", all.x = TRUE, sort = FALSE
  )
  pair_map_check[, mapping_consistent := group_id == group_gene_id]

  eligible_leads <- candidate_sgenes[
    sgene_fdr_significant & lead_pair_empirically_significant &
      is.finite(F_stat) & F_stat >= sq_config$minimum_F &
      grepl("^rs[0-9]+$", SNP)
  ]

  harmonisation_rows <- lapply(seq_len(nrow(eligible_leads)), function(i) {
    exposure_row <- eligible_leads[i]
    trait <- exposure_row$target_trait
    sq_harmonise_lead(exposure_row, gwas_list[[trait]], trait)
  })
  harmonisation <- bind_nonempty(harmonisation_rows)

  mr_results <- harmonisation[harmonisation_status == "included", .(
    gene_symbol,
    gene_ensembl,
    phenotype_id,
    trait,
    SNP,
    method = "Wald ratio",
    nsnp = 1L,
    beta = outcome_beta_aligned / exposure_beta,
    se = outcome_se / abs(exposure_beta),
    exposure_beta,
    exposure_se,
    outcome_beta = outcome_beta_aligned,
    outcome_se,
    F_stat = (exposure_beta / exposure_se)^2,
    exposure_effect_allele = effect_allele,
    exposure_other_allele = other_allele,
    exposure_eaf,
    outcome_eaf = outcome_eaf_aligned,
    af_difference
  )]
  if (nrow(mr_results) > 0L) {
    mr_results[, `:=`(
      ci_lower = beta - stats::qnorm(0.975) * se,
      ci_upper = beta + stats::qnorm(0.975) * se,
      pval = 2 * stats::pnorm(-abs(beta / se))
    )]
    mr_results[, fdr_targeted := stats::p.adjust(pval, method = "BH")]
    mr_results[, splicing_direction := data.table::fifelse(
      beta > 0,
      "higher normalized intron-excision phenotype associated with higher trait",
      "higher normalized intron-excision phenotype associated with lower trait"
    )]
    mr_results[, interpretation_limit := paste0(
      "LeafCutter normalized intron-excision phenotype; sign does not by itself ",
      "identify exon inclusion, transcript isoform, or protein abundance."
    )]
  }

  event_inventory <- candidate_pairs[, .(
    n_reported_significant_pairs = .N,
    n_reported_splice_events = data.table::uniqueN(phenotype_id),
    minimum_nominal_p = min(pval_nominal, na.rm = TRUE),
    n_pairs_p_lt_5e_8 = sum(conventional_p_lt_5e_8, na.rm = TRUE)
  ), by = gene_symbol]

  layer_summary <- data.table::data.table(
    gene_symbol = sq_config$candidate_genes,
    target_trait = unname(sq_config$target_trait[sq_config$candidate_genes])
  )
  layer_summary <- merge(
    layer_summary,
    candidate_sgenes[, .(
      gene_symbol, gene_ensembl, lead_phenotype_id = phenotype_id,
      lead_sQTL = SNP, sgene_qval = qval, lead_sQTL_pval = pval_nominal,
      lead_sQTL_F = F_stat, sgene_fdr_significant
    )],
    by = "gene_symbol", all.x = TRUE, sort = FALSE
  )
  layer_summary <- merge(
    layer_summary, event_inventory,
    by = "gene_symbol", all.x = TRUE, sort = FALSE
  )
  if (nrow(mr_results) > 0L) {
    layer_summary <- merge(
      layer_summary,
      mr_results[, .(
        gene_symbol, sqtl_mr_beta = beta, sqtl_mr_se = se,
        sqtl_mr_pval = pval, sqtl_mr_fdr = fdr_targeted
      )],
      by = "gene_symbol", all.x = TRUE, sort = FALSE
    )
  } else {
    layer_summary[, `:=`(
      sqtl_mr_beta = NA_real_, sqtl_mr_se = NA_real_,
      sqtl_mr_pval = NA_real_, sqtl_mr_fdr = NA_real_
    )]
  }
  for (column_name in c(
    "n_reported_significant_pairs", "n_reported_splice_events",
    "n_pairs_p_lt_5e_8"
  )) {
    layer_summary[is.na(get(column_name)), (column_name) := 0L]
  }
  layer_summary[, expression_layer :=
    "supported by the primary skeletal-muscle eQTL-MR/coloc Tier A analysis"]
  layer_summary[, splicing_layer := data.table::fcase(
    is.finite(sqtl_mr_fdr) & sqtl_mr_fdr < 0.05,
    "targeted sQTL-MR supported after FDR correction",
    sgene_fdr_significant & is.finite(sqtl_mr_pval),
    "skeletal-muscle sGene detected; targeted sQTL-MR not FDR-supported",
    sgene_fdr_significant,
    "skeletal-muscle sGene detected; targeted MR unavailable",
    default = "no skeletal-muscle sGene evidence at q <= 0.05"
  )]
  layer_summary[, integrated_inference := data.table::fifelse(
    is.finite(sqtl_mr_fdr) & sqtl_mr_fdr < 0.05,
    "expression and splicing layers both implicated",
    "expression layer implicated; independent splicing support not established"
  )]
  layer_summary[, sqtl_supported_order := as.integer(
    is.finite(sqtl_mr_fdr) & sqtl_mr_fdr < 0.05
  )]
  data.table::setorder(
    layer_summary, -sqtl_supported_order, sqtl_mr_fdr,
    gene_symbol, na.last = TRUE
  )
  layer_summary[, sqtl_supported_order := NULL]

  input_qc <- data.table::data.table(
    check = c(
      "candidate_genes_found_in_sGenes",
      "candidate_significant_pair_rows",
      "candidate_splice_events_in_significant_pairs",
      "pair_to_group_mapping_coverage",
      "pair_to_group_mapping_consistency",
      "eligible_gene_level_lead_sQTLs",
      "targeted_MR_tests_completed"
    ),
    observed = c(
      data.table::uniqueN(candidate_sgenes$gene_symbol),
      nrow(candidate_pairs),
      data.table::uniqueN(candidate_pairs$phenotype_id),
      mean(!is.na(pair_map_check$group_gene_id)),
      mean(pair_map_check$mapping_consistent, na.rm = TRUE),
      nrow(eligible_leads),
      nrow(mr_results)
    ),
    expected_or_rule = c(
      length(sq_config$candidate_genes),
      ">= 0", ">= 0", "1.0", "1.0",
      "q <= 0.05, lead p <= empirical threshold, F >= 10, valid rsID",
      "one prespecified gene-trait test per eligible lead sQTL"
    )
  )

  coloc_readiness <- data.table::data.table(
    gene_symbol = sq_config$candidate_genes,
    strict_sQTL_coloc_run = FALSE,
    available_local_source =
      "GTEx v11 sGenes plus sQTLs.signif_pairs only",
    reason_not_run = paste0(
      "Strict coloc requires dense cis-region summary statistics including ",
      "non-significant variants; a significant-pairs-only file is selected data."
    ),
    required_additional_data = paste0(
      "Complete variant-level cis statistics for the same skeletal-muscle ",
      "LeafCutter phenotype and the matched outcome locus."
    )
  )

  list(
    candidate_sgenes = candidate_sgenes,
    candidate_pairs = candidate_pairs,
    eligible_leads = eligible_leads,
    harmonisation = harmonisation,
    mr_results = mr_results,
    layer_summary = layer_summary,
    input_qc = input_qc,
    coloc_readiness = coloc_readiness
  )
}

log_step("Running GTEx v11 skeletal-muscle sQTL extension")
sq_results <- run_sqtl_extension()

sq_primary_finding <- if (
  nrow(sq_results$mr_results) > 0L &&
    any(sq_results$mr_results$fdr_targeted < 0.05, na.rm = TRUE)
) {
  supported <- sq_results$mr_results[
    fdr_targeted < 0.05,
    paste0(gene_symbol, "-", trait, " (FDR=", signif(fdr_targeted, 3), ")")
  ]
  paste("Targeted splicing-layer support:", paste(supported, collapse = "; "))
} else {
  "No targeted sQTL-MR result passed FDR correction."
}

sq_readme <- data.frame(
  item = c(
    "Purpose", "Candidate set", "Tissue and QTL release", "Exposure unit",
    "Instrument eligibility", "Outcome scope", "Harmonisation",
    "Multiple testing", "Primary finding", "Colocalisation boundary",
    "Protein-layer boundary", "Interpretation"
  ),
  value = c(
    "Assess whether the six primary Tier A genes also have skeletal-muscle splicing-layer evidence.",
    paste(sq_config$candidate_genes, collapse = ", "),
    "GTEx v11 Muscle_Skeletal sQTL; LeafCutter intron-excision phenotypes.",
    "One GTEx permutation-tested lead sQTL and lead splice phenotype per sGene; Wald-ratio MR.",
    "sGene q <= 0.05, lead nominal p <= GTEx empirical threshold, F >= 10, and valid dbSNP rsID.",
    "One prespecified primary trait per gene: ALM, GRIP, or WALK as defined in the main analysis.",
    "Matched by rsID; genomic REF/ALT allele pairs had to match exactly or by allele swap; allele-frequency difference <= 0.20.",
    "Benjamini-Hochberg FDR across the completed targeted sQTL-MR tests.",
    sq_primary_finding,
    "Not run from the significant-pairs archive because selection to significant variants violates dense-region coloc input requirements.",
    "Skeletal-muscle pQTL was not assessed because no pQTL summary-statistic source was supplied.",
    "A LeafCutter effect concerns normalized intron excision; it does not alone identify exon inclusion, a transcript isoform, protein abundance, or a drug direction."
  ),
  stringsAsFactors = FALSE
)

sq_sources <- data.frame(
  source = c("GTEx v11 adult QTL downloads", "LeafCutter method"),
  use = c(
    "Muscle_Skeletal sGenes, significant sQTL pairs, and phenotype-group mapping",
    "Definition and interpretation of intron-excision phenotypes"
  ),
  official_url = c(
    "https://gtexportal.org/home/downloads/adult-gtex/overview",
    "https://pubmed.ncbi.nlm.nih.gov/29229983/"
  ),
  stringsAsFactors = FALSE
)

sq_output_tables <- list(
  `00_README` = sq_readme,
  `01_input_QC` = as.data.frame(sq_results$input_qc),
  `02_candidate_sGenes` = as.data.frame(sq_results$candidate_sgenes),
  `03_significant_sQTL_pairs` = as.data.frame(sq_results$candidate_pairs),
  `04_eligible_lead_sQTL` = as.data.frame(sq_results$eligible_leads),
  `05_harmonisation` = as.data.frame(sq_results$harmonisation),
  `06_targeted_sQTL_MR` = if (nrow(sq_results$mr_results) > 0L) {
    as.data.frame(sq_results$mr_results)
  } else {
    sq_empty_table(c(
      "gene_symbol", "phenotype_id", "trait", "SNP", "method",
      "nsnp", "beta", "se", "pval", "fdr_targeted"
    ))
  },
  `07_layer_summary` = as.data.frame(sq_results$layer_summary),
  `08_coloc_readiness` = as.data.frame(sq_results$coloc_readiness),
  `09_sources` = sq_sources
)

sq_workbook <- openxlsx::createWorkbook()
sq_header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
sq_note_style <- openxlsx::createStyle(
  fgFill = "#D9EAF7", valign = "top", wrapText = TRUE
)
for (sheet_name in names(sq_output_tables)) {
  sheet_data <- sq_output_tables[[sheet_name]]
  openxlsx::addWorksheet(sq_workbook, sheet_name, gridLines = FALSE)
  openxlsx::writeDataTable(
    sq_workbook, sheet_name, sheet_data, tableStyle = "TableStyleMedium2"
  )
  if (ncol(sheet_data) > 0L) {
    openxlsx::addStyle(
      sq_workbook, sheet_name, sq_header_style,
      rows = 1L, cols = seq_len(ncol(sheet_data)), gridExpand = TRUE
    )
    openxlsx::freezePane(sq_workbook, sheet_name, firstRow = TRUE)
    widths <- vapply(seq_len(ncol(sheet_data)), function(column_index) {
      values <- head(as.character(sheet_data[[column_index]]), 500L)
      observed_lengths <- nchar(values)
      observed <- if (
        length(observed_lengths) == 0L || all(is.na(observed_lengths))
      ) 0L else max(observed_lengths, na.rm = TRUE)
      min(50, max(11, nchar(names(sheet_data)[column_index]) + 2L, observed + 2L))
    }, numeric(1))
    openxlsx::setColWidths(
      sq_workbook, sheet_name, cols = seq_len(ncol(sheet_data)), widths = widths
    )
  }
  if (sheet_name == "00_README" && nrow(sheet_data) > 0L) {
    openxlsx::addStyle(
      sq_workbook, sheet_name, sq_note_style,
      rows = 2:(nrow(sheet_data) + 1L), cols = 2L,
      gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setColWidths(sq_workbook, sheet_name, cols = 1:2, widths = c(26, 95))
  }
}
openxlsx::saveWorkbook(
  sq_workbook, sq_config$output_workbook, overwrite = TRUE
)
if (exists("repair_missing_ooxml_relationships", mode = "function")) {
  repair_missing_ooxml_relationships(sq_config$output_workbook)
}
sq_actual_sheets <- openxlsx::getSheetNames(sq_config$output_workbook)
if (!identical(sq_actual_sheets, names(sq_output_tables))) {
  stop("GTEx v11 sQTL workbook validation failed: sheet names differ.")
}
if (nrow(sq_results$layer_summary) != length(sq_config$candidate_genes)) {
  stop("GTEx v11 sQTL workbook validation failed: candidate rows missing.")
}

log_step(
  "GTEx v11 skeletal-muscle sQTL extension complete: ",
  sq_config$output_workbook,
  "; eligible lead sQTLs = ", nrow(sq_results$eligible_leads),
  "; targeted MR tests = ", nrow(sq_results$mr_results),
  "; FDR-supported tests = ",
  sum(sq_results$mr_results$fdr_targeted < 0.05, na.rm = TRUE)
)









####--------------------------------MAGMA-------------------------####

## MAGMA 01：准备输入文件
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

magma_dir <- file.path(config$paths$results, "MAGMA")
dir.create(magma_dir, recursive = TRUE, showWarnings = FALSE)


## ------------1. 冻结候选基因集-------------------
## 必须确认这里仍然是论文使用的771个基因

gene_membership <- as.data.table(gene_set$membership)
pathway_membership <- as.data.table(gene_set$pathways)

n_candidate <- uniqueN(gene_membership$gene_symbol)

if (n_candidate != 771L) {
  stop(
    "当前MSigDB产生了", n_candidate,
    "个候选基因，而不是论文中的771个。请使用原始分析保存的基因集，",
    "不要继续运行MAGMA。")}


## -------------2. SYMBOL转换为MAGMA基因位置文件使用的Entrez ID----------------------

symbol_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gene_membership$gene_symbol),
  keytype = "SYMBOL",
  columns = "ENTREZID")

symbol_entrez <- as.data.table(symbol_entrez)
symbol_entrez <- symbol_entrez[
  !is.na(SYMBOL) & !is.na(ENTREZID),
  .(gene_symbol = SYMBOL, ENTREZID = as.character(ENTREZID))
]
symbol_entrez <- unique(symbol_entrez)


##---------------------------- 3. 完整771基因集----------------------------


full_set <- merge(
  gene_membership[, .(gene_symbol)],
  symbol_entrez,
  by = "gene_symbol"
)[
  ,
  .(SET = "IR_FULL_771", ENTREZID)
]


##---------------------------- 4. 核心通路并集----------------------------

canonical_set <- merge(
  gene_membership[in_canonical == TRUE, .(gene_symbol)],
  symbol_entrez,
  by = "gene_symbol"
)[
  ,
  .(SET = "IR_CANONICAL_UNION", ENTREZID)
]


##---------------------------- 5. 各预设IR通路----------------------------

pathway_sets <- merge(
  pathway_membership[, .(SET = gs_name, gene_symbol)],
  symbol_entrez,
  by = "gene_symbol"
)[
  ,
  .(SET, ENTREZID)
]

magma_sets <- unique(
  rbindlist(
    list(full_set, canonical_set, pathway_sets),
    use.names = TRUE))

## MAGMA列式格式：第1列为基因集，第2列为Entrez基因ID
set_file <- file.path(magma_dir, "IR_gene_sets_entrez.txt")

fwrite(
  magma_sets,
  file = set_file,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE)

## 保存映射质控
mapping_qc <- data.table(
  metric = c(
    "original_IR_symbols",
    "mapped_IR_symbols",
    "mapped_IR_entrez_ids",
    "symbol_mapping_rate"
  ),
  value = c(
    uniqueN(gene_membership$gene_symbol),
    uniqueN(symbol_entrez$gene_symbol),
    uniqueN(symbol_entrez$ENTREZID),
    uniqueN(symbol_entrez$gene_symbol) /
      uniqueN(gene_membership$gene_symbol)))

fwrite(
  mapping_qc,
  file.path(magma_dir, "IR_gene_mapping_QC.tsv"),
  sep = "\t")


##------------------------- 6. 生成三个结局的SNP-P文件-------------------------
target_traits <- c("ALM", "GRIP", "WALK")

# [NEW: MAGMA reproducibility] The duplicate-rsID resolver uses the same
# chromosome, position and allele fields as the LD reference used by MAGMA.
ref_prefix <- file.path(config$project_dir, "1000G_EUR/1000G.EUR.QC")
bim_files <- paste0(ref_prefix, ".", 1:22, ".bim")
assert_files_exist(bim_files, "1000G EUR per-chromosome BIM")
magma_reference <- magma_read_bim_reference(bim_files)

magma_resolution_audits <- list()
magma_resolution_qc <- list()

# [NEW: MAGMA reproducibility] These counts are part of the frozen source-data
# audit. A changed count stops the run instead of silently changing the reported
# MAGMA universe or cross-phenotype comparisons.
magma_expected_resolution <- data.table(
  trait = target_traits,
  duplicate_rsids = c(26228L, 14712L, 14712L),
  replaced_first_records = c(2861L, 2562L, 2562L),
  excluded_duplicate_rsids = c(1L, 0L, 0L),
  pval_rows_written = c(17914304L, 9834830L, 9834830L)
)

for (trait in target_traits) {
  # Re-read the valid raw outcome rows so MAGMA does not inherit the
  # first-valid rsID rule used by the MR/eQTL branch above.
  gwas_dt <- vcf_to_gwas(
    config$files$outcomes[[trait]],
    trait,
    config$outcome_n[[trait]],
    deduplicate = FALSE,
    audit_scope = "MAGMA_raw"
  )
  resolved <- resolve_magma_duplicate_rsids(
    gwas_dt, magma_reference, trait = trait
  )
  magma_resolution_audits[[trait]] <- resolved$audit

  magma_pval <- resolved$rows[
    grepl("^rs", SNP) &
      is.finite(pval) &
      pval >= 0 &
      pval <= 1,
    .(
      SNP = as.character(SNP),
      P = pmax(as.numeric(pval), 1e-300)
    )
  ]
  if (anyDuplicated(magma_pval$SNP)) {
    stop("MAGMA duplicate-rsID resolution failed for ", trait)
  }
  if (any(!is.finite(magma_pval$P))) {
    stop("MAGMA P-value input contains non-finite values for ", trait)
  }

  magma_resolution_qc[[trait]] <- data.table(
    trait = trait,
    raw_valid_rows = nrow(gwas_dt),
    duplicate_rsids = nrow(resolved$audit),
    retained_rows = nrow(resolved$rows),
    excluded_duplicate_rsids = sum(
      grepl("^excluded_", resolved$audit$status)
    ),
    replaced_first_records = sum(
      resolved$audit$status == "replaced_first_record"
    ),
    absent_from_ld_reference = sum(
      resolved$audit$status == "absent_from_LD_reference"
    ),
    pval_rows_written = nrow(magma_pval)
  )

  fwrite(
    magma_pval,
    file.path(magma_dir, paste0(trait, "_MAGMA_pval.tsv")),
    sep = "\t",
    quote = FALSE
  )
}

magma_resolution_audit <- data.table::rbindlist(
  magma_resolution_audits, use.names = TRUE, fill = TRUE
)
magma_resolution_qc <- data.table::rbindlist(
  magma_resolution_qc, use.names = TRUE, fill = TRUE
)
qc_check <- merge(
  magma_resolution_qc[, .(trait, duplicate_rsids,
                          replaced_first_records,
                          excluded_duplicate_rsids,
                          pval_rows_written)],
  magma_expected_resolution,
  by = "trait",
  suffixes = c("", "_expected"),
  all = TRUE
)
qc_fields <- c("duplicate_rsids", "replaced_first_records",
               "excluded_duplicate_rsids", "pval_rows_written")
if (any(vapply(qc_fields, function(field) {
  !identical(qc_check[[field]], qc_check[[paste0(field, "_expected")]])
}, logical(1)))) {
  stop("MAGMA duplicate-rsID audit counts differ from the frozen source-data audit.")
}
openxlsx::write.xlsx(
  list(variant_resolution = as.data.frame(magma_resolution_audit)),
  file = result_path("MAGMA_variant_resolution_audit.xlsx"),
  overwrite = TRUE,
  asTable = TRUE
)
data.table::fwrite(
  magma_resolution_qc,
  result_path("MAGMA_variant_resolution_qc.tsv"),
  sep = "\t",
  quote = FALSE
)




magma_exe <- file.path(config$project_dir, "magma.exe")

magma_path <- function(path, mustWork = TRUE) {
  normal_path <- normalizePath(path, winslash = "\\", mustWork = mustWork)
  path_dir <- dirname(normal_path)
  file.path(utils::shortPathName(path_dir), basename(normal_path))}

ref_chr_prefix <- file.path(
  config$project_dir,
  "1000G_EUR",
  "1000G.EUR.QC.#CHR#")

gene_loc <- file.path(config$project_dir, "NCBI37.3.gene.loc")

annot_prefix <- file.path(magma_dir, "GRCh37_gene_body")
snp_loc <- file.path(magma_dir, "MAGMA_snp_loc.tsv")

snp_loc_tables <- lapply(bim_files, function(bim_file) {
  bim <- data.table::fread(
    bim_file,
    header = FALSE,
    select = c(1L, 2L, 4L),
    col.names = c("CHR", "SNP", "BP"),
    showProgress = FALSE
  )
  bim[, .(SNP, CHR, BP)]
})
data.table::fwrite(
  data.table::rbindlist(snp_loc_tables, use.names = TRUE),
  snp_loc,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE
)

status <- system2(
  magma_exe,
  args = c(
    "--annotate",
    "--snp-loc", magma_path(snp_loc),
    "--gene-loc", magma_path(gene_loc),
    "--out", magma_path(annot_prefix, mustWork = FALSE)
  )
)

if (status != 0L) {
  stop("MAGMA SNP–gene annotation failed.")
}





trait_n <- c(ALM  = 450243,GRIP = 461089,WALK = 459915)

set_file <- file.path(magma_dir, "IR_gene_sets_entrez.txt")
gene_annot <- paste0(annot_prefix, ".genes.annot.txt")

for (trait in names(trait_n)) {
  
  pval_file <- file.path(
    magma_dir,
    paste0(trait, "_MAGMA_pval.tsv")
  )
  
  gene_prefix <- file.path(
    magma_dir,
    paste0(trait, "_gene")
  )
  
  gsa_prefix <- file.path(
    magma_dir,
    paste0(trait, "_IR_competitive")
  )
  
  ## 基因层分析
  for (chromosome in 1:22) {
    status_gene <- system2(
      magma_exe,
      args = c(
        "--bfile", magma_path(ref_chr_prefix, mustWork = FALSE),
        "--pval", magma_path(pval_file),
        paste0("N=", trait_n[[trait]]),
        "--gene-annot", magma_path(gene_annot),
        "--batch", chromosome, "chr",
        "--out", magma_path(gene_prefix, mustWork = FALSE)
      )
    )

    if (status_gene != 0L) {
      stop("MAGMA gene analysis failed: ", trait, ", chromosome ", chromosome)
    }
  }

  status_merge <- system2(
    magma_exe,
    args = c(
      "--merge", magma_path(gene_prefix, mustWork = FALSE),
      "--out", magma_path(gene_prefix, mustWork = FALSE)
    )
  )
  
  if (status_merge != 0L) {
    stop("MAGMA gene batch merge failed: ", trait)
  }
  
  ## 竞争性基因集分析
  status_set <- system2(
    magma_exe,
    args = c(
      "--gene-results", magma_path(paste0(gene_prefix, ".genes.raw")),
      "--set-annot", magma_path(set_file), "col=2,1",
      "--settings", "gene-info",
      "--out", magma_path(gsa_prefix, mustWork = FALSE)
    )
  )
  
  if (status_set != 0L) {
    stop("MAGMA gene-set analysis failed: ", trait)
  }
}



read_magma_gsa <- function(trait) {
  
  file <- file.path(
    magma_dir,
    paste0(trait, "_IR_competitive.gsa.out.txt")
  )
  
  x <- fread(file, skip = "VARIABLE", header = TRUE)
  x[, trait := trait]
  x[]
}

magma_results <- rbindlist(
  lapply(c("ALM", "GRIP", "WALK"), read_magma_gsa),
  fill = TRUE
)

## 主检验和次要通路检验分开校正
magma_results[
  ,
  analysis_family := fifelse(
    VARIABLE == "IR_FULL_771",
    "primary_IR_total",
    "secondary_pathways"
  )
]

magma_results[
  ,
  FDR := p.adjust(P, method = "BH"),
  by = analysis_family
]

setcolorder(
  magma_results,
  c(
    "trait", "VARIABLE", "analysis_family",
    "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR"
  )
)

fwrite(
  magma_results,
  file.path(magma_dir, "MAGMA_competitive_results.tsv"),
  sep = "\t")


magma_results[ order(analysis_family, FDR, P)]
 










# [RESUME HERE AFTER MAIN ANALYSIS]
# If the time-consuming upstream analyses have already produced the
# results_public_release workbooks, rerun from this marker to regenerate only
# manuscript figures, tables and catalog outputs. In a fresh R session, make
# sure config has already been created or use the standalone manuscript-output
# source script with the same sheet-name fixes.
# [NEW: integrated manuscript outputs and formal MAGMA revision analyses]
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

project_dir <- config$project_dir
result_dir <- file.path(project_dir, "results_public_release")
output_dir <- file.path(project_dir, "manuscript_outputs_public_release")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  main = file.path(result_dir, "manuscript_analysis_tables.xlsx"),
  magic = file.path(result_dir, "MAGIC_IR_mechanistic_extension.xlsx"),
  snrna = file.path(result_dir, "GSE167186_snRNA_celltype_validation.xlsx"),
  clinical = file.path(result_dir, "clinical_translation_public_databases.xlsx"),
  sqtl = file.path(result_dir, "GTEx_v11_sQTL_TierA_analysis.xlsx")
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
gefos_mr_support <- read_sheet("main", "09_GEFOS_external_evidence")
gefos_coloc_support <- read_sheet("main", "09_GEFOS_coloc_all_priors")
gefos_coloc_summary <- gefos_coloc_support[, .(
  gefos_coloc_max_PP4 = max(PP.H4, na.rm = TRUE)
), by = gene_symbol]

evidence_base <- data.table(gene_symbol = candidate_order)
evidence_base <- Reduce(
  function(x, y) merge(x, y, by = "gene_symbol", all.x = TRUE, sort = FALSE),
  list(
    evidence_base,
    fusion[, .(gene_symbol, fusion_replicated = replicated)],
    susie_summary,
    sn_support,
    sqtl_support,
    gefos_mr_support[, .(
      gene_symbol,
      gefos_test_status = test_status,
      gefos_external_fdr = external_fdr_4gene,
      gefos_bonferroni_supported = bonferroni_directional_support
    )],
    gefos_coloc_summary,
    clinical[, .(gene_symbol, hpa_skeletal_muscle_protein_reported)]
  )
)
evidence_base[is.na(fusion_replicated), fusion_replicated := FALSE]
evidence_base[is.na(snrna_supported), snrna_supported := FALSE]
evidence_base[is.na(sqtl_supported), sqtl_supported := FALSE]
evidence_base[is.na(hpa_skeletal_muscle_protein_reported),
              hpa_skeletal_muscle_protein_reported := FALSE]
evidence_base[is.na(gefos_test_status), gefos_test_status := "not_applicable"]

make_evidence_row <- function(gene_symbol, evidence, group, classification) {
  data.table(
    gene_symbol = gene_symbol,
    evidence = evidence,
    evidence_group = group,
    classification = classification
  )
}

evidence_long <- rbindlist(list(
  evidence_base[, make_evidence_row(
    gene_symbol, "Discovery cis-MR\n+ coloc", "Discovery",
    "Statistical criterion met"
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "SuSiE\nshared signal", "Discovery",
    ifelse(is.finite(susie_max_PP4) & susie_max_PP4 >= 0.80,
           "Statistical criterion met", "Criterion not met")
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "FUSION eQTL-source\nrobustness", "Same-outcome\nrobustness",
    ifelse(fusion_replicated, "Statistical criterion met", "Criterion not met")
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "GEFOS MR\nassociation", "External outcome follow-up",
    ifelse(
      gene_symbol %in% gefos_mr_support$gene_symbol,
      ifelse(
        gefos_test_status == "not_tested", "Not evaluable",
        ifelse(gefos_bonferroni_supported %in% TRUE |
                 (is.finite(gefos_external_fdr) & gefos_external_fdr < 0.05),
               "Statistical criterion met", "Criterion not met")
      ),
      "Not applicable"
    )
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "GEFOS outcome\ncolocalization", "External outcome follow-up",
    ifelse(
      gene_symbol %in% gefos_mr_support$gene_symbol,
      ifelse(is.finite(gefos_coloc_max_PP4) & gefos_coloc_max_PP4 >= 0.80,
             "Statistical criterion met", "Criterion not met"),
      "Not applicable"
    )
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "snRNA age\nassociation", "Descriptive context",
    ifelse(snrna_supported, "Statistical criterion met", "Criterion not met")
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "sQTL-MR\nclue", "Descriptive context",
    ifelse(sqtl_supported, "Context only", "Criterion not met")
  )],
  evidence_base[, make_evidence_row(
    gene_symbol, "HPA protein\nrecord", "Descriptive context",
    ifelse(hpa_skeletal_muscle_protein_reported, "Context only", "Criterion not met")
  )]
))

evidence_levels <- c(
  "Statistical criterion met", "Context only", "Criterion not met",
  "Not evaluable", "Not applicable"
)
evidence_labels <- c(
  "Statistical criterion met" = "S",
  "Context only" = "C",
  "Criterion not met" = "",
  "Not evaluable" = "NE",
  "Not applicable" = "NA"
)
evidence_long[, `:=`(
  gene_symbol = factor(gene_symbol, levels = rev(candidate_order)),
  evidence = factor(evidence, levels = c(
    "Discovery cis-MR\n+ coloc", "SuSiE\nshared signal",
    "FUSION eQTL-source\nrobustness",
    "GEFOS MR\nassociation", "GEFOS outcome\ncolocalization",
    "snRNA age\nassociation", "sQTL-MR\nclue", "HPA protein\nrecord"
  )),
  evidence_group = factor(evidence_group, levels = c(
    "Discovery", "Same-outcome\nrobustness",
    "External outcome follow-up", "Descriptive context"
  )),
  classification = factor(classification, levels = evidence_levels),
  tile_label = unname(evidence_labels[as.character(classification)])
)]

figure4 <- ggplot(
  evidence_long,
  aes(x = evidence, y = gene_symbol, fill = classification)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = tile_label), size = 3.5, fontface = "bold") +
  facet_grid(
    . ~ evidence_group, scales = "free_x", space = "free_x", switch = "x"
  ) +
  scale_fill_manual(values = c(
    "Statistical criterion met" = "#4A86B8",
    "Context only" = "#E8C36A",
    "Criterion not met" = "#ECEFF2",
    "Not evaluable" = "#B8BEC6",
    "Not applicable" = "#F7F7F7"
  ), drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "Evidence classification") +
  theme_manuscript(10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background.x = element_rect(fill = "#F0F2F5", colour = "white"),
    strip.text.x = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9.5),
    legend.key.width = grid::unit(0.32, "in"),
    legend.key.height = grid::unit(0.22, "in"),
    panel.spacing.x = grid::unit(0.16, "in"),
    plot.margin = margin(10, 10, 10, 10)
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
    "24 archived IR-related pathways", "Human skeletal muscle",
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
      # [REVISED] Draw row indices, never Entrez identifiers.
      pool <- setdiff(pool_table$gene_index, selected)
    }
    chosen <- pool[sample.int(length(pool), needed, replace = length(pool) < needed)]
    selected <- c(selected, chosen)
  }

  selected <- unique(selected)
  if (length(selected) < n_target) {
    fill_pool <- setdiff(pool_table$gene_index, selected)
    selected <- c(
      selected,
      fill_pool[sample.int(length(fill_pool), n_target - length(selected), replace = FALSE)]
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
  # [REVISED] Membership is embedded in this public script; no archive dependency.
  candidate_archive_file <- file.path(project_dir, "IR_sarcopenia_public_release.R")
  candidate_archive_display <- "Embedded archived membership in public script"

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
      source_pathways = "all archived pathways",
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
      "Remove every gene annotated to one pathway from IR_FULL_771, including genes shared with retained pathways.",
      paste0(
        config$random_iterations,
        " random gene sets initially matched by chromosome, NSNP quintile and ",
        "gene-length quintile; inherited sparse-stratum handling relaxes chromosome ",
        "then uses the remaining pool. Thus matching is approximate."
      ),
      paste(
        "This rerun reads a fixed archived membership without outcome-based selection.",
        "File timestamps do not establish that the historical candidate set was",
        "frozen before every outcome analysis; no such chronological claim is made."
      ),
      paste0(
        "Public script mtime (not historical freeze evidence): ",
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
        "with MAGMA beta coefficients. BH correction uses five alternative set definitions ",
        "and a separate family of 24 pathway-deletion tests."
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
  `S14_GEFOS_ALM` = read_sheet("main", "09_GEFOS_external_evidence"),
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
    "Archived pathways and membership",
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


# [NEW: Reviewer 5/7] Deliver all-prior and non-interpretable results alongside
# the selected signal summary; keep the original S1-S18 numbering unchanged.
additional_evidence <- list(
  E1_SuSiE_QC = read_sheet("main", "04c_coloc_SuSiE_TierA_QC"),
  E2_SuSiE_all_priors = read_sheet("main", "04d_coloc_SuSiE_all_priors"),
  E3_SuSiE_credible_sets = read_sheet("main", "04e_coloc_SuSiE_credible_sets"),
  E4_SuSiE_CS_members = read_sheet("main", "04f_coloc_SuSiE_CS_members"),
  E5_TierA_prior_summary = read_sheet("main", "04g_TierA_coloc_prior_summary"),
  E6_ABCC8_stability = read_sheet("main", "04h_ABCC8_PP4_stability"),
  E7_GEFOS_all_priors = read_sheet("main", "09_GEFOS_coloc_all_priors"),
  E8_snRNA_sample_QC = read_sheet("snrna", "01_sample_QC"),
  E9_snRNA_pseudobulk_QC = read_sheet("snrna", "09_pseudobulk_QC"),
  E10_snRNA_model = read_sheet("snrna", "10_README"),
  # [FORMER ERROR POINT]
  # Excel sheet names are limited to 31 characters. The actual workbook sheets
  # are the truncated names below, not the longer descriptive names.
E11_SuSiE_complete_region_audit = read_sheet(
    "main", "04i_coloc_SuSiE_complete_region"
  ),
E12_SuSiE_prior_applicability = read_sheet(
    "main", "04j_coloc_SuSiE_prior_applicabi"
  )
)
supplementary_tables <- c(supplementary_tables, additional_evidence)
catalog <- rbind(catalog, data.frame(
  Item = paste("Additional evidence", seq_along(additional_evidence)),
  Chinese_name = names(additional_evidence),
  English_name = names(additional_evidence),
  Output_file_or_sheet = names(additional_evidence),
  Manuscript_location = "Supplementary material",
  stringsAsFactors = FALSE
))

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

# [NEW: final submission outputs] Build all submission tables and figures from
# the same completed analysis release. The module is also runnable on its own.
# [UPDATED 2026-09-17] Pass this run's source root explicitly; do not select by cwd.
local({
  previous_source <- Sys.getenv("IR_SUBMISSION_SOURCE_DIR", unset = NA_character_)
  on.exit(if (is.na(previous_source)) Sys.unsetenv("IR_SUBMISSION_SOURCE_DIR") else
    Sys.setenv(IR_SUBMISSION_SOURCE_DIR = previous_source))
  Sys.setenv(IR_SUBMISSION_SOURCE_DIR = config$project_dir)
  source(file.path(release_root, "scripts", "04_submission_outputs.R"))
})

# [NEW: Reviewer 2] Systematically evaluate expressed genes at all six loci.
sys.source(locus_to_gene_script, envir = new.env(parent = globalenv()))


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

# [NEW] Final runtime provenance, after all analysis modules have completed.
capture.output(sessionInfo(), file = file.path(result_dir, "sessionInfo.txt"))
message("Public-release workflow completed: ", output_dir)

