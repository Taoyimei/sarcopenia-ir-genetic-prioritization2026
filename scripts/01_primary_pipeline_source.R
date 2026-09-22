
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

config <- list(
  project_dir = ".",

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
    eqtlgen = "Significant_cis_eQTLs.gz"
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
    susie_max_snps = 8000L
  ),

  optional = list(
    run_fusion_replication = TRUE,
    run_eqtlgen_comparison = TRUE,
    run_go_enrichment = TRUE,
    make_locus_plots = FALSE,
    coordinates_harmonised_for_plots = FALSE
  ),

  random_seed = 20260719L
)

config$paths <- list(
  results = file.path(config$project_dir, "results_extended"),
  plots = file.path(config$project_dir, "results_extended", "plots"),
  cache = file.path(config$project_dir, "eqtl_cache_gtex_v8_local")
)

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

vcf_to_gwas <- function(path, trait, fixed_n) {
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
  unique(dt, by = "SNP")
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
  unique(dt, by = "SNP")
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

build_candidate_gene_set <- function() {
  log_step("Building the prespecified IR candidate-gene set")
  msig_all <- msigdbr::msigdbr(species = "Homo sapiens")

  canonical_names <- c(
    "KEGG_INSULIN_SIGNALING_PATHWAY",
    "KEGG_TYPE_II_DIABETES_MELLITUS",
    "KEGG_PPAR_SIGNALING_PATHWAY",
    "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE",
    "REACTOME_SIGNALLING_BY_INSULIN_RECEPTOR",
    "REACTOME_PI3K_CASCADE",
    "REACTOME_FOXO_MEDIATED_TRANSCRIPTION",
    "REACTOME_MTOR_SIGNALLING",
    "REACTOME_MTORC1_MEDIATED_SIGNALLING",
    "GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS",
    "GOBP_RESPONSE_TO_INSULIN_STIMULUS",
    "GOBP_GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN_STIMULUS",
    "GOBP_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS",
    "GOBP_GLUCOSE_HOMEOSTASIS",
    "GOBP_PI3K_SIGNALING")
  

  expansion_pattern <- paste(
    c(
      "INSULIN_SIGNALING", "INSULIN_RECEPTOR", "INSULIN_RESISTANCE",
      "GLUCOSE_IMPORT_IN_RESPONSE_TO_INSULIN", "GLUCOSE_HOMEOSTASIS",
      "SLC2A4", "GLUT4", "IRS1", "IRS2", "FOXO_MEDIATED", "MTORC1"
    ),
    collapse = "|")
  
  canonical <- msig_all |>
    dplyr::filter(gs_name %in% canonical_names) |>
    dplyr::mutate(layer = "canonical")

  expanded <- msig_all |>
    dplyr::filter(
      gs_subcollection %in% c(
        "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME", "GO:BP"
      ),
      stringr::str_detect(
        gs_name, stringr::regex(expansion_pattern, ignore_case = TRUE)
      ),
      !gs_name %in% canonical_names
    ) |>
    dplyr::mutate(layer = "expanded")

  pathways <- dplyr::bind_rows(canonical, expanded) |>
    dplyr::distinct(gs_name, gene_symbol, layer)

  membership <- pathways |>
    dplyr::group_by(gene_symbol) |>
    dplyr::summarise(
      in_canonical = any(layer == "canonical"),
      in_expanded = any(layer == "expanded"),
      gene_layer = ifelse(in_canonical, "canonical", "expanded_only"),
      n_pathways = dplyr::n_distinct(gs_name),
      .groups = "drop"
    )

  list(msig_all = msig_all, pathways = pathways, membership = membership)
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
    data.table(pos = e$pos, minus_log10_p = -log10(pmax(e$pval, 1e-300)),
               dataset = "Skeletal-muscle eQTL"),
    data.table(pos = g$pos, minus_log10_p = -log10(pmax(g$pval, 1e-300)),
               dataset = trait)
  )
  p <- ggplot(plot_dt, aes(pos / 1e6, minus_log10_p, colour = dataset)) +
    geom_point(alpha = 0.65, size = 1.3) +
    facet_wrap(~dataset, ncol = 1, scales = "free_y") +
    labs(
      title = paste(gene_symbol, trait, sep = " - "),
      x = "Genomic position (Mb)", y = expression(-log[10](P)), colour = NULL
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

select_susie_input_snps <- function(e, g, bim, max_snps) {
  if (nrow(bim) <= max_snps) {
    return(list(
      snps = bim$SNP,
      status = "complete_ld_region",
      rule = "All LD-reference allele-matched variants were used."
    ))
  }

  score <- data.table(
    SNP = bim$SNP,
    eqtl_p = e[match(bim$SNP, SNP), pval],
    outcome_p = g[match(bim$SNP, SNP), pval]
  )
  score[, `:=`(
    eqtl_rank = rank(eqtl_p, ties.method = "first"),
    outcome_rank = rank(outcome_p, ties.method = "first")
  )]
  score[, best_rank := pmin(eqtl_rank, outcome_rank)]
  score[, worst_rank := pmax(eqtl_rank, outcome_rank)]
  selected <- score[order(best_rank, worst_rank, eqtl_p, outcome_p, SNP)]
  list(
    snps = selected[seq_len(max_snps), SNP],
    status = "trimmed_to_susie_limit",
    rule = paste0(
      "Top ", max_snps,
      " variants by best rank across skeletal-muscle eQTL and outcome GWAS ",
      "P values after LD-reference allele matching."
    )
  )
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
      SNP = e$SNP, pos = e$pos, pval = e$pval,
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
      x = "Genomic position (Mb)",
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
  output_dir = "results_extended",
  output_workbook = file.path(
    "results_extended", "GSE167186_snRNA_celltype_validation.xlsx"
  ),
  output_umap_pdf = file.path(
    "results_extended", "GSE167186_snRNA_celltype_UMAP.pdf"
  ),
  output_umap_png = file.path(
    "results_extended", "GSE167186_snRNA_celltype_UMAP.png"
  ),
  output_dotplot_pdf = file.path(
    "results_extended", "GSE167186_TierA_celltype_expression.pdf"
  ),
  output_dotplot_png = file.path(
    "results_extended", "GSE167186_TierA_celltype_expression.png"
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
  umap.method = "uwot",
  metric = "cosine",
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

sn_umap_cluster <- Seurat::DimPlot(
  sn_object,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.06
) +
  ggplot2::labs(
    title = NULL, x = "UMAP 1", y = "UMAP 2", colour = "Louvain cluster"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(legend.position = "none")

sn_umap_celltype <- Seurat::DimPlot(
  sn_object,
  reduction = "umap",
  group.by = "cell_type",
  raster = TRUE,
  pt.size = 0.06
) +
  ggplot2::labs(
    title = NULL, x = "UMAP 1", y = "UMAP 2", colour = "Cell type"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 7)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))

sn_marker_dotplot <- Seurat::DotPlot(
  sn_object,
  features = unique(unlist(lapply(sn_marker_sets, head, 3L), use.names = FALSE)),
  group.by = "cell_type",
  assay = "RNA",
  dot.scale = 5
) +
  Seurat::RotatedAxis() +
  ggplot2::labs(
    x = "Canonical marker gene",
    y = NULL,
    colour = "Scaled average\nexpression",
    size = "Nuclei detected (%)"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_line(linewidth = 0.2)
  )

sn_umap <- (sn_umap_cluster | sn_umap_celltype) /
  sn_marker_dotplot +
  patchwork::plot_layout(heights = c(1.15, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  sn_config$output_umap_pdf, sn_umap,
  width = 12, height = 8.5, units = "in"
)
ggplot2::ggsave(
  sn_config$output_umap_png, sn_umap,
  width = 12, height = 8.5, units = "in", dpi = 320,
  device = grDevices::png, type = "cairo-png"
)

sn_localization_plot <- ggplot2::ggplot(
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

sn_candidate_de[, gene := factor(gene, levels = sn_targets)]
sn_candidate_de[, cell_type := factor(
  cell_type, levels = rev(names(sn_marker_sets))
)]
sn_age_effect_plot <- ggplot2::ggplot(
  sn_candidate_de,
  ggplot2::aes(x = gene, y = cell_type, fill = logFC_Old_vs_Young)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(
    data = sn_candidate_de[!is.na(candidate_FDR_global) & candidate_FDR_global < 0.05],
    label = "*", size = 4
  ) +
  ggplot2::geom_point(
    data = sn_candidate_de[is.na(PValue)],
    shape = 4, size = 2, stroke = 0.6
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3182BD", mid = "white", high = "#D24B40",
    midpoint = 0, na.value = "#D9D9D9"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "log2FC\nOld vs Young"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  )

sn_dotplot <- sn_localization_plot / sn_age_effect_plot +
  patchwork::plot_layout(heights = c(1.25, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  sn_config$output_dotplot_pdf, sn_dotplot,
  width = 10, height = 10, units = "in"
)
ggplot2::ggsave(
  sn_config$output_dotplot_png, sn_dotplot,
  width = 10, height = 10, units = "in", dpi = 320,
  device = grDevices::png, type = "cairo-png"
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
    "Log normalization, 2,000 variable genes, scaling, 10 PCs, R-native UWOT UMAP with cosine distance, and Louvain clustering at resolution 0.2 reproduce the published framework; the supplementary figure shows unsupervised clusters, final cell-type labels, and the canonical-marker dot plot. No batch correction is imposed because the source study reported no batch effect.",
    "Clusters are assigned to seven prespecified skeletal-muscle cell classes using canonical-marker module scores. A support-cell score of at least 0.25 overrides ambient myofiber signal; otherwise the higher fast/slow myofiber score is used. All scores are exported for audit.",
    "Only ABCC8, SMAD3, MAPK1, YWHAZ, RXRA, and ZBTB7B are used for candidate-set inference.",
    "Mean log-normalized expression and percentage of nuclei detected are descriptive localization summaries, not independent-sample significance tests; both age groups and all seven annotated cell types are displayed.",
    "Raw counts are summed within donor and cell type. edgeR quasi-likelihood models compare Old versus Young using donors, not nuclei, as biological replicates; cell types require at least 20 nuclei per donor and at least 3 donors per group. All 42 planned candidate-gene-by-cell-type comparisons, including untested low-expression entries, are displayed and exported.",
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
    "results_extended", "clinical_translation_public_databases.xlsx"
  ),
  discovery_workbook = file.path(
    "results_extended", "GEO_human_muscle_validation.xlsx"
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

for (trait in target_traits) {
  
  gwas_dt <- copy(gwas_list[[trait]])
  
  magma_pval <- gwas_dt[
    grepl("^rs", SNP) &
      is.finite(pval) &
      pval >= 0 &
      pval <= 1,
    .(
      SNP = as.character(SNP),
      P = pmax(as.numeric(pval), 1e-300)
    )
  ]
  
  ## 每个rsID仅保留P值最小的一条
  setorder(magma_pval, SNP, P)
  magma_pval <- unique(magma_pval, by = "SNP")
  
  fwrite(
    magma_pval,
    file.path(magma_dir, paste0(trait, "_MAGMA_pval.tsv")),
    sep = "\t",
    quote = FALSE)}




magma_exe <- file.path(config$project_dir, "magma.exe")

magma_path <- function(path, mustWork = TRUE) {
  normal_path <- normalizePath(path, winslash = "\\", mustWork = mustWork)
  path_dir <- dirname(normal_path)
  file.path(utils::shortPathName(path_dir), basename(normal_path))}

ref_prefix <- file.path(config$project_dir, "1000G_EUR/1000G.EUR.QC")

ref_chr_prefix <- file.path(
  config$project_dir,
  "1000G_EUR",
  "1000G.EUR.QC.#CHR#")

gene_loc <- file.path(config$project_dir, "NCBI37.3.gene.loc")

annot_prefix <- file.path(magma_dir, "GRCh37_gene_body")
snp_loc <- file.path(magma_dir, "MAGMA_snp_loc.tsv")
bim_files <- paste0(ref_prefix, ".", 1:22, ".bim")
assert_files_exist(bim_files, "1000G EUR per-chromosome BIM")

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
 








