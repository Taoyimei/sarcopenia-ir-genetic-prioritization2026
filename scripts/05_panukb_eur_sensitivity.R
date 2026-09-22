## Pan-UK Biobank European-ancestry sensitivity analysis
## Reviewer 3, major comment 1
##
## Scope
## - Retains the prespecified discovery GWAS as the primary analysis.
## - Repeats candidate-gene cis-MR for the two matching Pan-UKB phenotypes.
## - Repeats single-signal colocalisation for the 20 originally FDR-significant
##   GRIP/WALK gene-trait pairs, including SMAD3-GRIP and RXRA-WALK.
## - Uses only EUR-specific Pan-UKB estimates and excludes low-confidence EUR
##   variants. Pan-ancestry meta-analysis estimates are not mixed with the
##   European GTEx eQTL and 1000 Genomes EUR LD framework.
## - Converts GTEx v8 GRCh38 positions to Pan-UKB GRCh37 with the UCSC
##   hg38ToHg19 chain before positional and allele harmonisation.
## - Does not replace or silently overwrite the primary results.

required_packages <- c(
  "data.table", "openxlsx", "Rsamtools", "GenomicRanges", "IRanges", "rtracklayer",
  "TwoSampleMR", "ieugwasr", "coloc"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

options(timeout = 3600)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

project_root <- arg_value("project-root", ".")
analysis_tables <- arg_value(
  "analysis-tables",
  file.path(
    project_root, "★一修", "最终统一版_小意见核查修订_20260913",
    "04_GitHub复现包", "reproducibility_package", "processed_results",
    "manuscript_analysis_tables.xlsx"
  )
)
output_dir <- arg_value(
  "output-dir",
  file.path(
    project_root, "★一修", "最终统一版_小意见核查修订_20260913",
    "04_GitHub复现包", "reproducibility_package", "processed_results"
  )
)
output_xlsx <- file.path(output_dir, "PanUKB_EUR_GWAS_sensitivity.xlsx")
cache_dir <- arg_value("cache-dir", file.path(tempdir(), "panukb_eur_cache"))
chain_file <- arg_value("chain-file", file.path(project_root, "hg38ToHg19.over.chain.gz"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

eqtl_cache <- file.path(project_root, "eqtl_cache_gtex_v8_local")
plink_bin <- file.path(project_root, "plink", "plink", "plink.exe")
ld_prefix <- file.path(project_root, "1000G_EUR", "1000G.EUR.QC")

panukb <- data.table(
  trait = c("GRIP", "WALK"),
  phenocode = c("47", "924"),
  modifier = c("irnt", "none"),
  description = c("Hand grip strength (right)", "Usual walking pace"),
  n_eur = c(418827L, 417933L),
  phenotype_qc_eur = c("PASS", "PASS"),
  url = c(
    paste0(
      "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/",
      "continuous-47-both_sexes-irnt.tsv.bgz"
    ),
    paste0(
      "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/",
      "continuous-924-both_sexes.tsv.bgz"
    )
  )
)

panukb_columns <- c(
  "chr", "pos", "ref", "alt", "af_meta_hq", "beta_meta_hq",
  "se_meta_hq", "neglog10_pval_meta_hq",
  "neglog10_pval_heterogeneity_hq", "af_meta", "beta_meta", "se_meta",
  "neglog10_pval_meta", "neglog10_pval_heterogeneity", "af_AFR",
  "af_AMR", "af_CSA", "af_EAS", "af_EUR", "af_MID", "beta_AFR",
  "beta_AMR", "beta_CSA", "beta_EAS", "beta_EUR", "beta_MID", "se_AFR",
  "se_AMR", "se_CSA", "se_EAS", "se_EUR", "se_MID",
  "neglog10_pval_AFR", "neglog10_pval_AMR", "neglog10_pval_CSA",
  "neglog10_pval_EAS", "neglog10_pval_EUR", "neglog10_pval_MID",
  "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA",
  "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID"
)
panukb_columns_walk <- panukb_columns[
  !grepl("_CSA$", panukb_columns)
]

stopifnot(file.exists(analysis_tables), dir.exists(eqtl_cache), file.exists(chain_file))
stopifnot(file.exists(plink_bin), file.exists(paste0(ld_prefix, ".1.bed")))
chain_import_file <- chain_file
if (grepl("\\.gz$", chain_file, ignore.case = TRUE)) {
  chain_import_file <- file.path(cache_dir, "hg38ToHg19.over.chain")
  if (!file.exists(chain_import_file)) {
    input_connection <- gzfile(chain_file, "rb")
    output_connection <- file(chain_import_file, "wb")
    repeat {
      bytes <- readBin(input_connection, "raw", n = 1024L * 1024L)
      if (!length(bytes)) break
      writeBin(bytes, output_connection)
    }
    close(input_connection)
    close(output_connection)
  }
}
hg38_to_hg19 <- rtracklayer::import.chain(chain_import_file)

lift_points_to_hg19 <- function(x) {
  if (!nrow(x)) return(x)
  gr <- GRanges(
    seqnames = paste0("chr", x$chr),
    ranges = IRanges(start = x$pos, end = x$pos),
    row_id = seq_len(nrow(x))
  )
  lifted <- liftOver(gr, hg38_to_hg19)
  keep <- which(lengths(lifted) == 1L)
  if (!length(keep)) return(x[0])
  mapped <- unlist(lifted[keep], use.names = FALSE)
  out <- copy(x[keep])
  out[, `:=`(
    chr_hg38 = chr, pos_hg38 = pos,
    chr = sub("^chr", "", as.character(seqnames(mapped))),
    pos = start(mapped)
  )]
  out
}

gene_positions <- as.data.table(openxlsx::read.xlsx(
  analysis_tables, sheet = "01_candidate_gene_coordinates"
))
original_mr <- as.data.table(openxlsx::read.xlsx(
  analysis_tables, sheet = "02_cisMR_all_results"
))
original_sig <- as.data.table(openxlsx::read.xlsx(
  analysis_tables, sheet = "02_cisMR_FDR_significant"
))
original_coloc <- as.data.table(openxlsx::read.xlsx(
  analysis_tables, sheet = "04_coloc_primary_prior"
))

normal_p <- function(beta, se) {
  2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
}

read_and_clump_gene <- function(gene_row) {
  gene_id <- as.character(gene_row$gene_ensembl)
  gene_symbol <- as.character(gene_row$gene_symbol)
  chromosome <- as.character(gene_row$chr)
  cache_file <- file.path(eqtl_cache, paste0(gene_id, "_muscle.rds"))
  if (!file.exists(cache_file) || !chromosome %in% as.character(1:22)) return(NULL)
  x <- as.data.table(readRDS(cache_file))
  required <- c(
    "SNP", "beta", "se", "pval", "eaf", "maf", "effect_allele",
    "other_allele", "pos", "chr", "gene_ensembl"
  )
  if (!all(required %in% names(x))) return(NULL)
  x[, F_stat := (beta / se)^2]
  x <- x[
    is.finite(beta) & is.finite(se) & se > 0 & is.finite(eaf) &
      eaf > 0 & eaf < 1 & pval < 5e-8 & F_stat >= 10
  ]
  if (!nrow(x)) return(NULL)
  x <- unique(x, by = "SNP")
  if (nrow(x) > 1L) {
    clumped <- tryCatch(
      ieugwasr::ld_clump(
        data.frame(rsid = x$SNP, pval = x$pval, id = gene_symbol),
        clump_kb = 1000, clump_r2 = 0.001, clump_p = 5e-8,
        bfile = paste0(ld_prefix, ".", chromosome),
        plink_bin = plink_bin
      ),
      error = function(e) NULL
    )
    if (is.null(clumped) || !nrow(clumped)) return(NULL)
    x <- x[SNP %in% clumped$rsid]
  }
  x[, `:=`(gene_symbol = gene_symbol, gene_ensembl = gene_id)]
  x[]
}

instrument_cache <- file.path(cache_dir, "prespecified_instruments.rds")
if (file.exists(instrument_cache)) {
  message("Reading cached prespecified instruments...")
  instruments <- as.data.table(readRDS(instrument_cache))
} else {
  message("Reconstructing the prespecified skeletal-muscle instruments...")
  instrument_list <- lapply(
    seq_len(nrow(gene_positions)),
    function(i) read_and_clump_gene(gene_positions[i])
  )
  instruments <- rbindlist(instrument_list, fill = TRUE)
  saveRDS(instruments, instrument_cache)
}
if (!nrow(instruments)) stop("No instruments could be reconstructed.")
if (!all(c("chr_hg38", "pos_hg38") %in% names(instruments))) {
  instruments <- lift_points_to_hg19(instruments)
}
setcolorder(
  instruments,
  c("gene_ensembl", "gene_symbol", "SNP", "chr", "pos", "effect_allele",
    "other_allele", "eaf", "beta", "se", "pval", "F_stat")
)

make_query_ranges <- function(instruments, coloc_targets) {
  inst_gr <- GRanges(
    seqnames = as.character(instruments$chr),
    ranges = IRanges(start = instruments$pos, end = instruments$pos)
  )
  target_pos <- gene_positions[
    gene_ensembl %in% coloc_targets$gene_ensembl,
    .(chr = as.character(chr), start = start, end = end)
  ]
  target_gr38 <- GRanges(
    seqnames = paste0("chr", target_pos$chr),
    ranges = IRanges(target_pos$start, target_pos$end)
  )
  target_lifted <- unlist(liftOver(target_gr38, hg38_to_hg19), use.names = FALSE)
  coloc_gr <- GRanges(
    seqnames = sub("^chr", "", as.character(seqnames(target_lifted))),
    ranges = IRanges(
      pmax(1L, start(target_lifted) - 1000000L),
      end(target_lifted) + 1000000L
    )
  )
  reduce(c(inst_gr, coloc_gr), min.gapwidth = 50000L)
}

scan_panukb <- function(url, ranges, trait) {
  raw <- vector("list", length(ranges))
  for (k in seq_along(ranges)) {
    r <- ranges[k]
    cache_file <- file.path(
      cache_dir,
      sprintf(
        "%s_%s_%s_%s.rds", trait, as.character(seqnames(r)), start(r), end(r)
      )
    )
    if (file.exists(cache_file)) {
      raw[[k]] <- readRDS(cache_file)
    } else {
      query_result <- NULL
      for (attempt in 1:6) {
        query_result <- tryCatch(
          scanTabix(TabixFile(url), param = r)[[1]],
          error = function(e) NULL
        )
        if (!is.null(query_result)) break
        message("  retrying ", trait, " region ", k, " (attempt ",
                attempt + 1L, "/6)")
        Sys.sleep(5 * attempt)
      }
      if (is.null(query_result)) {
        stop("Pan-UKB query failed after six attempts: ", trait,
             " region ", k)
      }
      raw[[k]] <- query_result
      saveRDS(raw[[k]], cache_file)
    }
    if (k %% 10L == 0L || k == length(ranges)) {
      message("  ", trait, " regions read: ", k, "/", length(ranges))
    }
  }
  lines <- unique(unlist(raw, use.names = FALSE))
  if (!length(lines)) return(data.table())
  x <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t")
  expected_columns <- if (trait == "WALK") panukb_columns_walk else panukb_columns
  if (ncol(x) != length(expected_columns)) {
    stop("Unexpected Pan-UKB column count for ", trait, ": ", ncol(x))
  }
  setnames(x, expected_columns)
  x[, `:=`(
    chr = as.character(chr), pos = as.integer(pos),
    af_EUR = as.numeric(af_EUR), beta_EUR = as.numeric(beta_EUR),
    se_EUR = as.numeric(se_EUR),
    neglog10_pval_EUR = as.numeric(neglog10_pval_EUR),
    low_confidence_EUR = tolower(as.character(low_confidence_EUR))
  )]
  x <- x[
    is.finite(beta_EUR) & is.finite(se_EUR) & se_EUR > 0 &
      is.finite(af_EUR) & af_EUR > 0 & af_EUR < 1 &
      low_confidence_EUR %in% c("false", "0")
  ]
  unique(x, by = c("chr", "pos", "ref", "alt"))
}

harmonise_gene <- function(inst, outcome, trait, n_outcome) {
  z <- merge(
    inst, outcome,
    by = c("chr", "pos"), all = FALSE, allow.cartesian = TRUE,
    suffixes = c("_exposure", "_outcome")
  )
  z <- z[
    (effect_allele == alt & other_allele == ref) |
      (effect_allele == ref & other_allele == alt)
  ]
  if (!nrow(z)) return(NULL)
  z[, swapped := effect_allele == ref & other_allele == alt]
  z[, `:=`(
    beta_outcome_aligned = fifelse(swapped, -beta_EUR, beta_EUR),
    eaf_outcome_aligned = fifelse(swapped, 1 - af_EUR, af_EUR),
    pval_outcome = pmax(10^(-neglog10_pval_EUR), 1e-300)
  )]
  z <- unique(z, by = "SNP")
  exposure <- data.frame(
    SNP = z$SNP,
    beta.exposure = z$beta,
    se.exposure = z$se,
    effect_allele.exposure = z$effect_allele,
    other_allele.exposure = z$other_allele,
    eaf.exposure = z$eaf,
    pval.exposure = z$pval,
    samplesize.exposure = 706,
    exposure = z$gene_symbol,
    id.exposure = z$gene_symbol
  )
  outcome_df <- data.frame(
    SNP = z$SNP,
    beta.outcome = z$beta_outcome_aligned,
    se.outcome = z$se_EUR,
    effect_allele.outcome = z$effect_allele,
    other_allele.outcome = z$other_allele,
    eaf.outcome = z$eaf_outcome_aligned,
    pval.outcome = z$pval_outcome,
    samplesize.outcome = n_outcome,
    outcome = trait,
    id.outcome = trait
  )
  h <- suppressMessages(TwoSampleMR::harmonise_data(exposure, outcome_df, action = 2))
  h[h$mr_keep %in% TRUE, , drop = FALSE]
}

run_mr <- function(h, gene_id, gene_symbol, trait) {
  if (is.null(h) || !nrow(h)) return(NULL)
  methods <- if (nrow(h) == 1L) "mr_wald_ratio" else "mr_ivw"
  fit <- suppressMessages(TwoSampleMR::mr(h, method_list = methods))
  if (!nrow(fit)) return(NULL)
  data.table(
    gene_ensembl = gene_id,
    gene_symbol = gene_symbol,
    trait = trait,
    method = fit$method[[1]],
    nsnp = fit$nsnp[[1]],
    beta = fit$b[[1]],
    se = fit$se[[1]],
    pval = fit$pval[[1]],
    ci_lower = fit$b[[1]] - 1.96 * fit$se[[1]],
    ci_upper = fit$b[[1]] + 1.96 * fit$se[[1]],
    min_F = min(h$beta.exposure^2 / h$se.exposure^2),
    instrument_snps = paste(h$SNP, collapse = ";")
  )
}

run_coloc <- function(gene_id, gene_symbol, trait, outcome, n_outcome) {
  cache_file <- file.path(eqtl_cache, paste0(gene_id, "_muscle.rds"))
  if (!file.exists(cache_file)) return(NULL)
  e <- as.data.table(readRDS(cache_file))
  e <- lift_points_to_hg19(e)
  g <- outcome[chr == unique(e$chr) & pos >= min(e$pos) & pos <= max(e$pos)]
  if (!nrow(g)) return(NULL)
  z <- merge(
    e, g, by = c("chr", "pos"), all = FALSE, allow.cartesian = TRUE,
    suffixes = c("_eqtl", "_gwas")
  )
  z <- z[
    (effect_allele == alt & other_allele == ref) |
      (effect_allele == ref & other_allele == alt)
  ]
  if (nrow(z) < 50L) return(NULL)
  z[, swapped := effect_allele == ref & other_allele == alt]
  z[, `:=`(
    beta_gwas_aligned = fifelse(swapped, -beta_EUR, beta_EUR),
    maf_gwas = pmin(fifelse(swapped, 1 - af_EUR, af_EUR),
                    1 - fifelse(swapped, 1 - af_EUR, af_EUR)),
    snp_key = paste(chr, pos, ref, alt, sep = ":")
  )]
  z <- z[
    is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(maf) & maf >= 0.01 & maf <= 0.49 &
      is.finite(maf_gwas) & maf_gwas >= 0.01 & maf_gwas <= 0.49
  ]
  z <- unique(z, by = "snp_key")
  if (nrow(z) < 50L) return(NULL)
  d1 <- list(
    snp = z$snp_key, beta = z$beta, varbeta = z$se^2,
    MAF = z$maf, N = 706, type = "quant", sdY = 1
  )
  d2 <- list(
    snp = z$snp_key, beta = z$beta_gwas_aligned, varbeta = z$se_EUR^2,
    MAF = z$maf_gwas, N = n_outcome, type = "quant"
  )
  if (trait == "GRIP") {
    d2$sdY <- 1
  } else {
    d2$sdY <- sqrt(stats::median(
      2 * z$maf_gwas * (1 - z$maf_gwas) * n_outcome * z$se_EUR^2,
      na.rm = TRUE
    ))
  }
  if (!is.finite(d2$sdY) || d2$sdY <= 0) return(NULL)
  rbindlist(lapply(c(1e-6, 1e-5, 1e-4), function(p12) {
    fit <- tryCatch(
      coloc::coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- fit$summary
    data.table(
      gene_ensembl = gene_id, gene_symbol = gene_symbol, trait = trait,
      p12 = p12, n_shared_variants = nrow(z),
      PP0 = unname(s[["PP.H0.abf"]]), PP1 = unname(s[["PP.H1.abf"]]),
      PP2 = unname(s[["PP.H2.abf"]]), PP3 = unname(s[["PP.H3.abf"]]),
      PP4 = unname(s[["PP.H4.abf"]])
    )
  }), fill = TRUE)
}

mr_results <- list()
coloc_results <- list()
harmonised_rows <- list()

for (i in seq_len(nrow(panukb))) {
  trait <- panukb$trait[[i]]
  message("Reading Pan-UKB EUR regions for ", trait, "...")
  coloc_targets <- original_sig[trait == panukb$trait[[i]]]
  ranges <- make_query_ranges(instruments, coloc_targets)
  message("  merged query regions: ", length(ranges))
  outcome <- scan_panukb(panukb$url[[i]], ranges, trait)
  if (!nrow(outcome)) stop("No Pan-UKB records returned for ", trait)

  for (gene_id in unique(instruments$gene_ensembl)) {
    inst <- instruments[gene_ensembl == gene_id]
    gene_symbol <- inst$gene_symbol[[1]]
    h <- harmonise_gene(inst, outcome, trait, panukb$n_eur[[i]])
    fit <- run_mr(h, gene_id, gene_symbol, trait)
    if (!is.null(fit)) mr_results[[length(mr_results) + 1L]] <- fit
    if (!is.null(h) && nrow(h)) {
      harmonised_rows[[length(harmonised_rows) + 1L]] <- data.table(
        gene_ensembl = gene_id, gene_symbol = gene_symbol, trait = trait,
        SNP = h$SNP, beta_exposure = h$beta.exposure,
        se_exposure = h$se.exposure, beta_outcome = h$beta.outcome,
        se_outcome = h$se.outcome, effect_allele = h$effect_allele.exposure,
        other_allele = h$other_allele.exposure,
        eaf_exposure = h$eaf.exposure, eaf_outcome = h$eaf.outcome,
        F_stat = h$beta.exposure^2 / h$se.exposure^2
      )
    }
  }

  for (j in seq_len(nrow(coloc_targets))) {
    fit <- run_coloc(
      coloc_targets$gene_ensembl[[j]], coloc_targets$gene_symbol[[j]], trait,
      outcome, panukb$n_eur[[i]]
    )
    if (!is.null(fit) && nrow(fit)) {
      coloc_results[[length(coloc_results) + 1L]] <- fit
    }
  }
}

mr <- rbindlist(mr_results, fill = TRUE)
harmonised <- rbindlist(harmonised_rows, fill = TRUE)
coloc_all <- rbindlist(coloc_results, fill = TRUE)

mr[, fdr_within_trait := p.adjust(pval, method = "BH"), by = trait]
mr[, mr_fdr_significant := fdr_within_trait < 0.05]

comparison <- merge(
  original_mr[trait %in% c("GRIP", "WALK"), .(
    gene_ensembl, gene_symbol, trait, original_method = method,
    original_nsnp = nsnp, original_beta = beta, original_se = se,
    original_pval = pval, original_fdr = fdr_within_trait
  )],
  mr[, .(
    gene_ensembl, trait, panukb_method = method, panukb_nsnp = nsnp,
    panukb_beta = beta, panukb_se = se, panukb_pval = pval,
    panukb_fdr = fdr_within_trait, panukb_instrument_snps = instrument_snps
  )],
  by = c("gene_ensembl", "trait"), all = TRUE
)
comparison[, `:=`(
  same_direction = is.finite(original_beta) & is.finite(panukb_beta) &
    sign(original_beta) == sign(panukb_beta),
  original_fdr_significant = is.finite(original_fdr) & original_fdr < 0.05,
  panukb_fdr_significant = is.finite(panukb_fdr) & panukb_fdr < 0.05
)]

if (nrow(coloc_all)) {
  coloc_primary <- coloc_all[p12 == 1e-5]
  coloc_ranges <- coloc_all[, .(
    PP4_min = min(PP4), PP4_max = max(PP4),
    PP4_primary = PP4[p12 == 1e-5][1],
    PP3_primary = PP3[p12 == 1e-5][1],
    n_shared_variants = n_shared_variants[1]
  ), by = .(gene_ensembl, gene_symbol, trait)]
} else {
  coloc_primary <- data.table()
  coloc_ranges <- data.table()
}

tier_a_pairs <- data.table(
  gene_symbol = c("SMAD3", "RXRA"), trait = c("GRIP", "WALK")
)
tier_a <- merge(
  tier_a_pairs,
  comparison[, .(
    gene_symbol, trait, original_beta, original_pval, original_fdr,
    panukb_beta, panukb_pval, panukb_fdr, same_direction,
    panukb_instrument_snps
  )],
  by = c("gene_symbol", "trait"), all.x = TRUE
)
tier_a <- merge(tier_a, coloc_ranges, by = c("gene_symbol", "trait"), all.x = TRUE)
tier_a[, panukb_internal_tier_a :=
  is.finite(panukb_fdr) & panukb_fdr < 0.05 &
    is.finite(PP4_primary) & PP4_primary >= 0.80]

counts <- rbindlist(list(
  mr[, .(
    metric = c("eligible genes with valid Pan-UKB MR", "Pan-UKB MR FDR significant"),
    value = c(.N, sum(mr_fdr_significant)), denominator = c(.N, .N)
  ), by = trait],
  comparison[!is.na(original_beta) & !is.na(panukb_beta), .(
    metric = c("direction concordant among shared tests",
               "original FDR-significant pairs retained at Pan-UKB FDR"),
    value = c(sum(same_direction),
              sum(original_fdr_significant & panukb_fdr_significant)),
    denominator = c(.N, sum(original_fdr_significant))
  ), by = trait],
  data.table(
    trait = c("GRIP", "WALK"),
    metric = "original FDR-significant pairs with Pan-UKB coloc attempted",
    value = c(sum(coloc_ranges$trait == "GRIP"), sum(coloc_ranges$trait == "WALK")),
    denominator = c(sum(original_sig$trait == "GRIP"),
                    sum(original_sig$trait == "WALK"))
  )
), fill = TRUE)

scope <- data.table(
  item = c(
    "analysis role", "ancestry", "genome build", "effect allele",
    "variant QC", "MR testing families", "colocalisation family",
    "not performed", "resource chronology", "interpretation"
  ),
  value = c(
    "Updated UK Biobank sensitivity analysis; prespecified primary results retained",
    "Pan-UKB EUR-specific beta, SE and AF only",
    "GTEx v8 GRCh38 positions lifted to Pan-UKB GRCh37 before allele matching",
    "Pan-UKB ALT allele; aligned to the GTEx exposure effect allele",
    "Finite EUR beta/SE/AF; low_confidence_EUR == false; allele-pair match",
    "BH within GRIP and WALK across all valid candidate-gene MR tests",
    "20 originally FDR-significant GRIP/WALK pairs; three P12 priors; no P-value FDR on PP4",
    "Genome-wide MAGMA update; full files total approximately 4.3 GB and are not required for regional MR/coloc sensitivity",
    "Pan-UKB summary statistics first released in 2020 and subsequently updated; cited article published in 2025",
    "Sensitivity to an alternative analysis of overlapping UK Biobank participants, not independent replication"
  )
)

wb <- createWorkbook()
addWorksheet(wb, "Readme")
writeData(wb, "Readme", scope)
addWorksheet(wb, "Phenotype_compatibility")
writeData(wb, "Phenotype_compatibility", panukb)
addWorksheet(wb, "Analysis_counts")
writeData(wb, "Analysis_counts", counts)
addWorksheet(wb, "Instrument_QC")
writeData(wb, "Instrument_QC", instruments)
addWorksheet(wb, "Harmonised_instruments")
writeData(wb, "Harmonised_instruments", harmonised)
addWorksheet(wb, "PanUKB_MR_all")
writeData(wb, "PanUKB_MR_all", mr)
addWorksheet(wb, "MR_comparison")
writeData(wb, "MR_comparison", comparison)
addWorksheet(wb, "Coloc_all_priors")
writeData(wb, "Coloc_all_priors", coloc_all)
addWorksheet(wb, "Coloc_primary")
writeData(wb, "Coloc_primary", coloc_primary)
addWorksheet(wb, "TierA_sensitivity")
writeData(wb, "TierA_sensitivity", tier_a)

for (s in names(wb)) {
  freezePane(wb, s, firstRow = TRUE)
  setColWidths(wb, s, cols = 1:min(30, ncol(readWorkbook(wb, s))), widths = "auto")
}
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

message("Sensitivity analysis written to: ", output_xlsx)
message("Eligible genes after exposure reconstruction: ",
        uniqueN(instruments$gene_ensembl))
message("Valid Pan-UKB MR tests: ", nrow(mr))
message("Pan-UKB colocalisation posterior vectors: ", nrow(coloc_all))
