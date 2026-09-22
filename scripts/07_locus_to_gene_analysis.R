options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
  library(VariantAnnotation)
  library(rtracklayer)
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  library(AnnotationFilter)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(coloc)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "reproducibility_package/scripts/07_locus_to_gene_analysis.R"
root <- "."
config <- list(
  files = list(
    outcomes = c(
      ALM = file.path("outcome", "ebi-a-GCST90000025_ALM.vcf"),
      GRIP = file.path("outcome", "ukb-b-10215_right_grip.vcf"),
      WALK = file.path("outcome", "ukb-b-4711_walkpace.vcf")
    ),
    gtex_muscle = file.path("eQTL", "Muscle_Skeletal.tsv.gz")
  ),
  outcome_n = c(ALM = 450243, GRIP = 461089, WALK = 459915)
)
main_book <- file.path(root, "results_public_release", "manuscript_analysis_tables.xlsx")
output_file <- file.path(root, "results_public_release", "TierA_locus_to_gene_analysis.xlsx")

strip_version <- function(x) sub("\\..*$", "", as.character(x))
normal_p <- function(beta, se) 2 * pnorm(-abs(beta / se))
is_palindromic <- function(a1, a2) {
  paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
}

tier <- as.data.table(read.xlsx(main_book, "04c_coloc_SuSiE_TierA_QC"))
tier <- tier[, .(
  locus_gene = gene_symbol,
  locus_gene_ensembl = strip_version(gene_ensembl),
  trait,
  chromosome = as.character(chromosome),
  gene_start = as.integer(region_start),
  gene_end = as.integer(region_end)
)]
tier[, `:=`(
  locus_start_grch38 = pmax(1L, gene_start - 1000000L),
  locus_end_grch38 = gene_end + 1000000L,
  locus_id = paste(locus_gene, trait, sep = "_")
)]

all_genes <- ensembldb::genes(EnsDb.Hsapiens.v86, return.type = "GRanges")
all_gene_dt <- as.data.table(as.data.frame(all_genes))[, .(
  gene_ensembl = strip_version(gene_id),
  chromosome = as.character(seqnames),
  gene_start = as.integer(start),
  gene_end = as.integer(end),
  strand = as.character(strand),
  gene_biotype = as.character(gene_biotype)
)]
all_gene_dt[, tss := fifelse(strand == "-", gene_end, gene_start)]
gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(all_gene_dt$gene_ensembl),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)
gene_symbols <- as.data.table(gene_symbols)[!is.na(SYMBOL)]
gene_symbols[, ENSEMBL := strip_version(ENSEMBL)]
gene_symbols <- unique(gene_symbols, by = "ENSEMBL")
all_gene_dt <- merge(all_gene_dt, gene_symbols,
                     by.x = "gene_ensembl", by.y = "ENSEMBL", all.x = TRUE)
setnames(all_gene_dt, "SYMBOL", "gene_symbol")

chain_gz <- file.path(root, "hg38ToHg19.over.chain.gz")
chain_file <- tempfile(fileext = ".over.chain")
in_con <- gzfile(chain_gz, "rb")
out_con <- file(chain_file, "wb")
while (length(chunk <- readBin(in_con, "raw", n = 1024L * 1024L))) {
  writeBin(chunk, out_con)
}
close(in_con); close(out_con)
chain <- tryCatch(import.chain(chain_file), finally = unlink(chain_file))
locus_gr38 <- GRanges(
  seqnames = tier$chromosome,
  ranges = IRanges(tier$locus_start_grch38, tier$locus_end_grch38),
  locus_id = tier$locus_id
)
seqlevelsStyle(locus_gr38) <- "UCSC"
lifted <- liftOver(locus_gr38, chain)
tier[, `:=`(locus_start_grch37 = NA_integer_, locus_end_grch37 = NA_integer_)]
for (i in seq_along(lifted)) {
  if (length(lifted[[i]]) == 0L) next
  widest <- lifted[[i]][which.max(width(lifted[[i]]))]
  tier[i, `:=`(
    locus_start_grch37 = as.integer(start(widest)),
    locus_end_grch37 = as.integer(end(widest))
  )]
}
if (anyNA(tier$locus_start_grch37) || anyNA(tier$locus_end_grch37)) {
  stop("At least one Tier A locus could not be lifted from GRCh38 to GRCh37.")
}

gtex_columns <- c(
  "variant", "r2", "pvalue", "molecular_trait_object_id",
  "molecular_trait_id", "maf", "gene_id", "median_tpm", "beta", "se",
  "an", "ac", "chromosome", "position", "ref", "alt", "type", "rsid"
)

read_eqtl_locus <- function(chr, start38, end38) {
  query <- GRanges(as.character(chr), IRanges(start38, end38))
  lines <- scanTabix(config$files$gtex_muscle, param = query)[[1]]
  if (length(lines) == 0L) return(data.table())
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
              sep = "\t", showProgress = FALSE)
  setnames(dt, gtex_columns)
  dt[, gene_ensembl := strip_version(gene_id)]
  dt[, .(
    SNP = as.character(rsid),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pval = as.numeric(pvalue),
    eaf = as.numeric(ac) / as.numeric(an),
    maf = as.numeric(maf),
    effect_allele = toupper(as.character(alt)),
    other_allele = toupper(as.character(ref)),
    median_tpm = as.numeric(median_tpm),
    position_grch38 = as.integer(position),
    gene_ensembl
  )][grepl("^rs", SNP) & is.finite(beta) & is.finite(se) & se > 0 &
       is.finite(eaf) & eaf > 0 & eaf < 1]
}

scan_plain_vcf_snps <- function(file, wanted_snps, fixed_n) {
  con <- file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- NULL
  retained <- list()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    if (is.null(header)) {
      header_line <- lines[startsWith(lines, "#CHROM")]
      if (length(header_line)) header <- strsplit(header_line[1], "\t", fixed = TRUE)[[1]]
    }
    lines <- lines[!startsWith(lines, "#")]
    if (!length(lines)) next
    ids <- tstrsplit(lines, "\t", fixed = TRUE, keep = 3L)[[1]]
    keep <- ids %chin% wanted_snps
    if (any(keep)) retained[[length(retained) + 1L]] <- lines[keep]
  }
  lines <- unlist(retained, use.names = FALSE)
  if (!length(lines)) return(data.table())
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
              sep = "\t", showProgress = FALSE)
  setnames(dt, header)
  sample_col <- names(dt)[length(names(dt))]
  rows <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    fields <- strsplit(as.character(dt$FORMAT[i]), ":", fixed = TRUE)[[1]]
    values <- strsplit(as.character(dt[[sample_col]][i]), ":", fixed = TRUE)[[1]]
    value <- function(name) {
      j <- match(name, fields)
      if (is.na(j) || j > length(values)) return(NA_real_)
      suppressWarnings(as.numeric(values[j]))
    }
    rows[[i]] <- list(beta = value("ES"), se = value("SE"),
                      eaf = value("AF"), N = value("SS"))
  }
  geno <- rbindlist(rows)
  out <- data.table(
    SNP = as.character(dt$ID),
    beta = geno$beta,
    se = geno$se,
    eaf = geno$eaf,
    N = geno$N,
    effect_allele = toupper(as.character(dt$ALT)),
    other_allele = toupper(as.character(dt$REF)),
    chromosome = as.character(dt$`#CHROM`),
    position_grch37 = as.integer(dt$POS)
  )
  out[!is.finite(N), N := fixed_n]
  out[, `:=`(pval = normal_p(beta, se), maf = pmin(eaf, 1 - eaf))]
  out[grepl("^rs", SNP) & is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(eaf) & eaf > 0 & eaf < 1 & !duplicated(SNP)]
}

align_region <- function(e, g) {
  common <- intersect(e$SNP, g$SNP)
  if (length(common) < 50L) return(NULL)
  e <- e[match(common, SNP)]
  g <- g[match(common, SNP)]
  keep <- is.finite(e$maf) & e$maf >= 0.01 & e$maf < 0.5 &
    is.finite(g$maf) & g$maf >= 0.01 & g$maf < 0.5
  e <- e[keep]; g <- g[keep]
  same <- e$effect_allele == g$effect_allele & e$other_allele == g$other_allele
  flip <- e$effect_allele == g$other_allele & e$other_allele == g$effect_allele
  ambiguous <- is_palindromic(e$effect_allele, e$other_allele) &
    (e$maf > 0.42 | g$maf > 0.42)
  keep <- (same | flip) & !ambiguous
  e <- e[keep]; g <- g[keep]; flip <- flip[keep]
  if (any(flip)) {
    e[flip, `:=`(beta = -beta, eaf = 1 - eaf)]
  }
  if (nrow(e) < 50L) return(NULL)
  list(e = e, g = g)
}

run_coloc <- function(e, g, trait) {
  d1 <- list(snp = e$SNP, beta = e$beta, varbeta = e$se^2,
             MAF = e$maf, N = 706, type = "quant", sdY = 1)
  sdY <- sqrt(median(2 * g$maf * (1 - g$maf) * g$N * g$se^2,
                     na.rm = TRUE))
  if (!is.finite(sdY) || sdY <= 0) return(data.table())
  d2 <- list(snp = g$SNP, beta = g$beta, varbeta = g$se^2,
             MAF = g$maf, N = median(g$N, na.rm = TRUE),
             type = "quant", sdY = sdY)
  rbindlist(lapply(c(1e-6, 1e-5, 1e-4), function(p12) {
    invisible(capture.output(
      fit <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12)
    ))
    s <- fit$summary
    data.table(
      p12_requested = p12,
      n_shared_snps = nrow(e),
      PP0 = unname(s["PP.H0.abf"]),
      PP1 = unname(s["PP.H1.abf"]),
      PP2 = unname(s["PP.H2.abf"]),
      PP3 = unname(s["PP.H3.abf"]),
      PP4 = unname(s["PP.H4.abf"])
    )
  }))
}

all_results <- list()
eligibility <- list()
eqtl_by_locus <- lapply(seq_len(nrow(tier)), function(i) {
  loc <- tier[i]
  read_eqtl_locus(loc$chromosome, loc$locus_start_grch38,
                  loc$locus_end_grch38)
})
names(eqtl_by_locus) <- tier$locus_id
gwas_by_trait <- list()
for (trait_name in unique(tier$trait)) {
  message("Scanning ", trait_name, " GWAS once for all requested loci")
  trait_loci <- tier[trait == trait_name, locus_id]
  wanted <- unique(unlist(lapply(eqtl_by_locus[trait_loci], `[[`, "SNP")))
  gwas_by_trait[[trait_name]] <- scan_plain_vcf_snps(
    config$files$outcomes[[trait_name]], wanted,
    config$outcome_n[[trait_name]]
  )
}
for (i in seq_len(nrow(tier))) {
  loc <- tier[i]
  message("Locus ", loc$locus_id)
  eqtl_all <- eqtl_by_locus[[loc$locus_id]]
  genes_here <- all_gene_dt[
    chromosome == loc$chromosome &
      tss >= loc$locus_start_grch38 & tss <= loc$locus_end_grch38
  ]
  gwas <- gwas_by_trait[[loc$trait]][SNP %in% eqtl_all$SNP]

  for (j in seq_len(nrow(genes_here))) {
    gene <- genes_here[j]
    e <- eqtl_all[gene_ensembl == gene$gene_ensembl]
    median_tpm <- if (nrow(e)) median(e$median_tpm, na.rm = TRUE) else NA_real_
    n_eqtl <- nrow(e)
    aligned <- if (nrow(e)) align_region(e, gwas) else NULL
    n_shared <- if (is.null(aligned)) 0L else nrow(aligned$e)
    eligible <- is.finite(median_tpm) && median_tpm >= 1 && n_shared >= 50L
    reason <- if (!is.finite(median_tpm)) {
      "no_GTEx_v8_skeletal_muscle_eQTL_rows_in_fixed_locus"
    } else if (median_tpm < 1) {
      "median_TPM_below_1"
    } else if (n_shared < 50L) {
      "fewer_than_50_harmonised_SNPs"
    } else {
      "eligible"
    }
    eligibility[[length(eligibility) + 1L]] <- data.table(
      locus_id = loc$locus_id,
      locus_gene = loc$locus_gene,
      trait = loc$trait,
      chromosome = loc$chromosome,
      locus_start_grch38 = loc$locus_start_grch38,
      locus_end_grch38 = loc$locus_end_grch38,
      gene_symbol = gene$gene_symbol,
      gene_ensembl = gene$gene_ensembl,
      gene_biotype = gene$gene_biotype,
      tss_grch38 = gene$tss,
      median_tpm = median_tpm,
      n_eqtl_rows_in_fixed_locus = n_eqtl,
      n_harmonised_snps = n_shared,
      eligible = eligible,
      exclusion_reason = reason,
      nominated_tier_a_gene = gene$gene_ensembl == loc$locus_gene_ensembl
    )
    if (!eligible) next
    fit <- run_coloc(aligned$e, aligned$g, loc$trait)
    fit[, `:=`(
      locus_id = loc$locus_id,
      locus_gene = loc$locus_gene,
      trait = loc$trait,
      gene_symbol = gene$gene_symbol,
      gene_ensembl = gene$gene_ensembl,
      gene_biotype = gene$gene_biotype,
      median_tpm = median_tpm,
      lead_eqtl_snp = aligned$e$SNP[which.min(aligned$e$pval)],
      lead_eqtl_p = min(aligned$e$pval, na.rm = TRUE),
      nominated_tier_a_gene = gene$gene_ensembl == loc$locus_gene_ensembl
    )]
    all_results[[length(all_results) + 1L]] <- fit
  }
}

eligibility <- rbindlist(eligibility, fill = TRUE)
results <- rbindlist(all_results, fill = TRUE)
setcolorder(results, c(
  "locus_id", "locus_gene", "trait", "gene_symbol", "gene_ensembl",
  "gene_biotype", "median_tpm", "nominated_tier_a_gene", "lead_eqtl_snp",
  "lead_eqtl_p", "p12_requested", "n_shared_snps",
  "PP0", "PP1", "PP2", "PP3", "PP4"
))
primary <- results[p12_requested == 1e-5]
primary[, rank_PP4_within_locus := frank(-PP4, ties.method = "min"), by = locus_id]
primary[, `:=`(
  passes_PP4_0_80 = PP4 >= 0.80,
  PP4_greater_than_PP3 = PP4 > PP3
)]
summary <- primary[, .(
  genes_with_TPM_ge_1_and_50_SNPs = .N,
  genes_PP4_ge_0_80 = sum(PP4 >= 0.80),
  top_gene = gene_symbol[which.max(PP4)],
  top_gene_PP4 = max(PP4),
  nominated_gene = locus_gene[1],
  nominated_gene_PP4 = PP4[nominated_tier_a_gene][1],
  nominated_gene_rank = rank_PP4_within_locus[nominated_tier_a_gene][1],
  competing_gene_stronger_than_nominated = any(
    !nominated_tier_a_gene & PP4 > PP4[nominated_tier_a_gene][1], na.rm = TRUE
  )
), by = .(locus_id, trait)]

readme <- data.table(
  item = c(
    "Purpose", "Fixed locus", "Regional gene universe", "Expression threshold",
    "Variant threshold", "eQTL source", "Outcome source", "Colocalization model",
    "Primary prior", "Sensitivity priors", "Interpretation boundary"
  ),
  value = c(
    "Compare the nominated Tier A gene with all adequately expressed regional genes at each of six loci.",
    "Tier A gene body plus or minus 1 Mb in GRCh38; the same fixed interval is used for every gene within a locus.",
    "All EnsDb.Hsapiens.v86 genes whose transcription start site lies within the fixed locus.",
    "GTEx v8 skeletal-muscle median TPM at least 1.",
    "At least 50 rsID-matched, allele-harmonised SNPs with MAF at least 0.01 in the fixed locus.",
    "GTEx v8 skeletal-muscle full regional cis-eQTL summary statistics, GRCh38, n=706.",
    "The same GRCh37 GWAS and sample size used for the nominated Tier A gene.",
    "coloc.abf single-signal analysis with p1=p2=1e-4; this is a locus-to-gene sensitivity analysis, not proof of mediation.",
    "p12=1e-5.",
    "p12=1e-6 and 1e-4.",
    "A high PP4 identifies a transcript whose eQTL and outcome associations are compatible with a shared regional signal. It does not prove that the transcript is the causal mediator, and correlated eQTLs can support multiple genes."
  )
)

wb <- createWorkbook()
header_style <- createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF",
                            textDecoration = "bold", halign = "center")
sub_style <- createStyle(fgFill = "#D9EAF7", textDecoration = "bold")
add_sheet <- function(name, x, freeze = TRUE) {
  addWorksheet(wb, name)
  writeData(wb, name, x, headerStyle = header_style, withFilter = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(x)), widths = "auto")
  if (freeze) freezePane(wb, name, firstRow = TRUE)
  addStyle(wb, name, sub_style, rows = 1, cols = seq_len(ncol(x)), gridExpand = TRUE)
}
add_sheet("README", readme, FALSE)
add_sheet("Locus_summary", summary)
add_sheet("Primary_p12_1e-5", primary[order(locus_id, rank_PP4_within_locus)])
add_sheet("All_priors", results[order(locus_id, gene_symbol, p12_requested)])
add_sheet("Gene_eligibility", eligibility[order(locus_id, -nominated_tier_a_gene,
                                                 gene_symbol)])
add_sheet("Locus_definitions", tier)
saveWorkbook(wb, output_file, overwrite = TRUE)
message("Wrote ", output_file)
