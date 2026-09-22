# [NEW: MAGMA reproducibility] Resolve duplicate rsIDs against the supplied
# 1000G EUR BIM coordinates and unordered allele pairs.  The resolver is kept
# separate from the main release script so it can be tested on a small fixture
# without reading the full GWAS files.

magma_unordered_allele_key <- function(a, b) {
  a <- toupper(trimws(as.character(a)))
  b <- toupper(trimws(as.character(b)))
  invalid <- is.na(a) | is.na(b) | !nzchar(a) | !nzchar(b)
  out <- rep(NA_character_, length(a))
  valid <- !invalid
  out[valid] <- paste(
    pmin(a[valid], b[valid]),
    pmax(a[valid], b[valid]),
    sep = "|"
  )
  out
}

magma_read_bim_reference <- function(bim_files) {
  refs <- lapply(bim_files, function(bim_file) {
    bim <- data.table::fread(
      bim_file,
      header = FALSE,
      select = c(1L, 2L, 4L, 5L, 6L),
      col.names = c("chr", "SNP", "pos", "ref_a1", "ref_a2"),
      showProgress = FALSE
    )
    bim[, `:=`(
      chr = as.character(chr),
      SNP = as.character(SNP),
      pos = as.integer(pos),
      allele_key = magma_unordered_allele_key(ref_a1, ref_a2)
    )]
    bim[
      grepl("^rs", SNP) & is.finite(pos) & !is.na(allele_key),
      .(chr, SNP, pos, allele_key, ref_a1, ref_a2)
    ]
  })
  unique(
    data.table::rbindlist(refs, use.names = TRUE, fill = TRUE),
    by = c("SNP", "chr", "pos", "allele_key")
  )
}

magma_selected_record <- function(row) {
  if (is.null(row) || nrow(row) == 0L) return(NA_character_)
  paste(
    paste0("source_row=", row$source_row[[1L]]),
    paste0("chr=", row$chr[[1L]]),
    paste0("pos=", row$pos[[1L]]),
    paste0("effect_allele=", row$effect_allele[[1L]]),
    paste0("other_allele=", row$other_allele[[1L]]),
    paste0("pval=", format(row$pval[[1L]], digits = 17, scientific = TRUE)),
    sep = ";"
  )
}

resolve_magma_duplicate_rsids <- function(source_dt, reference_dt, trait) {
  source_dt <- data.table::copy(data.table::as.data.table(source_dt))
  reference_dt <- data.table::copy(data.table::as.data.table(reference_dt))
  required_source <- c(
    "SNP", "chr", "pos", "effect_allele", "other_allele", "pval"
  )
  missing_source <- setdiff(required_source, names(source_dt))
  if (length(missing_source) > 0L) {
    stop("MAGMA source data are missing columns: ",
         paste(missing_source, collapse = ", "))
  }
  required_reference <- c("SNP", "chr", "pos", "allele_key")
  missing_reference <- setdiff(required_reference, names(reference_dt))
  if (length(missing_reference) > 0L) {
    stop("MAGMA LD reference is missing columns: ",
         paste(missing_reference, collapse = ", "))
  }
  if (nrow(source_dt) == 0L) {
    return(list(
      rows = source_dt,
      audit = data.table::data.table(
        trait = character(), SNP = character(), n_source_records = integer(),
        n_reference_records = integer(), n_matching_records = integer(),
        status = character(), old_P = numeric(), new_P = numeric(),
        selected_record = character()
      )
    ))
  }

  source_dt[, `:=`(
    SNP = as.character(SNP),
    chr = as.character(chr),
    pos = as.integer(pos),
    source_row = seq_len(.N),
    source_allele_key = magma_unordered_allele_key(
      effect_allele, other_allele
    )
  )]
  reference_dt[, `:=`(
    SNP = as.character(SNP),
    chr = as.character(chr),
    pos = as.integer(pos),
    allele_key = as.character(allele_key)
  )]
  reference_dt <- unique(
    reference_dt, by = c("SNP", "chr", "pos", "allele_key")
  )

  duplicate_ids <- source_dt[, .N, by = SNP][N > 1L, SNP]
  if (length(duplicate_ids) == 0L) {
    source_dt[, c("source_row", "source_allele_key") := NULL]
    return(list(
      rows = source_dt,
      audit = data.table::data.table(
        trait = character(), SNP = character(), n_source_records = integer(),
        n_reference_records = integer(), n_matching_records = integer(),
        status = character(), old_P = numeric(), new_P = numeric(),
        selected_record = character()
      )
    ))
  }

  keep_rows <- list(source_dt[!SNP %in% duplicate_ids])
  audit_rows <- vector("list", length(duplicate_ids))

  for (i in seq_along(duplicate_ids)) {
    snp <- duplicate_ids[[i]]
    source_group <- source_dt[SNP == snp]
    reference_group <- reference_dt[SNP == snp]
    first_record <- source_group[1L]
    selected_record <- NULL
    n_matching <- 0L

    if (nrow(reference_group) == 0L) {
      selected_record <- first_record
      status <- "absent_from_LD_reference"
    } else if (nrow(reference_group) == 1L) {
      n_matching <- sum(
        source_group$chr == reference_group$chr[[1L]] &
          source_group$pos == reference_group$pos[[1L]] &
          source_group$source_allele_key == reference_group$allele_key[[1L]],
        na.rm = TRUE
      )
      matching_rows <- source_group[
        chr == reference_group$chr[[1L]] &
          pos == reference_group$pos[[1L]] &
          source_allele_key == reference_group$allele_key[[1L]]
      ]
      if (n_matching == 1L) {
        selected_record <- matching_rows[1L]
        status <- if (selected_record$source_row[[1L]] ==
                      first_record$source_row[[1L]]) {
          "retained_matching_record"
        } else {
          "replaced_first_record"
        }
      } else {
        status <- "excluded_conflicting_matched_records"
      }
    } else {
      status <- "excluded_ambiguous_reference_mapping"
    }

    if (!is.null(selected_record)) keep_rows[[length(keep_rows) + 1L]] <-
      selected_record
    audit_rows[[i]] <- data.table::data.table(
      trait = as.character(trait),
      SNP = snp,
      n_source_records = nrow(source_group),
      n_reference_records = nrow(reference_group),
      n_matching_records = as.integer(n_matching),
      status = status,
      old_P = as.numeric(first_record$pval[[1L]]),
      new_P = if (is.null(selected_record)) NA_real_ else
        as.numeric(selected_record$pval[[1L]]),
      selected_record = if (is.null(selected_record)) NA_character_ else
        magma_selected_record(selected_record)
    )
  }

  rows <- data.table::rbindlist(keep_rows, use.names = TRUE, fill = TRUE)
  data.table::setorder(rows, source_row)
  rows[, c("source_row", "source_allele_key") := NULL]
  audit <- data.table::rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
  data.table::setorder(audit, SNP)
  list(rows = rows, audit = audit)
}
