#' Run full coloc workflow from GWAS and local eQTL folder
#'
#' @param gwas_file GWAS summary statistics file.
#' @param eqtl_dir Directory containing eQTL summary files.
#' @param output_dir Output folder.
#' @param gwas_colmap Named list for GWAS column mapping.
#' @param eqtl_colmap Named list for eQTL column mapping.
#' @param gwas_type "quant" or "cc".
#' @param eqtl_type "quant" or "cc".
#' @param min_overlap Minimum overlapping SNP count to run coloc.
#' @param p1,p2,p12 Prior probabilities for coloc.abf.
#' @return data.frame of coloc summary per eQTL file.
#' @export
run_coloc <- function(
  gwas_file,
  eqtl_dir,
  output_dir = "results/coloc",
  gwas_colmap = list(),
  eqtl_colmap = list(),
  gwas_type = "quant",
  eqtl_type = "quant",
  min_overlap = 100,
  p1 = 1e-4,
  p2 = 1e-4,
  p12 = 1e-5
) {
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("Package 'coloc' is required. Install with install.packages('coloc').", call. = FALSE)
  }

  if (!dir.exists(eqtl_dir)) {
    stop("eQTL directory not found: ", eqtl_dir, call. = FALSE)
  }

  ensure_dir(output_dir)

  gwas <- read_sumstats(gwas_file, colmap = gwas_colmap)
  assert_columns(gwas, c("snp", "beta", "se", "p"), "GWAS")

  eqtl_files <- list.files(eqtl_dir, full.names = TRUE)
  eqtl_files <- eqtl_files[file.info(eqtl_files)$isdir %in% FALSE]
  if (length(eqtl_files) == 0) {
    stop("No eQTL files found in: ", eqtl_dir, call. = FALSE)
  }

  results <- vector("list", length(eqtl_files))

  for (i in seq_along(eqtl_files)) {
    eqtl_file <- eqtl_files[[i]]
    eqtl <- read_sumstats(eqtl_file, colmap = eqtl_colmap)
    assert_columns(eqtl, c("snp", "beta", "se", "p"), basename(eqtl_file))

    merged <- merge(gwas, eqtl, by = "snp", suffixes = c("_gwas", "_eqtl"))
    if (nrow(merged) < min_overlap) {
      results[[i]] <- data.frame(
        eqtl_file = basename(eqtl_file),
        nsnp = nrow(merged),
        pp_h0 = NA_real_,
        pp_h1 = NA_real_,
        pp_h2 = NA_real_,
        pp_h3 = NA_real_,
        pp_h4 = NA_real_,
        status = sprintf("skip: overlap < %d", min_overlap),
        stringsAsFactors = FALSE
      )
      next
    }

    d1 <- list(
      beta = merged$beta_gwas,
      varbeta = merged$se_gwas^2,
      snp = merged$snp,
      pvalues = merged$p_gwas,
      N = if ("n_gwas" %in% names(merged)) merged$n_gwas else NULL,
      type = gwas_type
    )

    d2 <- list(
      beta = merged$beta_eqtl,
      varbeta = merged$se_eqtl^2,
      snp = merged$snp,
      pvalues = merged$p_eqtl,
      N = if ("n_eqtl" %in% names(merged)) merged$n_eqtl else NULL,
      type = eqtl_type
    )

    fit <- coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = p1, p2 = p2, p12 = p12)
    pp <- fit$summary

    results[[i]] <- data.frame(
      eqtl_file = basename(eqtl_file),
      nsnp = nrow(merged),
      pp_h0 = as.numeric(pp["PP.H0.abf"]),
      pp_h1 = as.numeric(pp["PP.H1.abf"]),
      pp_h2 = as.numeric(pp["PP.H2.abf"]),
      pp_h3 = as.numeric(pp["PP.H3.abf"]),
      pp_h4 = as.numeric(pp["PP.H4.abf"]),
      status = "ok",
      stringsAsFactors = FALSE
    )

    utils::write.table(
      fit$results,
      file = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(eqtl_file)), "_snp_posterior.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  out <- do.call(rbind, results)
  utils::write.table(out, file.path(output_dir, "coloc_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out
}
