#' Standardize mapped GWAS for post-analysis workflows
#'
#' @param gwas_mapped_file File from rsid mapping workflow.
#' @param colmap Optional column mapping.
#' @return Standard validated mapped GWAS data.frame.
#' @export
prepare_gwas_post_input <- function(gwas_mapped_file, colmap = list()) {
  dat <- read_mapped_gwas(gwas_mapped_file, colmap = colmap)
  validate_mapped_gwas(dat)
  dat
}

#' Prepare mapped GWAS for coloc
#'
#' @param gwas_mapped_file File from rsid mapping workflow.
#' @param colmap Optional column mapping.
#' @return data.frame with columns snp, beta, se, p, n.
#' @export
prepare_gwas_for_coloc <- function(gwas_mapped_file, colmap = list()) {
  dat <- prepare_gwas_post_input(gwas_mapped_file, colmap = colmap)
  out <- dat[, c("rsid", "beta", "se", "pval", "n")]
  names(out) <- c("snp", "beta", "se", "p", "n")
  out
}

#' Prepare mapped GWAS for LDSC
#'
#' @param gwas_mapped_file File from rsid mapping workflow.
#' @param colmap Optional column mapping.
#' @return data.frame with columns SNP, A1, A2, P, N, Z.
#' @export
prepare_gwas_for_ldsc <- function(gwas_mapped_file, colmap = list()) {
  dat <- prepare_gwas_post_input(gwas_mapped_file, colmap = colmap)
  out <- dat[, c("rsid", "effect_allele", "other_allele", "pval", "n", "z")]
  names(out) <- c("SNP", "A1", "A2", "P", "N", "Z")
  out
}

#' Prepare mapped GWAS for HDL/LCV
#'
#' @param gwas_mapped_file File from rsid mapping workflow.
#' @param colmap Optional column mapping.
#' @return validated mapped GWAS data.frame.
#' @export
prepare_gwas_for_hdl_lcv <- function(gwas_mapped_file, colmap = list()) {
  prepare_gwas_post_input(gwas_mapped_file, colmap = colmap)
}

#' Prepare mapped GWAS for MAPGEN
#'
#' @param gwas_mapped_file File from rsid mapping workflow.
#' @param colmap Optional column mapping.
#' @return validated mapped GWAS data.frame.
#' @export
prepare_gwas_for_mapgen <- function(gwas_mapped_file, colmap = list()) {
  prepare_gwas_post_input(gwas_mapped_file, colmap = colmap)
}
