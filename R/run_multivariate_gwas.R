#' Run multivariate GWAS workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param multivariate_script Script with `run_multivariate_gwas_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_multivariate_gwas <- function(
  gwas_mapped_files,
  multivariate_script,
  colmaps = NULL,
  output_dir = "results/multivariate_gwas",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "mvgwas")

  source(multivariate_script)
  if (exists("run_multivariate_gwas_pipeline", mode = "function")) {
    run_multivariate_gwas_pipeline(gwas_files = prepped, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(multivariate_script, "--gwas-list", paste(prepped, collapse = ","), "--out", output_dir, extra_args))
  }
  output_dir
}
