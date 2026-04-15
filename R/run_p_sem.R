#' Run P-SEM workflow
#'
#' @param gwas_mapped_files Vector of rsid-mapped GWAS file paths.
#' @param p_sem_script Script with `run_p_sem_pipeline()` or CLI entry.
#' @param colmaps Optional list of column mappings aligned to files.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_p_sem <- function(
  gwas_mapped_files,
  p_sem_script,
  colmaps = NULL,
  output_dir = "results/p_sem",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "psem")

  source(p_sem_script)
  if (exists("run_p_sem_pipeline", mode = "function")) {
    run_p_sem_pipeline(gwas_files = prepped, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(p_sem_script, "--gwas-list", paste(prepped, collapse = ","), "--out", output_dir, extra_args))
  }
  output_dir
}
