#' Run parallel LAVA workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param lava_script Script with `run_parallel_lava_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param cores Number of cores.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_parallel_lava <- function(
  gwas_mapped_files,
  lava_script,
  colmaps = NULL,
  cores = 4,
  output_dir = "results/parallel_lava",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "lava")

  source(lava_script)
  if (exists("run_parallel_lava_pipeline", mode = "function")) {
    run_parallel_lava_pipeline(gwas_files = prepped, out_dir = output_dir, cores = cores)
  } else {
    run_system_cmd("Rscript", c(lava_script, "--gwas-list", paste(prepped, collapse = ","), "--cores", as.character(cores), "--out", output_dir, extra_args))
  }
  output_dir
}
