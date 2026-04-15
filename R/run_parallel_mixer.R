#' Run parallel MiXeR workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param mixer_script Script with `run_parallel_mixer_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param cores Number of cores.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_parallel_mixer <- function(
  gwas_mapped_files,
  mixer_script,
  colmaps = NULL,
  cores = 4,
  output_dir = "results/parallel_mixer",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "mixer")

  source(mixer_script)
  if (exists("run_parallel_mixer_pipeline", mode = "function")) {
    run_parallel_mixer_pipeline(gwas_files = prepped, out_dir = output_dir, cores = cores)
  } else {
    run_system_cmd("Rscript", c(mixer_script, "--gwas-list", paste(prepped, collapse = ","), "--cores", as.character(cores), "--out", output_dir, extra_args))
  }
  output_dir
}
