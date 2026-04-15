#' Run CCGWAS workflow
#'
#' @param case_gwas_mapped_file Case trait mapped GWAS file.
#' @param control_gwas_mapped_file Control trait mapped GWAS file.
#' @param ccgwas_script Script with `run_ccgwas_pipeline()` or CLI entry.
#' @param case_colmap Optional case colmap.
#' @param control_colmap Optional control colmap.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_ccgwas <- function(
  case_gwas_mapped_file,
  control_gwas_mapped_file,
  ccgwas_script,
  case_colmap = list(),
  control_colmap = list(),
  output_dir = "results/ccgwas",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  case_dat <- prepare_gwas_post_input(case_gwas_mapped_file, case_colmap)
  ctrl_dat <- prepare_gwas_post_input(control_gwas_mapped_file, control_colmap)

  case_file <- file.path(output_dir, "ccgwas_case.tsv")
  ctrl_file <- file.path(output_dir, "ccgwas_control.tsv")
  utils::write.table(case_dat, case_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(ctrl_dat, ctrl_file, sep = "\t", row.names = FALSE, quote = FALSE)

  source(ccgwas_script)
  if (exists("run_ccgwas_pipeline", mode = "function")) {
    run_ccgwas_pipeline(case_gwas = case_file, control_gwas = ctrl_file, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(ccgwas_script, "--case", case_file, "--control", ctrl_file, "--out", output_dir, extra_args))
  }
  output_dir
}
