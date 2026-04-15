#' Run FLAMES workflow on rsid-mapped GWAS input
#'
#' @param gwas_mapped_file GWAS file in rsid mapping format.
#' @param flames_script R script path containing `run_flames_pipeline()` or CLI entry.
#' @param gwas_colmap Optional column mapping.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args when script is executed by Rscript.
#' @return Path to FLAMES output directory.
#' @export
run_flames <- function(
  gwas_mapped_file,
  flames_script,
  gwas_colmap = list(),
  output_dir = "results/flames",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  gwas <- prepare_gwas_post_input(gwas_mapped_file, colmap = gwas_colmap)
  gwas_input <- file.path(output_dir, "flames_input.tsv")
  utils::write.table(gwas, gwas_input, sep = "\t", row.names = FALSE, quote = FALSE)

  source(flames_script)
  if (exists("run_flames_pipeline", mode = "function")) {
    run_flames_pipeline(gwas_file = gwas_input, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(flames_script, "--gwas", gwas_input, "--out", output_dir, extra_args))
  }

  output_dir
}
