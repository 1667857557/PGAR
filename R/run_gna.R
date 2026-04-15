#' Run GNA workflow
#'
#' @param gwas1_mapped_file Trait1 mapped GWAS file.
#' @param gwas2_mapped_file Trait2 mapped GWAS file.
#' @param gna_script Script with `run_gna_pipeline()` or CLI entry.
#' @param gwas1_colmap Optional trait1 colmap.
#' @param gwas2_colmap Optional trait2 colmap.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_gna <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  gna_script,
  gwas1_colmap = list(),
  gwas2_colmap = list(),
  output_dir = "results/gna",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  g1 <- prepare_gwas_post_input(gwas1_mapped_file, gwas1_colmap)
  g2 <- prepare_gwas_post_input(gwas2_mapped_file, gwas2_colmap)
  f1 <- file.path(output_dir, "gna_trait1.tsv")
  f2 <- file.path(output_dir, "gna_trait2.tsv")
  utils::write.table(g1, f1, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(g2, f2, sep = "\t", row.names = FALSE, quote = FALSE)

  source(gna_script)
  if (exists("run_gna_pipeline", mode = "function")) {
    run_gna_pipeline(trait1 = f1, trait2 = f2, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(gna_script, "--gwas1", f1, "--gwas2", f2, "--out", output_dir, extra_args))
  }

  output_dir
}
