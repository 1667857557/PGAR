#' Run scDRS workflow
#'
#' @param gwas_mapped_file Mapped GWAS file.
#' @param scdrs_cmd scDRS executable (e.g. python -m scdrs).
#' @param output_dir Output directory.
#' @param gwas_colmap Optional colmap.
#' @param python_bin Python executable.
#' @param extra_args Extra args passed to scDRS.
#' @return Output directory.
#' @export
run_scdrs <- function(
  gwas_mapped_file,
  scdrs_cmd = "scdrs",
  output_dir = "results/scdrs",
  gwas_colmap = list(),
  python_bin = "python",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  dat <- prepare_gwas_post_input(gwas_mapped_file, colmap = gwas_colmap)

  scdrs_input <- file.path(output_dir, "scdrs_gwas.tsv")
  utils::write.table(dat, scdrs_input, sep = "\t", row.names = FALSE, quote = FALSE)

  if (grepl("\\s", scdrs_cmd)) {
    cmd <- strsplit(scdrs_cmd, "\\s+")[[1]]
    run_system_cmd(cmd[1], c(cmd[-1], "--gwas", scdrs_input, "--out", output_dir, extra_args))
  } else if (scdrs_cmd == "scdrs") {
    run_system_cmd(python_bin, c("-m", "scdrs", "--gwas", scdrs_input, "--out", output_dir, extra_args))
  } else {
    run_system_cmd(scdrs_cmd, c("--gwas", scdrs_input, "--out", output_dir, extra_args))
  }

  output_dir
}
