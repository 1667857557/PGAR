#' Run clumping and fine-mapping workflow
#'
#' @param gwas_mapped_file Mapped GWAS file.
#' @param plink_bin Path to PLINK executable.
#' @param ref_bfile Reference panel PLINK bfile prefix.
#' @param output_dir Output directory.
#' @param clump_p1 Clump p1 threshold.
#' @param clump_r2 Clump r2 threshold.
#' @param clump_kb Clump kb window.
#' @param finemap_cmd Optional fine-mapping executable.
#' @param gwas_colmap Optional column mapping.
#' @param extra_args Extra args for fine-mapping command.
#' @return List with clump and fine-mapping outputs.
#' @export
run_clump_fine_mapping <- function(
  gwas_mapped_file,
  plink_bin,
  ref_bfile,
  output_dir = "results/clump_fine_mapping",
  clump_p1 = 5e-8,
  clump_r2 = 0.1,
  clump_kb = 250,
  finemap_cmd = NULL,
  gwas_colmap = list(),
  extra_args = character(0)
) {
  ensure_dir(output_dir)

  dat <- prepare_gwas_for_ldsc(gwas_mapped_file, colmap = gwas_colmap)
  dat$P <- as.numeric(dat$P)
  dat$N <- as.numeric(dat$N)

  clump_input <- file.path(output_dir, "clump_input.tsv")
  utils::write.table(dat[, c("SNP", "P")], clump_input, sep = "\t", row.names = FALSE, quote = FALSE)

  clump_prefix <- file.path(output_dir, "plink_clump")
  run_system_cmd(plink_bin, c(
    "--bfile", ref_bfile,
    "--clump", clump_input,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", as.character(clump_p1),
    "--clump-r2", as.character(clump_r2),
    "--clump-kb", as.character(clump_kb),
    "--out", clump_prefix
  ))

  fine_out <- NULL
  if (!is.null(finemap_cmd)) {
    fine_input <- file.path(output_dir, "fine_mapping_input.tsv")
    utils::write.table(dat, fine_input, sep = "\t", row.names = FALSE, quote = FALSE)
    fine_out <- file.path(output_dir, "fine_mapping")
    run_system_cmd(finemap_cmd, c("--sumstats", fine_input, "--out", fine_out, extra_args))
  }

  list(
    clump = paste0(clump_prefix, ".clumped"),
    fine_mapping = fine_out
  )
}
