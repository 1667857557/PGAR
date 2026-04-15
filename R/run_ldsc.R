#' Run LDSC on rsid-mapped GWAS input files
#'
#' @param gwas1_mapped_file First trait GWAS (rsid mapping format).
#' @param gwas2_mapped_file Second trait GWAS (rsid mapping format).
#' @param munge_py Path to munge_sumstats.py.
#' @param ldsc_py Path to ldsc.py
#' @param snplist Path to HapMap3 snplist used by munge.
#' @param ld_ref Path prefix for LD reference panel.
#' @param w_ld Path prefix for LDSC regression weights.
#' @param gwas1_colmap Optional column mapping for trait1.
#' @param gwas2_colmap Optional column mapping for trait2.
#' @param out_prefix Output prefix.
#' @param python_bin Python executable.
#' @param extra_args Extra CLI args for ldsc.
#' @return Output log path.
#' @export
run_ldsc <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  munge_py,
  ldsc_py,
  snplist,
  ld_ref,
  w_ld,
  gwas1_colmap = list(),
  gwas2_colmap = list(),
  out_prefix = "results/ldsc/ldsc",
  python_bin = "python",
  extra_args = character(0)
) {
  out_dir <- dirname(out_prefix)
  ensure_dir(out_dir)

  x1 <- prepare_gwas_for_ldsc(gwas1_mapped_file, colmap = gwas1_colmap)
  x2 <- prepare_gwas_for_ldsc(gwas2_mapped_file, colmap = gwas2_colmap)

  prep_for_ldsc <- function(x, prefix) {
    raw <- file.path(out_dir, paste0(prefix, "_raw.tsv"))
    utils::write.table(x, raw, sep = "\t", row.names = FALSE, quote = FALSE)
    raw
  }

  raw1 <- prep_for_ldsc(x1, "trait1")
  raw2 <- prep_for_ldsc(x2, "trait2")

  munge_out1 <- file.path(out_dir, "trait1")
  munge_out2 <- file.path(out_dir, "trait2")

  run_system_cmd(python_bin, c(
    munge_py,
    "--sumstats", raw1,
    "--out", munge_out1,
    "--merge-alleles", snplist,
    "--N-col", "N",
    "--a1", "A1",
    "--a2", "A2",
    "--snp", "SNP",
    "--p", "P"
  ))

  run_system_cmd(python_bin, c(
    munge_py,
    "--sumstats", raw2,
    "--out", munge_out2,
    "--merge-alleles", snplist,
    "--N-col", "N",
    "--a1", "A1",
    "--a2", "A2",
    "--snp", "SNP",
    "--p", "P"
  ))

  run_system_cmd(python_bin, c(
    ldsc_py,
    "--rg", paste0(munge_out1, ".sumstats.gz,", munge_out2, ".sumstats.gz"),
    "--ref-ld-chr", ld_ref,
    "--w-ld-chr", w_ld,
    "--out", out_prefix,
    extra_args
  ))

  paste0(out_prefix, ".log")
}
