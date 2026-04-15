#' Run LDSC via external Python script
#'
#' @param munged_sumstats_1 LDSC-ready trait1 sumstats.gz
#' @param munged_sumstats_2 LDSC-ready trait2 sumstats.gz
#' @param ldsc_py Path to ldsc.py
#' @param ld_ref Path prefix for LD reference panel.
#' @param w_ld Path prefix for LDSC regression weights.
#' @param out_prefix Output prefix.
#' @param extra_args Extra CLI args.
#' @return Output log path.
#' @export
run_ldsc <- function(
  munged_sumstats_1,
  munged_sumstats_2,
  ldsc_py,
  ld_ref,
  w_ld,
  out_prefix = "results/ldsc/ldsc",
  extra_args = character(0)
) {
  ensure_dir(dirname(out_prefix))

  args <- c(
    ldsc_py,
    "--rg", paste(munged_sumstats_1, munged_sumstats_2, sep = ","),
    "--ref-ld-chr", ld_ref,
    "--w-ld-chr", w_ld,
    "--out", out_prefix,
    extra_args
  )

  run_system_cmd("python", args)
  paste0(out_prefix, ".log")
}
