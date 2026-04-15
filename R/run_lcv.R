#' Run LCV analysis
#'
#' @param gwas1_file First GWAS summary file.
#' @param gwas2_file Second GWAS summary file.
#' @param lcv_script Optional path to script that defines RunLCV.
#' @param output_dir Output directory.
#' @return LCV result object.
#' @export
run_lcv <- function(
  gwas1_file,
  gwas2_file,
  lcv_script = NULL,
  output_dir = "results/lcv"
) {
  ensure_dir(output_dir)

  if (!is.null(lcv_script)) {
    source(lcv_script)
  }

  if (exists("RunLCV", mode = "function")) {
    res <- RunLCV(
      trait1 = read_sumstats(gwas1_file),
      trait2 = read_sumstats(gwas2_file)
    )
    saveRDS(res, file.path(output_dir, "lcv_result.rds"))
    return(res)
  }

  if (requireNamespace("lcv", quietly = TRUE) && exists("RunLCV", where = asNamespace("lcv"), mode = "function")) {
    runlcv <- get("RunLCV", envir = asNamespace("lcv"))
    res <- runlcv(
      trait1 = read_sumstats(gwas1_file),
      trait2 = read_sumstats(gwas2_file)
    )
    saveRDS(res, file.path(output_dir, "lcv_result.rds"))
    return(res)
  }

  stop("未找到 RunLCV 函数。请加载 LCV 脚本或安装并加载对应 R 包。", call. = FALSE)
}
