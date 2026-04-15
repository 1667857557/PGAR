#' Run LCV analysis on rsid-mapped GWAS inputs
#'
#' @param gwas1_mapped_file First trait GWAS (rsid mapping format).
#' @param gwas2_mapped_file Second trait GWAS (rsid mapping format).
#' @param gwas1_colmap Optional column mapping for trait1.
#' @param gwas2_colmap Optional column mapping for trait2.
#' @param lcv_script Optional path to script that defines RunLCV.
#' @param output_dir Output directory.
#' @return LCV result object.
#' @export
run_lcv <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  lcv_script = NULL,
  output_dir = "results/lcv",
  gwas1_colmap = list(),
  gwas2_colmap = list()
) {
  ensure_dir(output_dir)

  gwas1 <- prepare_gwas_for_hdl_lcv(gwas1_mapped_file, colmap = gwas1_colmap)
  gwas2 <- prepare_gwas_for_hdl_lcv(gwas2_mapped_file, colmap = gwas2_colmap)

  if (!is.null(lcv_script)) {
    source(lcv_script)
  }

  if (exists("RunLCV", mode = "function")) {
    res <- RunLCV(trait1 = gwas1, trait2 = gwas2)
    saveRDS(res, file.path(output_dir, "lcv_result.rds"))
    return(res)
  }

  if (requireNamespace("lcv", quietly = TRUE) && exists("RunLCV", where = asNamespace("lcv"), mode = "function")) {
    runlcv <- get("RunLCV", envir = asNamespace("lcv"))
    res <- runlcv(trait1 = gwas1, trait2 = gwas2)
    saveRDS(res, file.path(output_dir, "lcv_result.rds"))
    return(res)
  }

  stop("未找到 RunLCV 函数。请加载 LCV 脚本或安装并加载对应 R 包。", call. = FALSE)
}
