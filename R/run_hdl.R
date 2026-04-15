#' Run HDL analysis on rsid-mapped GWAS inputs
#'
#' @param gwas1_mapped_file First trait GWAS (rsid mapping format).
#' @param gwas2_mapped_file Second trait GWAS (rsid mapping format).
#' @param piecewise_path HDL reference panel path.
#' @param hdl_script Optional path to HDL source script.
#' @param gwas1_colmap Optional column mapping for trait1.
#' @param gwas2_colmap Optional column mapping for trait2.
#' @param output_dir Output directory.
#' @return Result object from HDL.
#' @export
run_hdl <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  piecewise_path,
  hdl_script = NULL,
  gwas1_colmap = list(),
  gwas2_colmap = list(),
  output_dir = "results/hdl"
) {
  ensure_dir(output_dir)

  gwas1 <- prepare_gwas_for_hdl_lcv(gwas1_mapped_file, colmap = gwas1_colmap)
  gwas2 <- prepare_gwas_for_hdl_lcv(gwas2_mapped_file, colmap = gwas2_colmap)

  if (!is.null(hdl_script)) {
    source(hdl_script)
  }

  if (!exists("HDL.rg", mode = "function")) {
    stop("未找到 HDL.rg 函数。请传入 hdl_script（包含 HDL.rg 定义）或先加载 HDL 环境。", call. = FALSE)
  }

  out <- HDL.rg(
    gwas.df1 = gwas1,
    gwas.df2 = gwas2,
    LD.path = piecewise_path
  )
  saveRDS(out, file.path(output_dir, "hdl_result.rds"))
  out
}
