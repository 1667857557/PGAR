#' Run HDL analysis on rsid-mapped GWAS inputs
#'
#' @param gwas1_mapped_file First trait GWAS (rsid mapping format).
#' @param gwas2_mapped_file Second trait GWAS (rsid mapping format).
#' @param hdl_script Optional path to HDL source script.
#' @param piecewise_path HDL reference panel path.
#' @param output_dir Output directory.
#' @return Result object from HDL.
#' @export
run_hdl <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  hdl_script = NULL,
  piecewise_path,
  output_dir = "results/hdl"
) {
  ensure_dir(output_dir)

  gwas1 <- read_mapped_gwas(gwas1_mapped_file)
  gwas2 <- read_mapped_gwas(gwas2_mapped_file)
  validate_mapped_gwas(gwas1)
  validate_mapped_gwas(gwas2)

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
