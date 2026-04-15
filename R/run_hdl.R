#' Run HDL analysis
#'
#' @param gwas1_file First GWAS summary file.
#' @param gwas2_file Second GWAS summary file.
#' @param hdl_script Optional path to HDL source script.
#' @param piecewise_path HDL reference panel path.
#' @param output_dir Output directory.
#' @return Result object from HDL if available, otherwise written output path.
#' @export
run_hdl <- function(
  gwas1_file,
  gwas2_file,
  hdl_script = NULL,
  piecewise_path,
  output_dir = "results/hdl"
) {
  ensure_dir(output_dir)

  if (!is.null(hdl_script)) {
    source(hdl_script)
  }

  if (exists("HDL.rg", mode = "function")) {
    out <- HDL.rg(
      gwas.df1 = read_sumstats(gwas1_file),
      gwas.df2 = read_sumstats(gwas2_file),
      LD.path = piecewise_path
    )
    saveRDS(out, file.path(output_dir, "hdl_result.rds"))
    return(out)
  }

  stop(
    "未找到 HDL.rg 函数。请传入 hdl_script（包含 HDL.rg 定义）或先加载 HDL 环境。",
    call. = FALSE
  )
}
