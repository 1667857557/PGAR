#' Run MAPGEN workflow on rsid-mapped GWAS input
#'
#' @param gwas_mapped_file GWAS文件（rsid mapping格式）。
#' @param gwas_colmap Optional column mapping for GWAS input.
#' @param mapgen_bin MAPGEN executable path.
#' @param output_dir Output directory.
#' @param extra_args Additional command args.
#' @return Character vector of command output.
#' @export
run_Mapgen <- function(
  gwas_mapped_file,
  mapgen_bin = "mapgen",
  output_dir = "results/mapgen",
  gwas_colmap = list(),
  extra_args = character(0)
) {
  ensure_dir(output_dir)

  dat <- prepare_gwas_for_mapgen(gwas_mapped_file, colmap = gwas_colmap)

  mapgen_input <- file.path(output_dir, "mapgen_input.tsv")
  utils::write.table(dat, mapgen_input, sep = "\t", row.names = FALSE, quote = FALSE)

  args <- c("--input", mapgen_input, "--out", output_dir, extra_args)
  run_system_cmd(mapgen_bin, args)
}
