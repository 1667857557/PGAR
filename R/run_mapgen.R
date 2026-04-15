#' Run MAPGEN workflow on rsid-mapped GWAS input
#'
#' @param gwas_mapped_file GWAS文件（rsid mapping格式）。
#' @param mapgen_bin MAPGEN executable path.
#' @param output_dir Output directory.
#' @param extra_args Additional command args.
#' @return Character vector of command output.
#' @export
run_Mapgen <- function(
  gwas_mapped_file,
  mapgen_bin = "mapgen",
  output_dir = "results/mapgen",
  extra_args = character(0)
) {
  ensure_dir(output_dir)

  dat <- read_mapped_gwas(gwas_mapped_file)
  validate_mapped_gwas(dat)

  mapgen_input <- file.path(output_dir, "mapgen_input.tsv")
  utils::write.table(dat, mapgen_input, sep = "\t", row.names = FALSE, quote = FALSE)

  args <- c("--input", mapgen_input, "--out", output_dir, extra_args)
  run_system_cmd(mapgen_bin, args)
}
