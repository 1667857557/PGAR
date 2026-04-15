#' Run MAPGEN workflow using command-line executable
#'
#' @param input_file Input file for MAPGEN.
#' @param mapgen_bin MAPGEN executable path.
#' @param output_dir Output directory.
#' @param extra_args Additional command args.
#' @return Character vector of command output.
#' @export
run_Mapgen <- function(
  input_file,
  mapgen_bin = "mapgen",
  output_dir = "results/mapgen",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  args <- c("--input", input_file, "--out", output_dir, extra_args)
  run_system_cmd(mapgen_bin, args)
}
