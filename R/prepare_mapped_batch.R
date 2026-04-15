#' Prepare multiple mapped GWAS files for batch workflows
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param colmaps Optional list of column maps.
#' @param output_dir Output directory.
#' @param prefix Output prefix.
#' @return Vector of normalized GWAS file paths.
#' @export
prepare_mapped_batch <- function(gwas_mapped_files, colmaps = NULL, output_dir, prefix) {
  if (length(gwas_mapped_files) == 0) {
    stop("gwas_mapped_files 不能为空", call. = FALSE)
  }
  if (is.null(colmaps)) {
    colmaps <- replicate(length(gwas_mapped_files), list(), simplify = FALSE)
  }
  if (length(colmaps) != length(gwas_mapped_files)) {
    stop("colmaps 长度必须与 gwas_mapped_files 相同", call. = FALSE)
  }

  out_files <- character(length(gwas_mapped_files))
  for (i in seq_along(gwas_mapped_files)) {
    dat <- prepare_gwas_post_input(gwas_mapped_files[[i]], colmap = colmaps[[i]])
    out_file <- file.path(output_dir, sprintf("%s_trait_%02d.tsv", prefix, i))
    utils::write.table(dat, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    out_files[[i]] <- out_file
  }
  out_files
}
