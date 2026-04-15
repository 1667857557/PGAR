#' Read and normalize GWAS summary statistics
#'
#' @param file Path to summary statistics file.
#' @param sep Field separator. If NULL, guessed from extension.
#' @param colmap Named list used to map user columns to standard names.
#'   Supported standard names: SNP, chr, pos, effect_allele, other_allele,
#'   beta, se, p, n, eaf, z.
#' @return A data.frame with standardized lower-case columns.
#' @export
read_sumstats <- function(file, sep = NULL, colmap = list()) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }

  if (is.null(sep)) {
    ext <- tolower(tools::file_ext(file))
    sep <- if (ext %in% c("tsv", "txt", "bgz", "gz")) "\t" else ","
  }

  dat <- utils::read.table(
    file = file,
    sep = sep,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )

  if (length(colmap) > 0) {
    for (std_name in names(colmap)) {
      raw_name <- colmap[[std_name]]
      if (!raw_name %in% names(dat)) {
        stop("Column mapping failed. Column not found: ", raw_name, call. = FALSE)
      }
      names(dat)[names(dat) == raw_name] <- std_name
    }
  }

  names(dat) <- tolower(names(dat))

  if (!"snp" %in% names(dat)) {
    stop("Input summary statistics must include SNP column (or map one to SNP).", call. = FALSE)
  }

  if (!"z" %in% names(dat) && all(c("beta", "se") %in% names(dat))) {
    dat$z <- dat$beta / dat$se
  }

  dat
}
