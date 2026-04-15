#' Read GWAS data processed by rsid mapping workflow
#'
#' 标准输入应来自 `rsid mapping.R` 的输出。
#' 期望列名（不区分大小写）:
#' rsid, chr, pos, effect_allele, other_allele, beta, se, pval, n
#'
#' @param file GWAS文件路径。
#' @param sep 分隔符，默认根据扩展名自动识别。
#' @param colmap 可选的列名映射（命名方式: 标准名=原始名）。
#' @return 标准化后的 data.frame。
#' @export
read_mapped_gwas <- function(file, sep = NULL, colmap = list()) {
  dat <- read_sumstats(file = file, sep = sep, colmap = colmap)

  alias_map <- c(
    snp = "rsid",
    p = "pval",
    a1 = "effect_allele",
    ea = "effect_allele",
    alt = "effect_allele",
    a2 = "other_allele",
    nea = "other_allele",
    ref = "other_allele"
  )

  for (from in names(alias_map)) {
    to <- alias_map[[from]]
    if (from %in% names(dat) && !to %in% names(dat)) {
      names(dat)[names(dat) == from] <- to
    }
  }

  required <- c("rsid", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval", "n")
  assert_columns(dat, required, "rsid mapping GWAS")

  if (!"z" %in% names(dat)) {
    dat$z <- dat$beta / dat$se
  }

  dat
}

#' Validate mapped GWAS format
#'
#' @param dat data.frame returned by `read_mapped_gwas()`.
#' @return logical TRUE (invisible) if valid, otherwise error.
#' @export
validate_mapped_gwas <- function(dat) {
  required <- c("rsid", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval", "n", "z")
  assert_columns(dat, required, "mapped GWAS")

  if (anyNA(dat$rsid) || any(dat$rsid == "")) {
    stop("mapped GWAS: rsid 存在缺失或空值", call. = FALSE)
  }
  if (anyNA(dat$beta) || anyNA(dat$se) || anyNA(dat$pval)) {
    stop("mapped GWAS: beta/se/pval 存在缺失", call. = FALSE)
  }
  invisible(TRUE)
}
