assert_columns <- function(dat, required, label) {
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("%s 缺少必需列: %s", label, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

run_system_cmd <- function(command, args, workdir = NULL) {
  out <- system2(command = command, args = args, stdout = TRUE, stderr = TRUE, wd = workdir)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste0(
        "External command failed (", command, ") with status ", status, "\n",
        paste(out, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  out
}
