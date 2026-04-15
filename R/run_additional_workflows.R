#' Run FLAMES workflow on rsid-mapped GWAS input
#'
#' @param gwas_mapped_file GWAS file in rsid mapping format.
#' @param flames_script R script path containing `run_flames_pipeline()` or CLI entry.
#' @param gwas_colmap Optional column mapping.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args when script is executed by Rscript.
#' @return Path to FLAMES output directory.
#' @export
run_flames <- function(
  gwas_mapped_file,
  flames_script,
  gwas_colmap = list(),
  output_dir = "results/flames",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  gwas <- prepare_gwas_post_input(gwas_mapped_file, colmap = gwas_colmap)
  gwas_input <- file.path(output_dir, "flames_input.tsv")
  utils::write.table(gwas, gwas_input, sep = "\t", row.names = FALSE, quote = FALSE)

  source(flames_script)
  if (exists("run_flames_pipeline", mode = "function")) {
    run_flames_pipeline(gwas_file = gwas_input, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(flames_script, "--gwas", gwas_input, "--out", output_dir, extra_args))
  }

  output_dir
}

#' Run T-SEM workflow
#'
#' @param gwas_mapped_files Vector of rsid-mapped GWAS file paths.
#' @param t_sem_script Script with `run_t_sem_pipeline()` or CLI entry.
#' @param colmaps Optional list of column mappings aligned to files.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_t_sem <- function(
  gwas_mapped_files,
  t_sem_script,
  colmaps = NULL,
  output_dir = "results/t_sem",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- .prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "tsem")

  source(t_sem_script)
  if (exists("run_t_sem_pipeline", mode = "function")) {
    run_t_sem_pipeline(gwas_files = prepped, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(t_sem_script, "--gwas-list", paste(prepped, collapse = ","), "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run P-SEM workflow
#'
#' @param gwas_mapped_files Vector of rsid-mapped GWAS file paths.
#' @param p_sem_script Script with `run_p_sem_pipeline()` or CLI entry.
#' @param colmaps Optional list of column mappings aligned to files.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_p_sem <- function(
  gwas_mapped_files,
  p_sem_script,
  colmaps = NULL,
  output_dir = "results/p_sem",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- .prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "psem")

  source(p_sem_script)
  if (exists("run_p_sem_pipeline", mode = "function")) {
    run_p_sem_pipeline(gwas_files = prepped, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(p_sem_script, "--gwas-list", paste(prepped, collapse = ","), "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run GNA workflow
#'
#' @param gwas1_mapped_file Trait1 mapped GWAS file.
#' @param gwas2_mapped_file Trait2 mapped GWAS file.
#' @param gna_script Script with `run_gna_pipeline()` or CLI entry.
#' @param gwas1_colmap Optional trait1 colmap.
#' @param gwas2_colmap Optional trait2 colmap.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_gna <- function(
  gwas1_mapped_file,
  gwas2_mapped_file,
  gna_script,
  gwas1_colmap = list(),
  gwas2_colmap = list(),
  output_dir = "results/gna",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  g1 <- prepare_gwas_post_input(gwas1_mapped_file, gwas1_colmap)
  g2 <- prepare_gwas_post_input(gwas2_mapped_file, gwas2_colmap)
  f1 <- file.path(output_dir, "gna_trait1.tsv")
  f2 <- file.path(output_dir, "gna_trait2.tsv")
  utils::write.table(g1, f1, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(g2, f2, sep = "\t", row.names = FALSE, quote = FALSE)

  source(gna_script)
  if (exists("run_gna_pipeline", mode = "function")) {
    run_gna_pipeline(trait1 = f1, trait2 = f2, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(gna_script, "--gwas1", f1, "--gwas2", f2, "--out", output_dir, extra_args))
  }

  output_dir
}

#' Run parallel LAVA workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param lava_script Script with `run_parallel_lava_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param cores Number of cores.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_parallel_lava <- function(
  gwas_mapped_files,
  lava_script,
  colmaps = NULL,
  cores = 4,
  output_dir = "results/parallel_lava",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- .prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "lava")

  source(lava_script)
  if (exists("run_parallel_lava_pipeline", mode = "function")) {
    run_parallel_lava_pipeline(gwas_files = prepped, out_dir = output_dir, cores = cores)
  } else {
    run_system_cmd("Rscript", c(lava_script, "--gwas-list", paste(prepped, collapse = ","), "--cores", as.character(cores), "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run parallel MiXeR workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param mixer_script Script with `run_parallel_mixer_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param cores Number of cores.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_parallel_mixer <- function(
  gwas_mapped_files,
  mixer_script,
  colmaps = NULL,
  cores = 4,
  output_dir = "results/parallel_mixer",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- .prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "mixer")

  source(mixer_script)
  if (exists("run_parallel_mixer_pipeline", mode = "function")) {
    run_parallel_mixer_pipeline(gwas_files = prepped, out_dir = output_dir, cores = cores)
  } else {
    run_system_cmd("Rscript", c(mixer_script, "--gwas-list", paste(prepped, collapse = ","), "--cores", as.character(cores), "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run CCGWAS workflow
#'
#' @param case_gwas_mapped_file Case trait mapped GWAS file.
#' @param control_gwas_mapped_file Control trait mapped GWAS file.
#' @param ccgwas_script Script with `run_ccgwas_pipeline()` or CLI entry.
#' @param case_colmap Optional case colmap.
#' @param control_colmap Optional control colmap.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_ccgwas <- function(
  case_gwas_mapped_file,
  control_gwas_mapped_file,
  ccgwas_script,
  case_colmap = list(),
  control_colmap = list(),
  output_dir = "results/ccgwas",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  case_dat <- prepare_gwas_post_input(case_gwas_mapped_file, case_colmap)
  ctrl_dat <- prepare_gwas_post_input(control_gwas_mapped_file, control_colmap)

  case_file <- file.path(output_dir, "ccgwas_case.tsv")
  ctrl_file <- file.path(output_dir, "ccgwas_control.tsv")
  utils::write.table(case_dat, case_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(ctrl_dat, ctrl_file, sep = "\t", row.names = FALSE, quote = FALSE)

  source(ccgwas_script)
  if (exists("run_ccgwas_pipeline", mode = "function")) {
    run_ccgwas_pipeline(case_gwas = case_file, control_gwas = ctrl_file, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(ccgwas_script, "--case", case_file, "--control", ctrl_file, "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run multivariate GWAS workflow
#'
#' @param gwas_mapped_files Vector of mapped GWAS files.
#' @param multivariate_script Script with `run_multivariate_gwas_pipeline()` or CLI entry.
#' @param colmaps Optional list of colmaps.
#' @param output_dir Output directory.
#' @param extra_args Extra CLI args.
#' @return Output directory.
#' @export
run_multivariate_gwas <- function(
  gwas_mapped_files,
  multivariate_script,
  colmaps = NULL,
  output_dir = "results/multivariate_gwas",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  prepped <- .prepare_mapped_batch(gwas_mapped_files, colmaps = colmaps, output_dir = output_dir, prefix = "mvgwas")

  source(multivariate_script)
  if (exists("run_multivariate_gwas_pipeline", mode = "function")) {
    run_multivariate_gwas_pipeline(gwas_files = prepped, out_dir = output_dir)
  } else {
    run_system_cmd("Rscript", c(multivariate_script, "--gwas-list", paste(prepped, collapse = ","), "--out", output_dir, extra_args))
  }
  output_dir
}

#' Run clumping and fine-mapping workflow
#'
#' @param gwas_mapped_file Mapped GWAS file.
#' @param plink_bin Path to PLINK executable.
#' @param ref_bfile Reference panel PLINK bfile prefix.
#' @param output_dir Output directory.
#' @param clump_p1 Clump p1 threshold.
#' @param clump_r2 Clump r2 threshold.
#' @param clump_kb Clump kb window.
#' @param finemap_cmd Optional fine-mapping executable.
#' @param gwas_colmap Optional column mapping.
#' @param extra_args Extra args for fine-mapping command.
#' @return List with clump and fine-mapping outputs.
#' @export
run_clump_fine_mapping <- function(
  gwas_mapped_file,
  plink_bin,
  ref_bfile,
  output_dir = "results/clump_fine_mapping",
  clump_p1 = 5e-8,
  clump_r2 = 0.1,
  clump_kb = 250,
  finemap_cmd = NULL,
  gwas_colmap = list(),
  extra_args = character(0)
) {
  ensure_dir(output_dir)

  dat <- prepare_gwas_for_ldsc(gwas_mapped_file, colmap = gwas_colmap)
  dat$P <- as.numeric(dat$P)
  dat$N <- as.numeric(dat$N)

  clump_input <- file.path(output_dir, "clump_input.tsv")
  utils::write.table(dat[, c("SNP", "P")], clump_input, sep = "\t", row.names = FALSE, quote = FALSE)

  clump_prefix <- file.path(output_dir, "plink_clump")
  run_system_cmd(plink_bin, c(
    "--bfile", ref_bfile,
    "--clump", clump_input,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", as.character(clump_p1),
    "--clump-r2", as.character(clump_r2),
    "--clump-kb", as.character(clump_kb),
    "--out", clump_prefix
  ))

  fine_out <- NULL
  if (!is.null(finemap_cmd)) {
    fine_input <- file.path(output_dir, "fine_mapping_input.tsv")
    utils::write.table(dat, fine_input, sep = "\t", row.names = FALSE, quote = FALSE)
    fine_out <- file.path(output_dir, "fine_mapping")
    run_system_cmd(finemap_cmd, c("--sumstats", fine_input, "--out", fine_out, extra_args))
  }

  list(
    clump = paste0(clump_prefix, ".clumped"),
    fine_mapping = fine_out
  )
}

#' Run scDRS workflow
#'
#' @param gwas_mapped_file Mapped GWAS file.
#' @param scdrs_cmd scDRS executable (e.g. python -m scdrs).
#' @param output_dir Output directory.
#' @param gwas_colmap Optional colmap.
#' @param python_bin Python executable.
#' @param extra_args Extra args passed to scDRS.
#' @return Output directory.
#' @export
run_scdrs <- function(
  gwas_mapped_file,
  scdrs_cmd = "scdrs",
  output_dir = "results/scdrs",
  gwas_colmap = list(),
  python_bin = "python",
  extra_args = character(0)
) {
  ensure_dir(output_dir)
  dat <- prepare_gwas_post_input(gwas_mapped_file, colmap = gwas_colmap)

  scdrs_input <- file.path(output_dir, "scdrs_gwas.tsv")
  utils::write.table(dat, scdrs_input, sep = "\t", row.names = FALSE, quote = FALSE)

  if (grepl("\\s", scdrs_cmd)) {
    cmd <- strsplit(scdrs_cmd, "\\s+")[[1]]
    run_system_cmd(cmd[1], c(cmd[-1], "--gwas", scdrs_input, "--out", output_dir, extra_args))
  } else if (scdrs_cmd == "scdrs") {
    run_system_cmd(python_bin, c("-m", "scdrs", "--gwas", scdrs_input, "--out", output_dir, extra_args))
  } else {
    run_system_cmd(scdrs_cmd, c("--gwas", scdrs_input, "--out", output_dir, extra_args))
  }

  output_dir
}

.prepare_mapped_batch <- function(gwas_mapped_files, colmaps = NULL, output_dir, prefix) {
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
