# PGAR

一个用于 **GWAS-post** 分析的 R 包脚手架，核心约定是：

> 所有 post 分析函数都使用 **`rsid mapping.R` 处理后的 GWAS 数据格式** 作为输入。

## 1) 统一的 GWAS 输入格式（rsid mapping 输出）

列名（不区分大小写）应可映射到：

- `rsid`
- `chr`
- `pos`
- `effect_allele`
- `other_allele`
- `beta`
- `se`
- `pval`
- `n`

包内 `read_mapped_gwas()` / `validate_mapped_gwas()` 会进行校验。

## 2) 主要函数

> 代码组织：每个流程函数放在独立的 `R/run_*.R` 文件中，便于按功能维护。

### 2.1 预处理函数（对应你上传代码的各步骤封装）

- `prepare_gwas_post_input()`
- `prepare_gwas_for_coloc()`
- `prepare_gwas_for_ldsc()`
- `prepare_gwas_for_hdl_lcv()`
- `prepare_gwas_for_mapgen()`

### 2.2 分析函数

- `run_coloc(gwas_mapped_file, eqtl_dir, ...)`
- `run_hdl(gwas1_mapped_file, gwas2_mapped_file, ...)`
- `run_ldsc(gwas1_mapped_file, gwas2_mapped_file, ...)`
- `run_Mapgen(gwas_mapped_file, ...)`
- `run_lcv(gwas1_mapped_file, gwas2_mapped_file, ...)`
- `run_flames(gwas_mapped_file, ...)`
- `run_t_sem(gwas_mapped_files, ...)`
- `run_p_sem(gwas_mapped_files, ...)`
- `run_gna(gwas1_mapped_file, gwas2_mapped_file, ...)`
- `run_parallel_lava(gwas_mapped_files, ...)`
- `run_parallel_mixer(gwas_mapped_files, ...)`
- `run_ccgwas(case_gwas_mapped_file, control_gwas_mapped_file, ...)`
- `run_multivariate_gwas(gwas_mapped_files, ...)`
- `run_clump_fine_mapping(gwas_mapped_file, ...)`
- `run_scdrs(gwas_mapped_file, ...)`

## 3) 示例

```r
library(PGAR)

# 读取并校验 rsid mapping 输出
x <- read_mapped_gwas("data/trait1.rsid_mapped.tsv")
validate_mapped_gwas(x)

# coloc: mapped GWAS + eQTL目录
res_coloc <- run_coloc(
  gwas_mapped_file = "data/trait1.rsid_mapped.tsv",
  eqtl_dir = "data/eqtl"
)

# HDL/LCV: 两个mapped GWAS
res_hdl <- run_hdl(
  gwas1_mapped_file = "data/trait1.rsid_mapped.tsv",
  gwas2_mapped_file = "data/trait2.rsid_mapped.tsv",
  hdl_script = "scripts/HDL.R",
  piecewise_path = "refs/UKB_imputed_SVD_eigen99_extraction"
)

res_lcv <- run_lcv(
  gwas1_mapped_file = "data/trait1.rsid_mapped.tsv",
  gwas2_mapped_file = "data/trait2.rsid_mapped.tsv",
  lcv_script = "scripts/LCV.R"
)

# LDSC: mapped GWAS -> munge -> ldsc
ldsc_log <- run_ldsc(
  gwas1_mapped_file = "data/trait1.rsid_mapped.tsv",
  gwas2_mapped_file = "data/trait2.rsid_mapped.tsv",
  munge_py = "~/tools/ldsc/munge_sumstats.py",
  ldsc_py = "~/tools/ldsc/ldsc.py",
  snplist = "~/tools/ldsc/w_hm3.snplist",
  ld_ref = "~/tools/ldsc/eur_w_ld_chr/",
  w_ld = "~/tools/ldsc/eur_w_ld_chr/"
)
```


## 4) 更多流程

```r
run_flames(
  gwas_mapped_file = "data/trait1.rsid_mapped.tsv",
  flames_script = "scripts/flames_pipeline.R"
)

run_parallel_lava(
  gwas_mapped_files = c("data/trait1.rsid_mapped.tsv", "data/trait2.rsid_mapped.tsv"),
  lava_script = "scripts/lava_parallel.R",
  cores = 8
)

run_clump_fine_mapping(
  gwas_mapped_file = "data/trait1.rsid_mapped.tsv",
  plink_bin = "plink",
  ref_bfile = "refs/1kg_eur"
)
```

> 注意：新增流程函数默认会优先调用你本地脚本中的 `run_*_pipeline`，若不存在则回退为 `Rscript` / 外部命令调用。
