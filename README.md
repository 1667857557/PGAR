# PGAR

一个用于 **GWAS-post** 分析的 R 包脚手架，统一封装多种分析流程：

- `run_coloc`: 从原始 GWAS summary + 本地 eQTL 文件夹批量跑共定位。
- `run_hdl`: 跑 HDL 遗传相关分析。
- `run_ldsc`: 调用 `ldsc.py` 跑 LDSC rg。
- `run_Mapgen`: 调用 MAPGEN 可执行程序。
- `run_lcv`: 跑 LCV 因果方向分析。

## 安装

```r
devtools::install_local(".")
```

## 示例

```r
library(PGAR)

res_coloc <- run_coloc(
  gwas_file = "data/gwas_sumstats.tsv",
  eqtl_dir = "data/eqtl",
  output_dir = "results/coloc",
  gwas_colmap = list(SNP = "rsid", beta = "BETA", se = "SE", p = "P"),
  eqtl_colmap = list(SNP = "snp", beta = "beta", se = "se", p = "pval")
)

res_hdl <- run_hdl(
  gwas1_file = "data/trait1.tsv",
  gwas2_file = "data/trait2.tsv",
  hdl_script = "scripts/HDL.R",
  piecewise_path = "refs/UKB_imputed_SVD_eigen99_extraction"
)

ldsc_log <- run_ldsc(
  munged_sumstats_1 = "data/trait1.sumstats.gz",
  munged_sumstats_2 = "data/trait2.sumstats.gz",
  ldsc_py = "~/tools/ldsc/ldsc.py",
  ld_ref = "~/tools/ldsc/eur_w_ld_chr/",
  w_ld = "~/tools/ldsc/eur_w_ld_chr/"
)
```

> 注意：`run_hdl` / `run_lcv` 默认依赖你本地已加载的脚本函数（例如 `HDL.rg`, `RunLCV`）；`run_ldsc` / `run_Mapgen` 依赖外部命令行工具。
