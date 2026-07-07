# GSE135251 - Step 3: 상위 전사인자(TF) 활성 추론 (decoupleR + CollecTRI)
# 진행단계별 TF activity heatmap + 그룹간 유의성 검정(ANOVA, BH-보정)

set.seed(42)

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)
source("shared/R/pub_theme.R")
source("shared/R/tf_activity.R")

data_dir = "GSE135251_human_MASLD/data"
fig_dir = "GSE135251_human_MASLD/figures"

qc = readRDS(file.path(data_dir, "qc_ready.rds"))
vsd = readRDS(file.path(data_dir, "vsd.rds"))
meta_final = qc$meta
group_levels = qc$group_levels

expr = assay(vsd)

# CollecTRI 네트워크는 gene symbol 기준 -> Ensembl ID를 symbol로 변환
gene_map = bitr(rownames(expr), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
gene_map = gene_map[!duplicated(gene_map$ENSEMBL), ]
expr_sym = expr[gene_map$ENSEMBL, ]
rownames(expr_sym) = gene_map$SYMBOL
expr_sym = expr_sym[!duplicated(rownames(expr_sym)), ]

tf_result = run_tf_activity(expr_sym, organism = "human", cache_dir = "shared/cache")
message(sprintf("TF activity 계산 완료: %d TFs x %d samples", nrow(tf_result$wide), ncol(tf_result$wide)))

tf_group = summarize_tf_by_group(tf_result$wide, meta_final$analysis_group, group_levels)
tf_stats = test_tf_across_groups(tf_result$wide, meta_final$analysis_group, group_levels, ordered = TRUE)

n_sig = sum(tf_stats$padj < 0.05)
message(sprintf("그룹간(Normal->Advanced) 유의한 TF (BH padj<0.05): %d / %d", n_sig, nrow(tf_stats)))

if (n_sig >= 5) {
  top_tf = tf_stats %>% filter(padj < 0.05) %>% slice_head(n = 30) %>% pull(TF)
  heatmap_title = sprintf("TF activity across MASLD progression (BH padj<0.05, n=%d)", length(top_tf))
} else {
  # 통계적으로 유의한 TF가 부족하면, 유의하다고 없는 것을 있는 척하지 않고
  # nominal p-value 상위만 "참고용"으로 표시함을 제목에 명시.
  message("BH padj<0.05를 통과한 TF가 5개 미만 -> nominal p-value 상위 20개를 참고용으로만 표시")
  top_tf = tf_stats %>% slice_head(n = 20) %>% pull(TF)
  heatmap_title = "TF activity (top 20 by nominal p; not BH-significant)"
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(fig_dir, "03_tf_activity_heatmap.pdf"), width = 7, height = 9)
pheatmap(tf_group[top_tf, , drop = FALSE],
         color = pub_heatmap_colors(),
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         border_color = NA,
         fontsize_row = 8,
         main = heatmap_title)
dev.off()

saveRDS(list(tf_result = tf_result, tf_group = tf_group, tf_stats = tf_stats),
        file.path(data_dir, "tf_activity.rds"))
write_csv(tf_stats, file.path(data_dir, "tf_stats.csv"))

message("TF activity 분석 완료.")
print(head(tf_stats, 10))
