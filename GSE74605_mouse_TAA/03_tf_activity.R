# GSE74605 - Step 3: 상위 전사인자(TF) 활성 추론 (decoupleR + CollecTRI, mouse)
# GSE135251/GSE222576과 동일한 run_tf_activity() 함수 재사용

set.seed(42)

library(tidyverse)
library(limma)
library(pheatmap)
source("shared/R/pub_theme.R")
source("shared/R/tf_activity.R")

data_dir = "GSE74605_mouse_TAA/data"
fig_dir = "GSE74605_mouse_TAA/figures"

qc = readRDS(file.path(data_dir, "qc_ready.rds"))
deg = readRDS(file.path(data_dir, "deg_results.rds"))
expr = qc$expr
meta = qc$meta
group_levels = qc$group_levels
probe_map = deg$probe_map

common_probes = intersect(rownames(expr), probe_map$PROBEID)
expr_sub = expr[common_probes, ]
symbol_vec = probe_map$SYMBOL[match(common_probes, probe_map$PROBEID)]
expr_sym = avereps(expr_sub, ID = symbol_vec)

tf_result = run_tf_activity(expr_sym, organism = "mouse", cache_dir = "shared/cache")
message(sprintf("TF activity 계산 완료: %d TFs x %d samples", nrow(tf_result$wide), ncol(tf_result$wide)))

tf_group = summarize_tf_by_group(tf_result$wide, meta$group, group_levels)
tf_stats = test_tf_across_groups(tf_result$wide, meta$group, group_levels, ordered = TRUE)

n_sig = sum(tf_stats$padj < 0.05)
message(sprintf("시점간(Naive->6w) 유의한 TF (BH padj<0.05): %d / %d", n_sig, nrow(tf_stats)))

if (n_sig >= 5) {
  top_tf = tf_stats %>% filter(padj < 0.05) %>% slice_head(n = 30) %>% pull(TF)
  heatmap_title = sprintf("TF activity across TAA time course (BH padj<0.05, n=%d)", length(top_tf))
} else {
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
