# Cross-dataset TF 메타분석 (3개 독립 데이터셋)
# 1) GSE135251 - 인체 MASLD, 대사성 만성질환 축 (n=216, RNA-seq)
# 2) GSE222576 - 마우스 CCl4, 화학독성 축 #1 (n=18, microarray)
# 3) GSE74605  - 마우스 TAA,  화학독성 축 #2 (n=13, microarray) - CCl4와는 다른 독성기전의
#    화학적 간독성물질로, "화학독성 축" 결과가 CCl4 한 물질에 국한된 우연이 아닌지 검증
#
# 세 데이터셋은 플랫폼/종/독성기전이 모두 달라 raw expression 레벨에서 합치지 않음.
# 각자 독립적으로 계산한 TF activity의 그룹간 유의성(p-value)을 Fisher's method로 결합하되,
# 1차 선별 기준은 "세 데이터셋 모두에서 독립적으로 유의(padj<0.05) + 방향까지 일치"라는
# 더 엄격한 기준을 사용함 (Fisher's method는 한쪽의 극단적 p-value에 끌려갈 수 있어서;
# 2-데이터셋 분석 때 YAP1: human padj~1e-18, mouse padj=0.74인데 결합만 유의했던 사례 참고).

set.seed(42)

library(tidyverse)
library(pheatmap)
source("shared/R/pub_theme.R")
source("shared/R/tf_activity.R")
source("shared/R/ortholog_map.R")

fig_dir = "cross_dataset_metaanalysis/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

datasets = list(
  human = list(dir = "GSE135251_human_MASLD", label = "Human MASLD (GSE135251)", organism = "human"),
  ccl4  = list(dir = "GSE222576_mouse_CCl4",   label = "Mouse CCl4 (GSE222576)",  organism = "mouse"),
  taa   = list(dir = "GSE74605_mouse_TAA",     label = "Mouse TAA (GSE74605)",    organism = "mouse")
)

# 마우스 심볼을 사람 ortholog 심볼로 변환 (정식 Orthology.eg.db 매칭 - toupper() 대소문자 변환이 아님).
# 사람 데이터는 이미 HGNC 표준 심볼이라 그대로 사용. 하나의 사람 ortholog로 여러 마우스 심볼이
# 매핑되는 모호한 경우(paralog 등)는 잘못된 병합을 막기 위해 제외.
canonicalize_symbols = function(symbols, organism) {
  if (organism == "human") return(symbols)
  mapped = mouse_symbol_to_human_symbol(symbols)
  mapped[duplicated(mapped) | duplicated(mapped, fromLast = TRUE)] = NA
  unname(mapped)
}

# 각 데이터셋의 TF별 p-value/padj + "진행에 따른 방향"(baseline -> 가장 진행된 그룹, activity 평균 차이의 부호)
load_dataset_tf = function(spec) {
  stats = read_csv(file.path(spec$dir, "data", "tf_stats.csv"), show_col_types = FALSE) %>%
    mutate(TF = canonicalize_symbols(TF, spec$organism)) %>%
    filter(!is.na(TF)) %>%
    dplyr::select(TF, p_value, padj)
  act = readRDS(file.path(spec$dir, "data", "tf_activity.rds"))
  tf_group = act$tf_group
  rownames(tf_group) = canonicalize_symbols(rownames(tf_group), spec$organism)
  tf_group = tf_group[!is.na(rownames(tf_group)), , drop = FALSE]
  direction = sign(tf_group[, ncol(tf_group)] - tf_group[, 1])
  stats %>% filter(TF %in% names(direction)) %>% mutate(direction = direction[TF])
}

per_dataset = map(datasets, load_dataset_tf)

merged = reduce2(
  per_dataset, names(per_dataset),
  function(acc, df, nm) {
    df = rename_with(df, ~ paste0(., "_", nm), -TF)
    if (is.null(acc)) df else inner_join(acc, df, by = "TF")
  },
  .init = NULL
)

message(sprintf("세 데이터셋 공통 평가 TF: %d개", nrow(merged)))

p_cols = paste0("p_value_", names(datasets))
padj_cols = paste0("padj_", names(datasets))
dir_cols = paste0("direction_", names(datasets))

combined = merged %>%
  rowwise() %>%
  mutate(p_fisher = fisher_combine(c_across(all_of(p_cols)))) %>%
  ungroup() %>%
  mutate(
    padj_fisher = p.adjust(p_fisher, method = "BH"),
    all_significant = if_all(all_of(padj_cols), ~ . < 0.05),
    direction_concordant = (rowSums(across(all_of(dir_cols)) == 1) == length(dir_cols)) |
                            (rowSums(across(all_of(dir_cols)) == -1) == length(dir_cols)),
    true_consensus = all_significant & direction_concordant
  ) %>%
  arrange(desc(true_consensus), padj_fisher)

write_csv(combined, "cross_dataset_metaanalysis/combined_tf_results.csv")

message(sprintf("Fisher 결합 유의(BH padj<0.05): %d / %d", sum(combined$padj_fisher < 0.05), nrow(combined)))
message(sprintf("세 데이터셋 모두 유의: %d / %d", sum(combined$all_significant), nrow(combined)))
message(sprintf("세 데이터셋 모두 유의 + 방향 일치 (최종 consensus): %d / %d <- 논문 핵심 후보군",
                 sum(combined$true_consensus), nrow(combined)))

top_consensus = combined %>% filter(true_consensus) %>% slice_head(n = 15) %>% pull(TF)
print(combined %>% dplyr::select(TF, all_of(padj_cols), padj_fisher, direction_concordant, true_consensus) %>% head(15))

# Figure A: 상위 consensus TF들의 -log10(padj)를 데이터셋별로 비교하는 grouped bar plot
# (2D scatter는 3개 축을 동시에 못 보여주므로, 3-way 비교에는 grouped bar가 더 명확함)
bar_df = combined %>%
  filter(TF %in% top_consensus) %>%
  dplyr::select(TF, all_of(padj_cols)) %>%
  pivot_longer(-TF, names_to = "dataset", values_to = "padj") %>%
  mutate(
    dataset = str_remove(dataset, "^padj_"),
    dataset_label = map_chr(dataset, ~ datasets[[.]]$label),
    TF = factor(TF, levels = top_consensus)
  )

p_bar = ggplot(bar_df, aes(x = TF, y = -log10(padj), fill = dataset_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = NULL) +
  labs(x = NULL, y = "-log10(padj)",
       title = "Top consensus TFs: significant + direction-concordant across all 3 datasets") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_pub(file.path(fig_dir, "01_consensus_significance_bar.pdf"), p_bar, width = 10, height = 6)

# Figure B: 상위 consensus TF activity를 세 데이터셋 나란히 (스터디별 개별 z-score, 절대 스케일 차이 방지)
tf_groups = map(datasets, function(spec) {
  act = readRDS(file.path(spec$dir, "data", "tf_activity.rds"))
  m = act$tf_group
  rownames(m) = canonicalize_symbols(rownames(m), spec$organism)
  m[!is.na(rownames(m)), , drop = FALSE]
})

zscore_rows = function(mat) {
  m = rowMeans(mat, na.rm = TRUE)
  s = apply(mat, 1, sd, na.rm = TRUE)
  s[s == 0] = 1
  sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

avail_tf = purrr::reduce(map(tf_groups, rownames), intersect, .init = top_consensus)

z_blocks = map2(tf_groups, names(datasets), function(m, nm) {
  z = zscore_rows(m[avail_tf, , drop = FALSE])
  colnames(z) = paste0(nm, "_", colnames(m))
  z
})
combined_mat = purrr::reduce(z_blocks, cbind)
gaps = cumsum(map_int(z_blocks, ncol))[-length(z_blocks)]

pdf(file.path(fig_dir, "02_consensus_tf_heatmap.pdf"), width = 11, height = 6)
pheatmap(combined_mat,
         color = pub_heatmap_colors(),
         scale = "none",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         gaps_col = gaps,
         border_color = NA,
         fontsize_row = 9,
         main = "Consensus TF activity across 3 independent liver fibrosis datasets (top 15)")
dev.off()

message("Cross-dataset 메타분석(3개 데이터셋) 완료.")
