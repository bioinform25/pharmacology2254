# GSE135251 - Step 2: DESeq2 기반 DEG 분석
# 유의기준은 분석 전에 고정: padj(BH) < 0.05 & |log2FC| > 1. 결과가 안 나와도 사후에 완화하지 않음
# (과거 GSE202379 분석에서 padj 기준 0개가 나오자 p<0.01로 낮춘 사례가 있었음 - 이번엔 반복하지 않음).

set.seed(42)

library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
source("shared/R/pub_theme.R")

SIG_PADJ = 0.05
SIG_LFC = 1

data_dir = "GSE135251_human_MASLD/data"
fig_dir = "GSE135251_human_MASLD/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

qc = readRDS(file.path(data_dir, "qc_ready.rds"))
count_final = qc$counts
meta_final = qc$meta
group_levels = qc$group_levels

dds = DESeqDataSetFromMatrix(countData = count_final, colData = meta_final, design = ~ analysis_group)
dds = dds[rowSums(counts(dds)) > 0, ]

vsd = vst(dds, blind = FALSE)
saveRDS(vsd, file.path(data_dir, "vsd.rds"))

# QC: PCA plot으로 그룹 분리가 실제로 발현 패턴에 반영되는지 확인
pcaData = plotPCA(vsd, intgroup = "analysis_group", returnData = TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

p_pca = ggplot(pcaData, aes(PC1, PC2, color = analysis_group)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = pub_palette_progression, name = "Progression") +
  theme_pub() +
  ggtitle("PCA: MASLD Progression Stages (GSE135251)")
ggsave_pub(file.path(fig_dir, "01_pca.pdf"), p_pca)

# DEG: Steatosis_Only를 기준(reference)으로 섬유화 진행 단계와 비교
dds$analysis_group = relevel(dds$analysis_group, ref = "Steatosis_Only")
dds = DESeq(dds)

contrasts_list = list(
  Early_vs_Steatosis    = "analysis_group_Early_NASH_Fibrosis_vs_Steatosis_Only",
  Advanced_vs_Steatosis = "analysis_group_Advanced_Fibrosis_vs_Steatosis_Only"
)

deg_results = list()
for (nm in names(contrasts_list)) {
  res_shrunk = lfcShrink(dds, coef = contrasts_list[[nm]], type = "apeglm")
  res_df = as.data.frame(res_shrunk) %>%
    rownames_to_column("GeneID") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  n_sig = sum(res_df$padj < SIG_PADJ & abs(res_df$log2FoldChange) > SIG_LFC)
  message(sprintf("[%s] 유의 유전자 (padj<%.2f & |log2FC|>%.1f): %d개 / 검정된 %d개 중",
                   nm, SIG_PADJ, SIG_LFC, n_sig, nrow(res_df)))

  deg_results[[nm]] = res_df

  p_volcano = EnhancedVolcano(res_df,
    lab = res_df$GeneID, x = "log2FoldChange", y = "padj",
    title = nm, subtitle = "GSE135251 (human MASLD)",
    pCutoff = SIG_PADJ, FCcutoff = SIG_LFC,
    pointSize = 2.5, labSize = 0,  # 유전자가 많아 라벨은 생략, TF 분석 단계에서 별도 라벨링
    col = c("grey70", "#009E73", "#0072B2", "#D55E00"),
    legendLabels = c("NS", "Log2FC", "padj", "padj & Log2FC"))
  ggsave_pub(file.path(fig_dir, sprintf("02_volcano_%s.pdf", nm)), p_volcano, width = 7, height = 7)
}

saveRDS(deg_results, file.path(data_dir, "deg_results.rds"))
message("DEG 분석 완료.")
