# GSE74605 - Step 2: limma 기반 DEG 분석 (각 TAA 시점 vs Naive)
# 유의기준은 다른 두 데이터셋과 동일하게 사전 고정: padj(BH) < 0.05 & |log2FC| > 1

set.seed(42)

library(limma)
library(tidyverse)
library(EnhancedVolcano)
library(illuminaMousev2.db)
library(AnnotationDbi)
source("shared/R/pub_theme.R")

SIG_PADJ = 0.05
SIG_LFC = 1

data_dir = "GSE74605_mouse_TAA/data"
fig_dir = "GSE74605_mouse_TAA/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

qc = readRDS(file.path(data_dir, "qc_ready.rds"))
expr = qc$expr
meta = qc$meta
group_levels = qc$group_levels

# Illumina probe ID(ILMN_xxxx) -> gene symbol
probe_map = AnnotationDbi::select(illuminaMousev2.db,
                                   keys = rownames(expr),
                                   columns = "SYMBOL", keytype = "PROBEID")
probe_map = probe_map %>% filter(!is.na(SYMBOL), !duplicated(PROBEID))

# QC: PCA
pca = prcomp(t(expr), scale. = FALSE)
pca_df = data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = meta$group)
percentVar = round(100 * summary(pca)$importance[2, 1:2])

p_pca = ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_pub() +
  ggtitle("PCA: TAA Time Course (GSE74605)")
ggsave_pub(file.path(fig_dir, "01_pca.pdf"), p_pca)

# limma: Naive를 기준(reference)으로 각 시점 비교
design = model.matrix(~ group, data = meta)
fit = lmFit(expr, design)
fit = eBayes(fit)

timepoint_coefs = colnames(design)[-1]
deg_results = list()

for (coef in timepoint_coefs) {
  nm = str_remove(coef, "^group")
  tt = topTable(fit, coef = coef, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("PROBEID") %>%
    left_join(probe_map, by = "PROBEID") %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(adj.P.Val)

  n_sig = sum(tt$adj.P.Val < SIG_PADJ & abs(tt$logFC) > SIG_LFC)
  message(sprintf("[%s vs Naive] 유의 유전자 (padj<%.2f & |log2FC|>%.1f): %d개 / 검정된 %d개 중",
                   nm, SIG_PADJ, SIG_LFC, n_sig, nrow(tt)))

  deg_results[[nm]] = tt

  p_volcano = EnhancedVolcano(tt,
    lab = tt$SYMBOL, x = "logFC", y = "adj.P.Val",
    title = paste(nm, "vs Naive"), subtitle = "GSE74605 (mouse TAA)",
    pCutoff = SIG_PADJ, FCcutoff = SIG_LFC,
    pointSize = 2.5, labSize = 3.0)
  ggsave_pub(file.path(fig_dir, sprintf("02_volcano_%s.pdf", nm)), p_volcano, width = 7, height = 7)
}

saveRDS(list(deg_results = deg_results, probe_map = probe_map),
        file.path(data_dir, "deg_results.rds"))
message("DEG 분석 완료.")
