# GSE74605 - 마우스 TAA(thioacetamide) 유발 간섬유화 시계열 (Illumina MouseRef-8 v2.0, GPL6885)
# 3번째 독립 데이터셋: GSE222576(CCl4)과는 다른 화학적 간독성물질로, "화학독성 축" 결과가
# CCl4 한 물질에만 국한된 우연이 아닌지 검증하기 위한 견고성(robustness) 검증용.
# Step 1: 비정규화(non-normalized) signal 다운로드 + quantile 정규화 + 메타데이터 정리

set.seed(42)

library(GEOquery)
library(limma)
library(tidyverse)

data_dir = "GSE74605_mouse_TAA/data"
raw_dir = file.path(data_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

if (length(list.files(raw_dir, pattern = "non-normalized")) == 0) {
  getGEOSuppFiles("GSE74605", baseDir = raw_dir, makeDirectory = FALSE)
}
nn_file = list.files(raw_dir, pattern = "non-normalized", full.names = TRUE)[1]

# 파일 구조: 상단 4줄은 주석, 5번째 줄이 헤더.
# 컬럼: ID_REF, (SampleName, "Detection Pval") 쌍이 샘플 수만큼 반복.
raw = read.delim(nn_file, skip = 4, header = TRUE, check.names = FALSE)
signal_cols = colnames(raw)[colnames(raw) != "ID_REF" & colnames(raw) != "Detection Pval"]

signal = raw[, signal_cols]
rownames(signal) = raw$ID_REF
signal = as.matrix(signal)
storage.mode(signal) = "numeric"

# 샘플명(예: "Naive A", "Week 1 A", "Week 4 D") -> 진행단계 그룹
group_levels = c("Naive", "TAA_1w", "TAA_4w", "TAA_6w")
sample_meta = data.frame(sample = signal_cols) %>%
  mutate(
    group = case_when(
      str_detect(sample, "^Naive") ~ "Naive",
      str_detect(sample, "^Week 1 ") ~ "TAA_1w",
      str_detect(sample, "^Week 4 ") ~ "TAA_4w",
      str_detect(sample, "^Week 6 ") ~ "TAA_6w"
    ),
    group = factor(group, levels = group_levels)
  )
rownames(sample_meta) = sample_meta$sample
stopifnot(!any(is.na(sample_meta$group)))

# log2 변환 + quantile 정규화 (Illumina non-normalized signal에 대한 표준 처리 방식;
# neqc()에 필요한 negative control probe 정보가 이 처리된 파일에는 없어 quantile 정규화로 대체)
signal[signal <= 0] = NA
log_signal = log2(signal)
norm_expr = normalizeBetweenArrays(log_signal, method = "quantile")
norm_expr = norm_expr[complete.cases(norm_expr), ]

message(sprintf("정규화 완료: %d probes x %d samples", nrow(norm_expr), ncol(norm_expr)))

saveRDS(list(expr = norm_expr, meta = sample_meta, group_levels = group_levels),
        file.path(data_dir, "qc_ready.rds"))

print(table(sample_meta$group))
