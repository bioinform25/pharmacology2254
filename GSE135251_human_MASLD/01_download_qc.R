# GSE135251 - 인체 MASLD 스펙트럼 bulk RNA-seq (Govaere et al., n=216)
# Step 1: 다운로드 + 메타데이터 정리 + count 행렬 구성 + 저발현 필터링 + 샘플/메타 정합성 확인
# 작업 디렉토리는 pharmacology2254.Rproj를 연 상태(레포 루트)를 가정.

set.seed(42)

library(GEOquery)
library(tidyverse)

data_dir = "GSE135251_human_MASLD/data"
raw_dir = file.path(data_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# 이미 로컬에 받아둔 RAW.tar가 있으면 재사용 (216샘플 재다운로드 방지). 없으면 새로 받음.
existing_tar = "C:/Users/SAMSUNG/Desktop/5-1 1~2/GSE135251/GSE135251/GSE135251_RAW.tar"
local_tar = file.path(raw_dir, "GSE135251_RAW.tar")

if (!file.exists(local_tar)) {
  if (file.exists(existing_tar)) {
    file.copy(existing_tar, local_tar)
  } else {
    getGEOSuppFiles("GSE135251", baseDir = raw_dir, makeDirectory = FALSE)
  }
}

extract_dir = file.path(raw_dir, "GSE135251_RAW")
if (!dir.exists(extract_dir) || length(list.files(extract_dir)) == 0) {
  dir.create(extract_dir, showWarnings = FALSE)
  untar(local_tar, exdir = extract_dir)
}

# 메타데이터
gse = getGEO("GSE135251", destdir = data_dir, getGPL = TRUE)
metadata = pData(gse[[1]])
metadata = metadata[, c("disease:ch1", "fibrosis stage:ch1", "group in paper:ch1", "nas score:ch1", "Stage:ch1")]
colnames(metadata) = c("disease", "fibrosis_stage", "group", "NAS", "stage_category")

group_levels = c("Normal", "Steatosis_Only", "Early_NASH_Fibrosis", "Advanced_Fibrosis")

meta_final = metadata %>%
  mutate(analysis_group = case_when(
    disease == "Control" ~ "Normal",
    group == "NAFL" ~ "Steatosis_Only",
    group %in% c("NASH_F0-F1", "NASH_F2") ~ "Early_NASH_Fibrosis",
    group %in% c("NASH_F3", "NASH_F4") ~ "Advanced_Fibrosis",
    TRUE ~ "Others"
  )) %>%
  filter(analysis_group != "Others") %>%
  mutate(analysis_group = factor(analysis_group, levels = group_levels)) %>%
  arrange(analysis_group)

# 개별 샘플 count 파일 병합 (첫 컬럼: 유전자 ID, 두 번째 컬럼: count)
file_list = list.files(extract_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
stopifnot(length(file_list) > 0)

count_list = file_list %>% map(~ read.table(.x, header = TRUE, row.names = 1))
merged_counts = reduce(count_list, cbind)
sample_names = basename(file_list) %>% str_extract("GSM[0-9]+")
colnames(merged_counts) = sample_names

# 저발현 필터링: DESeq2 권장 방식 (최소 10개 샘플에서 count >= 5)
keep = rowSums(merged_counts >= 5) >= 10
merged_counts_filtered = merged_counts[keep, ]
message(sprintf("저발현 필터링: %d genes -> %d genes", nrow(merged_counts), nrow(merged_counts_filtered)))

if (any(duplicated(rownames(merged_counts_filtered)))) {
  merged_counts_filtered$gene_sum = rowSums(merged_counts_filtered)
  merged_counts_filtered = merged_counts_filtered %>%
    rownames_to_column("GeneID") %>%
    group_by(GeneID) %>%
    filter(gene_sum == max(gene_sum)) %>%
    ungroup() %>%
    select(-gene_sum) %>%
    column_to_rownames("GeneID")
}

common_samples = intersect(colnames(merged_counts_filtered), rownames(meta_final))
count_final = merged_counts_filtered[, common_samples]
meta_final = meta_final[common_samples, ]
count_final = count_final[, rownames(meta_final)]
stopifnot(all(colnames(count_final) == rownames(meta_final)))

saveRDS(list(counts = count_final, meta = meta_final, group_levels = group_levels),
        file.path(data_dir, "qc_ready.rds"))

message(sprintf("QC 완료: %d genes x %d samples", nrow(count_final), ncol(count_final)))
print(table(meta_final$analysis_group))
