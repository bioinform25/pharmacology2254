# GSE222576 - 마우스 CCl4 반복투여 간섬유화 시계열 (Affymetrix MoGene-2.0-st, n=3/군, 6개 시점)
# Step 1: CEL 파일 다운로드 + RMA 정규화 + 메타데이터 정리
# 작업 디렉토리는 pharmacology2254.Rproj를 연 상태(레포 루트)를 가정.

set.seed(42)

library(GEOquery)
library(oligo)
library(tidyverse)

data_dir = "GSE222576_mouse_CCl4/data"
raw_dir = file.path(data_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# CEL 원자료 다운로드 (RAW.tar 안에 샘플별 CEL.gz)
tar_file = file.path(raw_dir, "GSE222576_RAW.tar")
if (!file.exists(tar_file)) {
  getGEOSuppFiles("GSE222576", baseDir = raw_dir, makeDirectory = FALSE)
}

extract_dir = file.path(raw_dir, "CEL")
if (!dir.exists(extract_dir) || length(list.files(extract_dir)) == 0) {
  dir.create(extract_dir, showWarnings = FALSE)
  untar(tar_file, exdir = extract_dir)
}

# 메타데이터: title 컬럼에 그룹/시점 정보가 있음
# (예: "Control at 10 weeks, biological rep1", "CCl4-induced liver fibrosis for 2 weeks, biological rep1")
gse = getGEO("GSE222576", destdir = data_dir, GSEMatrix = TRUE, getGPL = FALSE)
pd = pData(gse[[1]])

group_levels = c("Control", "CCl4_2w", "CCl4_4w", "CCl4_6w", "CCl4_8w", "CCl4_10w")

meta = pd %>%
  rownames_to_column("gsm") %>%
  mutate(
    weeks = as.integer(str_extract(title, "(?<=for )\\d+(?= weeks)")),
    group = if_else(str_detect(title, "^Control"), "Control", paste0("CCl4_", weeks, "w")),
    group = factor(group, levels = group_levels)
  ) %>%
  select(gsm, title, group, weeks)
rownames(meta) = meta$gsm

# CEL 파일 경로를 GSM 순서에 맞춰 정렬
cel_files = list.files(extract_dir, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
cel_gsm = str_extract(basename(cel_files), "GSM[0-9]+")
names(cel_files) = cel_gsm
cel_files = cel_files[meta$gsm]
stopifnot(all(!is.na(cel_files)))

# RMA 정규화 (Gene ST array이므로 core/transcript-cluster 레벨로 요약)
raw_data = read.celfiles(cel_files)
eset = oligo::rma(raw_data, target = "core")
expr = exprs(eset)
colnames(expr) = meta$gsm

message(sprintf("RMA 정규화 완료: %d transcript clusters x %d samples", nrow(expr), ncol(expr)))

saveRDS(list(expr = expr, meta = meta, group_levels = group_levels),
        file.path(data_dir, "qc_ready.rds"))

print(table(meta$group))
