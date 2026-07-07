# 상위 전사인자(TF) 활성 추론 공용 래퍼.
# decoupleR + CollecTRI (Badia-i-Mompel et al. 2023) 기반, Bioconductor 공식 bulk RNA-seq
# 워크플로우(https://saezlab.github.io/decoupleR/articles/tf_bk.html)를 그대로 따름.
# GSE135251(RNA-seq VST)과 GSE222576(microarray RMA) 둘 다, "정규화된 발현행렬"만 넣으면
# 동일한 방식으로 TF activity를 계산하도록 설계 -> 두 데이터셋 스크립트가 이 함수 하나를 공유.
#
# 주의: decoupleR::get_collectri()는 내부적으로 OmnipathR::get_db("organisms")를 거치는데,
# 현재 설치된 OmnipathR(3.18.4)에서 이 organism-resolution 로직에 버그가 있어 항상 에러가 남
# (full_join 시 ncbi_tax_id 컬럼 누락). 대신 OmniPath REST API를 직접 호출해서
# decoupleR::get_collectri()와 동일한 구성 규칙(mor = is_stimulation을 +1/-1로,
# COMPLEX 소스는 AP1/NFKB로 축약)을 그대로 재현함 -> 결과는 동일하되 버그를 우회.

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(decoupleR)

ORGANISM_TAXID = c(human = 9606, mouse = 10090)

# decoupleR::get_collectri()와 동일한 로직으로 CollecTRI 네트워크를 직접 REST API에서 구성.
build_collectri_direct = function(organism = "human") {
  taxid = if (is.numeric(organism)) organism else ORGANISM_TAXID[[organism]]
  url = sprintf(
    "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&organisms=%d&fields=is_stimulation,is_inhibition",
    taxid
  )
  raw = read_tsv(url, show_col_types = FALSE)

  cols = c("source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
  interactions = raw[!str_detect(raw$source, "COMPLEX"), cols]
  complexes = raw[str_detect(raw$source, "COMPLEX"), cols] %>%
    mutate(source_genesymbol = case_when(
      str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
      str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB"
    ))

  # 일부 mouse 유전자는 CollecTRI/Omnipath에 정식 gene symbol이 없어 UniProt accession이
  # source_genesymbol/target_genesymbol에 그대로 남는 경우가 있음 (예: "A0A0R4J082").
  # 이런 항목이 "TF"로 그림/표에 노출되면 논문 신뢰도를 해치므로 accession 형태는 제외.
  uniprot_pattern = "^[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-[0-9]+)?$"

  rbind(interactions, complexes) %>%
    distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) %>%
    mutate(mor = case_when(is_stimulation == 1 ~ 1, is_stimulation == 0 ~ -1)) %>%
    filter(!is.na(mor), !is.na(source_genesymbol)) %>%
    filter(!str_detect(source_genesymbol, uniprot_pattern),
           !str_detect(target_genesymbol, uniprot_pattern)) %>%
    select(source = source_genesymbol, target = target_genesymbol, mor)
}

# CollecTRI 네트워크를 매번 새로 받지 않도록 로컬 캐시.
get_collectri_cached = function(organism = "human",
                                 cache_dir = "shared/cache") {
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_file = file.path(cache_dir, paste0("collectri_", organism, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  net = build_collectri_direct(organism)
  saveRDS(net, cache_file)
  net
}

# expr_matrix: gene symbol(row) x sample(column), 정규화된 값(VST/logCPM/RMA 등).
# organism: "human" 또는 "mouse" (get_collectri가 오솔로그 매핑을 내부적으로 처리).
# minsize: 최소 몇 개의 target 유전자가 데이터에 존재해야 그 TF를 평가할지 (decoupleR 기본 권장 5).
run_tf_activity = function(expr_matrix, organism = "human", minsize = 5,
                            cache_dir = "shared/cache") {
  net = get_collectri_cached(organism, cache_dir)

  acts = decoupleR::run_ulm(
    mat = as.matrix(expr_matrix),
    net = net,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = minsize
  )

  acts_wide = acts %>%
    select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()

  list(long = acts, wide = acts_wide, network = net)
}

# 그룹(진행단계/timepoint)별 평균 TF activity 행렬 -> heatmap 입력용.
# group_levels로 순서를 명시적으로 고정해서 heatmap 컬럼 순서가 항상 생물학적 순서를 따르게 함.
summarize_tf_by_group = function(acts_wide, group_vector, group_levels) {
  stopifnot(length(group_vector) == ncol(acts_wide))
  group_vector = factor(group_vector, levels = group_levels)
  sapply(group_levels, function(g) rowMeans(acts_wide[, group_vector == g, drop = FALSE], na.rm = TRUE))
}

# 그룹 간 TF activity 차이에 대한 실제 유의성 검정 (진행단계 순서를 가진 그룹이면 ordered ANOVA,
# 아니면 일반 ANOVA). "top N을 그냥 눈으로 골랐다"가 아니라 BH-보정 p-value로 랭킹을 정당화하기 위함.
test_tf_across_groups = function(acts_wide, group_vector, group_levels, ordered = TRUE) {
  group_vector = factor(group_vector, levels = group_levels, ordered = ordered)
  pvals = apply(acts_wide, 1, function(tf_scores) {
    fit = tryCatch(aov(tf_scores ~ group_vector), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    summary(fit)[[1]][["Pr(>F)"]][1]
  })
  data.frame(TF = rownames(acts_wide), p_value = pvals) %>%
    filter(!is.na(p_value)) %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj)
}

# 두 독립 데이터셋(플랫폼이 달라 raw 레벨 병합이 부적절한 경우)에서 나온 TF별 p-value를
# 표준 메타분석 방식(Fisher's / Stouffer's method)으로 결합. TF 매칭은 공통 gene symbol 기준.
fisher_combine = function(pvalues) {
  pvalues = pvalues[!is.na(pvalues) & pvalues > 0]
  stat = -2 * sum(log(pvalues))
  pchisq(stat, df = 2 * length(pvalues), lower.tail = FALSE)
}

stouffer_combine = function(pvalues, weights = NULL) {
  pvalues = pmin(pmax(pvalues, 1e-300), 1 - 1e-16)
  z = qnorm(1 - pvalues)
  if (is.null(weights)) weights = rep(1, length(pvalues))
  z_comb = sum(weights * z) / sqrt(sum(weights^2))
  1 - pnorm(z_comb)
}
