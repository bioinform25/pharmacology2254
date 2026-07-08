# 사람-마우스 유전자 심볼을 정식 ortholog 테이블로 매칭.
# 이전 버전은 toupper()로 대소문자만 맞춰서 매칭했는데, 이는 "사람-마우스 심볼이 대소문자만
# 다르다"는 가정이 항상 참이 아니라서 틀린 경우가 있음 (예: 사람 TP53의 마우스 ortholog는
# "Tp53"이 아니라 "Trp53" — p53 계열 전체가 마우스에서 "Trp" 접두사를 씀).
# NCBI Gene ortholog 데이터를 담은 Bioconductor 패키지 Orthology.eg.db를 사용 -
# biomaRt처럼 매 실행마다 Ensembl 서버에 실시간 쿼리하지 않아도 되어 재현성이 더 좋음.
#
# 매칭 체인: 마우스 SYMBOL -> 마우스 ENTREZID(org.Mm.eg.db) -> 사람 ENTREZID(Orthology.eg.db)
#           -> 사람 SYMBOL(org.Hs.eg.db). 사람 쪽은 이미 HGNC 표준 심볼이라 변환 불필요.

library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

# 마우스 gene symbol 벡터 -> 사람 ortholog gene symbol 벡터 (매칭 안 되면 NA).
mouse_symbol_to_human_symbol = function(mouse_symbols) {
  mm_entrez = AnnotationDbi::mapIds(org.Mm.eg.db, keys = mouse_symbols, column = "ENTREZID",
                                     keytype = "SYMBOL", multiVals = "first")

  valid = !is.na(mm_entrez)
  ortho = AnnotationDbi::select(Orthology.eg.db, keys = as.character(mm_entrez[valid]),
                                 keytype = "Mus_musculus", columns = "Homo_sapiens")
  lut = setNames(as.character(ortho$Homo_sapiens), as.character(ortho$Mus_musculus))
  hs_entrez = unname(lut[as.character(mm_entrez)])

  hs_sym = rep(NA_character_, length(mouse_symbols))
  has_entrez = !is.na(hs_entrez)
  hs_sym[has_entrez] = AnnotationDbi::mapIds(org.Hs.eg.db, keys = hs_entrez[has_entrez],
                                              column = "SYMBOL", keytype = "ENTREZID",
                                              multiVals = "first")
  names(hs_sym) = mouse_symbols
  hs_sym
}
