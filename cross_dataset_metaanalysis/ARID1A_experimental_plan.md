# ARID1A / THAP11 실험 검증 계획

> 두 후보 모두 재검증을 거쳤고 각자 다른 강점/약점이 있어(ARID1A: 기전적 근거 있으나 방향 해석이 복잡함 / THAP11: 문헌상 가장 깨끗하지만 druggability 없음) 어느 한쪽에 베팅하기보다 **Phase 1(세포주 수준)에서 둘 다 병행 검증**하고, 그 결과로 Phase 2 이후 어느 쪽에 자원을 집중할지 결정하는 구조로 계획을 짰다.

## 배경 가설 (2차 문헌 재확인 후 정정됨)
ARID1A는 SMARCA4(BRG1)와 함께 SWI/SNF(BAF) 크로마틴 리모델링 복합체를 구성한다. SMARCA4는 이미 간 HSC-근섬유아세포 전환의 직접적 드라이버로 검증되어 있다(Li et al., Cell Death Dis 2023 — Lrat-Cre;Smarca4^fl/fl 마우스에서 HSC-특이적 결손이 여러 CCl4/BDL 모델에서 섬유화를 완화).

**중요 정정**: 처음에는 "ARID1A는 간에서 문헌이 전혀 없다"고 판단했으나, 검색어를 바꿔(BAF250a, TGF-β/Smad, cholangiocarcinoma) 재확인한 결과 이는 틀렸다. 실제로는:
- **Hepatology 2019** (PMID 30584660): 간세포-특이적 Arid1a KO 마우스가 고지방식이 유도 지방간 **및 섬유화**가 더 심해짐 → ARID1A는 간세포에서 **보호적** 역할.
- **Cell Reports 2022** (PMC9808599): Arid1a 돌연변이가 담관세포에서 보호적인 TGF-β-Smad4 경로를 억제 → 담관암 유발 → 여기서도 ARID1A는 **보호적** 역할.

두 논문 모두 ARID1A **소실**이 나쁜 결과로 이어진다(보호적 역할)는 공통점이 있다. 그런데 저희 bulk-tissue 데이터를 직접 확인한 결과, ARID1A activity는 3개 데이터셋 전부에서 섬유화 진행에 따라 **증가**하는 방향이었다(direction=1, 3/3 일치) — 즉 표면적으로는 기존 문헌과 반대 방향처럼 보인다.

**이 불일치를 설명할 수 있는 두 가지 가능성**:
1. **세포구성 혼입(composition confound)**: bulk 데이터는 조직 전체의 평균이라, 섬유화가 진행될수록 조직 내 활성화된 HSC 비율 자체가 늘어난다. ARID1A가 HSC에서 원래도 상대적으로 높게 발현/활성화되는 유전자라면, 진짜 조절 변화 없이 "HSC 비율 증가"만으로도 bulk 신호가 올라간 것처럼 보일 수 있다.
2. **세포종류-특이적 반대 기능**: ARID1A가 간세포/담관세포(상피세포)에서는 보호적이지만, HSC(중간엽계 세포)에서는 반대로 섬유화를 촉진하는 세포종류-특이적 이중 기능을 가질 가능성 — 크로마틴 리모델러가 세포 맥락에 따라 정반대 기능을 하는 사례는 드물지 않음.

**두 기존 논문 모두 HSC 자체에서 ARID1A를 직접 테스트하지 않았다** — 따라서 "ARID1A가 HSC-자율적으로(cell-autonomous) 섬유화를 촉진하는가"라는 질문은 여전히 열려 있다. 이게 바로 Phase 1 실험이 필요한 이유이자, 이 실험이 답해야 할 첫 번째 질문이 됐다: **방향 불일치가 (1)구성 혼입 때문인지 (2)진짜 세포종류-반대 기능 때문인지 구분하는 것**.

**핵심 가설 (재구성)**: ARID1A는 간세포/담관세포에서는 보호적이지만, HSC에서는 SMARCA4와 함께(같은 BAF 복합체 안에서) 활성화를 촉진하는 세포종류-특이적 반대 기능을 가진다. 이 hypotheses가 맞다면, LX-2(HSC)에서의 knockdown 결과는 기존 간세포 논문과 반대 방향(즉 ARID1A KD가 섬유화 마커를 "감소"시켜야 함)으로 나와야 한다.

## Phase 1 — In vitro 발현/기능 확인 (우선순위 1)

두 세포주(LX-2 = HSC, AML12 = 간세포) × 두 타겟(ARID1A, THAP11)을 병행. THAP11은 문헌이 전혀 없으므로 "보호적/촉진적"에 대한 기존 가설이 없다 — Phase 1은 THAP11에 대해서는 순수 탐색(방향조차 모름), ARID1A에 대해서는 가설 검증(간세포=보호적 vs HSC=촉진적 반대 기능)이라는 성격 차이가 있다.

1. **발현 프로파일링**: LX-2를 quiescent → TGF-β1 자극 활성화 시간경과로 배양, ARID1A/THAP11 mRNA·단백질(qPCR, Western)이 활성화와 함께 증가하는지 확인 — bioinformatics 예측(둘 다 진행에 따라 증가) 검증.
2. **세포종류 대조 실험 (우선순위 최상)**: 같은 knockdown(ARID1A, THAP11 각각)을 LX-2(HSC)와 AML12(간세포) 양쪽에서 동시에 수행. 간세포 쪽은 ARID1A는 Hepatology 2019 지표(지방생성/지방산산화 유전자), THAP11은 사전 지표가 없으므로 일반 손상/스트레스 마커(예: CHOP, BAX)로 시작. HSC 쪽은 섬유화 마커(ACTA2/α-SMA, COL1A1, TIMP1, SERPINE1). 방향이 세포종류별로 갈리는지가 이 실험 전체의 성패를 가르므로 다른 실험보다 먼저 진행.
3. **Loss-of-function 표현형**: siRNA로 ARID1A, THAP11 각각 knockdown(대조군: scrambled siRNA) → TGF-β1 자극 유무 각각에서 섬유화 마커 mRNA+단백질, 증식(CCK-8/BrdU), 이동능(transwell) 측정.
   - *SMARCA4 2023 논문과 동일한 readout 패널을 그대로 사용 — ARID1A 결과를 SMARCA4 결과와 나란히 비교할 수 있어 "같은 복합체, 같은 표현형"이라는 주장이 훨씬 설득력 있어짐.*
4. **Gain-of-function**: 과발현만으로 quiescent LX-2에서 활성화 마커가 유도되는지 확인 (충분조건 테스트, 두 타겟 각각).

### Phase 1 일정 (총 8주 — 2개 타겟 x 2개 세포주로 확장되어 기존안(4-6주)보다 다소 늘어남)

| 주차 | 작업 |
|---|---|
| 1주 | siRNA(ARID1A/THAP11/scrambled), 과발현 플라스미드, 항체, 프라이머 주문 · LX-2/AML12 세포주 확보 및 계대배양 시작 |
| 2주 | 세포주 안정화, 항체/프라이머 예비검증(양성대조군으로 확인) |
| 3-4주 | 실험 1(발현 프로파일링) + 실험 2(세포종류 대조, 우선 진행) 수행 |
| 5-6주 | 실험 3(Loss-of-function 표현형: 마커/증식/이동능) 수행 |
| 7주 | 실험 4(Gain-of-function) 수행 |
| 8주 | 데이터 정리·통계·그림화, Phase 2 진행 여부(ARID1A/THAP11 중 우선순위) 결정 |

### Phase 1 예산 (개략 추정 — 실제 벤더 견적/원내 코어퍼실리티 단가로 재확인 필요)

| 항목 | 세부 | 추정 비용 (KRW) |
|---|---|---|
| 세포주 | LX-2, AML12 구입(원내 미보유 시) | 60만-100만 |
| 세포배양 소모품 | 배지, FBS, 플레이트/플라스크, 6-8주분 (2세포주) | 150만-250만 |
| siRNA + 형질주입시약 | ARID1A/THAP11/scrambled 각 2-3종 + Lipofectamine RNAiMAX 등 | 100만-180만 |
| 과발현 플라스미드 | ARID1A/THAP11 ORF 클론(Addgene/OriGene) + 미디프렙 | 80만-150만 |
| qPCR 시약 | 프라이머(타겟 2종 + 마커 8-10종) + master mix | 80만-120만 |
| Western blot 항체 | ARID1A, THAP11, α-SMA, COL1A1, GAPDH 등 1차 항체 6-8종 (보통 가장 비싼 항목) | 300만-500만 |
| 증식/이동능 assay 키트 | CCK-8 또는 BrdU, transwell insert | 40만-70만 |
| 일반 소모품 | 튜브, 팁, 기타 | 50만-80만 |
| **합계 (2타겟 병행)** | | **약 860만-1,450만원** |

- 위 표는 국내 통상 시약 단가 기준 개략 추정치이며, 실제로는 Bioneer(국내 프라이머/siRNA 합성), Santa Cruz/CST(항체), Thermo/Corning(배양 소모품) 등 실제 견적을 받아 확정해야 함.
- 세포주·기본 장비(세포배양기, qPCR기, Western 장비)가 이미 랩실에 있다는 전제 하의 소모품 위주 추정치 — 장비 신규 구매는 포함 안 됨.
- 항체 비용이 전체 예산의 1/3 이상을 차지 — 예산이 빠듯하면 우선 ARID1A/THAP11 항체만 먼저 검증하고 마커 항체는 가장 핵심적인 것(α-SMA, COL1A1)부터 시작하는 식으로 단계적 축소 가능.

## Phase 2 — 복합체/파트너 수준 기전 검증 (약 2-3개월, Phase 1 결과로 ARID1A/THAP11 중 우선순위 결정 후 진행)

Phase 1 결과에서 어느 쪽이 더 강한/일관된 표현형을 보이는지에 따라 아래 중 하나(또는 둘 다 자원이 되면 병행)로 진행. 두 타겟은 파트너가 다르므로 실험 설계도 다름:

5. **[ARID1A] Epistasis 실험**: ARID1A 단독 KD, SMARCA4 단독 KD, 이중 KD를 비교 — 같은 복합체 안에서 작동한다면 이중 KD가 단독 KD보다 크게 강하지 않아야 함(non-additive). "우연히 같이 뜬 유전자"와 "진짜 같은 기능 단위" 가설을 구분하는 핵심 실험.
6. **[ARID1A] 크로마틴 접근성**: ARID1A KD 전후로 SMARCA4의 알려진 타겟 유전자좌(IGFBP5 등, 2023 논문에서 확인된 좌위) ATAC-seq 또는 ChIP-qPCR(H3K27ac, BRG1 occupancy) — ARID1A가 SMARCA4-의존적 크로마틴 개방에 필요한지 직접 검증.
7. **[THAP11] 파트너 관계 검증**: THAP11은 알려진 결합 파트너가 HCF-1(HCFC1)이다 — 저희 consensus 리스트에도 HCFC1이 Tier 3로 들어가 있어(druggability 낮아 우선순위는 낮게 뒀지만) 같은 논리로 재검토 가치가 있음. THAP11-HCF-1 이중 KD epistasis, HCF-1 co-IP로 HSC에서도 상호작용이 유지되는지 확인.

## Phase 3 — 랩실 고유 독성물질 축과 연결 (Phase 1과 병행 가능)

8. TGF-β1 대신(또는 추가로) **CCl4/TAA 활성대사물 또는 처리된 간세포 조건배지**로 LX-2/AML12 자극 → ARID1A/THAP11 유도 여부 확인. 랩실이 이미 다루는 PHMG-P/CMIT-MIT를 간세포-HSC 공동배양에 적용해서 "폐 섬유화 유발 물질이 간 HSC의 이 축도 건드리는가"까지 확장하면 랩실 정체성과 직접 연결됨.

## Phase 4 — In vivo 확인 (약 6개월 이상, 자원 허용 시 — Phase 1/2에서 우선순위가 정해진 타겟 1개로 진행 권장)

9. **HSC-특이적 조건부 KO**: Lrat-Cre(이미 공개/확보 가능, Friedman lab 유래) × Arid1a-flox(Jackson Labs에서 확보 가능) 교배 — 새로 만들 필요 없이 기존 마우스 라인 조합으로 가능(THAP11-flox는 별도 확보 필요 여부 확인 필요). CCl4 또는 TAA로 섬유화 유도(저희가 이미 bioinformatics로 검증한 두 모델 그대로 재사용) 후 Sirius red, hydroxyproline, α-SMA IHC로 섬유화 정량, WT littermate와 비교.
10. **(장기) 약물 재창출 검증**: ARID1A 경로면 FHD-286류 BAF 복합체 억제제(임상 보류 이력이 있어 직접 사용은 제약이 있을 수 있음 — 접근 가능한 academic tool compound 대안 확인 필요)로 CCl4/TAA 모델에서 섬유화 완화 여부 테스트 — 랩실이 원래 하던 "독성 완화 물질 발굴" 논문 구조와 정확히 맞아떨어짐. THAP11 경로면 이 단계 전에 억제제 자체를 새로 찾아야 함(더 이른 단계의 신약개발 과제).

## 실현 가능성 메모
- LX-2, AML12 세포주와 ARID1A/THAP11 대상 siRNA·항체는 둘 다 표준적으로 구하기 쉬움 — Phase 1은 일반적인 약리학 랩 셋업으로 바로 시작 가능.
- THAP11은 문헌이 거의 없는 만큼 검증된 항체/siRNA 종류 자체가 ARID1A보다 적을 수 있음 — 주문 전 Cell Signaling/Sigma 카탈로그에서 실제 재고·검증 데이터(WB validated 여부) 확인 필요 (1주차 일정에 이미 포함).
- Lrat-Cre, Arid1a-flox 마우스는 기존 공개 라인이라 새로 만들 필요 없음 — Phase 4가 ARID1A 경로로 가면 진입장벽이 낮음. THAP11-flox 마우스는 존재 여부를 별도 확인해야 함 (없으면 Phase 4는 사실상 ARID1A 전용이 되고, THAP11은 Phase 1-2 in vitro 수준에서 논문 기여를 맺는 방향이 현실적).
- ATAC-seq(Phase 2)은 core facility 의뢰가 필요할 수 있음 — 예산/일정에 따라 ChIP-qPCR로 축소 가능.
- 예산표는 Phase 1(2타겟 병행 in vitro)만 포함 — Phase 2 이후(ATAC-seq, 동물실험)는 별도 예산 산정 필요하며 Phase 1 결과를 보고 나서 산정하는 게 낭비가 적음.
