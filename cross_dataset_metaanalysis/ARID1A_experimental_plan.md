# ARID1A 실험 검증 계획

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

## Phase 1 — In vitro 발현/기능 확인 (약 4-6주, 우선순위 1)

1. **발현 프로파일링**: LX-2(인체 HSC 세포주, 표준적으로 가장 많이 쓰임) 또는 1차 마우스 HSC를 quiescent → TGF-β1 자극 활성화 시간경과로 배양, ARID1A mRNA/단백질(qPCR, Western)이 활성화와 함께 증가하는지 확인 — 저희 bioinformatics 예측(진행에 따라 증가) 검증.
2. **Loss-of-function**: siRNA 또는 CRISPRi로 ARID1A knockdown → TGF-β1 자극 유무 각각에서 섬유화 마커(ACTA2/α-SMA, COL1A1, TIMP1, SERPINE1) mRNA+단백질, 증식(CCK-8/BrdU), 이동능(transwell) 측정.
   - *SMARCA4 2023 논문과 동일한 readout 패널을 그대로 사용 — 두 결과를 나란히 비교할 수 있어 "같은 복합체, 같은 표현형"이라는 주장이 훨씬 설득력 있어짐.*
3. **Gain-of-function**: ARID1A 과발현만으로 quiescent LX-2에서 활성화 마커가 유도되는지 확인 (충분조건 테스트).
4. **세포종류 대조 실험 (신규, 우선순위 최상)**: 같은 ARID1A knockdown을 HSC(LX-2)와 간세포주(AML12 또는 1차 마우스 간세포) 양쪽에서 동시에 수행하고 지방증/손상 마커(간세포 쪽은 Hepatology 2019 논문의 지표인 지방생성/지방산산화 유전자, HSC 쪽은 위 2번의 섬유화 마커)를 비교. 기존 문헌(간세포=보호적)과 저희 가설(HSC=촉진적)이 실제로 반대 방향으로 나오는지가 이 프로젝트 전체 스토리의 성패를 가르는 실험 — Phase 1의 다른 실험보다 먼저 하는 게 좋음.

## Phase 2 — 복합체 수준 기전 검증 (약 2-3개월)

4. **Epistasis 실험**: ARID1A 단독 KD, SMARCA4 단독 KD, 이중 KD를 비교 — 같은 복합체 안에서 작동한다면 이중 KD가 단독 KD보다 크게 강하지 않아야 함(non-additive). 이게 "우연히 같이 뜬 유전자"와 "진짜 같은 기능 단위" 가설을 구분하는 핵심 실험.
5. **크로마틴 접근성**: ARID1A KD 전후로 SMARCA4의 알려진 타겟 유전자좌(IGFBP5 등, 2023 논문에서 확인된 좌위) ATAC-seq 또는 ChIP-qPCR(H3K27ac, BRG1 occupancy) — ARID1A가 SMARCA4-의존적 크로마틴 개방에 필요한지 직접 검증.

## Phase 3 — 랩실 고유 독성물질 축과 연결 (Phase 1과 병행 가능)

6. TGF-β1 대신(또는 추가로) **CCl4/TAA 활성대사물 또는 처리된 간세포 조건배지**로 LX-2/1차 HSC 자극 → ARID1A 유도 여부 확인. 랩실이 이미 다루는 PHMG-P/CMIT-MIT를 간세포-HSC 공동배양에 적용해서 "폐 섬유화 유발 물질이 간 HSC의 ARID1A 축도 건드리는가"까지 확장하면 랩실 정체성과 직접 연결됨.

## Phase 4 — In vivo 확인 (약 6개월 이상, 자원 허용 시)

7. **HSC-특이적 조건부 KO**: Lrat-Cre(이미 공개/확보 가능, Friedman lab 유래) × Arid1a-flox(Jackson Labs에서 확보 가능) 교배 — 새로 만들 필요 없이 기존 마우스 라인 조합으로 가능. CCl4 또는 TAA로 섬유화 유도(저희가 이미 bioinformatics로 검증한 두 모델 그대로 재사용) 후 Sirius red, hydroxyproline, α-SMA IHC로 섬유화 정량, WT littermate와 비교.
8. **(장기) 약물 재창출 검증**: FHD-286류 BAF 복합체 억제제(임상 보류 이력이 있어 직접 사용은 제약이 있을 수 있음 — 접근 가능한 academic tool compound 대안 확인 필요)로 CCl4/TAA 모델에서 섬유화 완화 여부 테스트 — 랩실이 원래 하던 "독성 완화 물질 발굴" 논문 구조와 정확히 맞아떨어짐.

## 실현 가능성 메모
- LX-2 세포주, siRNA/CRISPRi 시약은 표준적으로 구하기 쉬움 — Phase 1은 일반적인 약리학 랩 셋업으로 바로 시작 가능. AML12(마우스 간세포주)도 마찬가지로 표준적으로 구하기 쉬워 세포종류 대조 실험 추가 진입장벽은 낮음.
- Lrat-Cre, Arid1a-flox 마우스 모두 기존에 공개된 라인이라 새로 만들 필요 없음 — Phase 4의 진입장벽이 생각보다 낮음.
- ATAC-seq(Phase 2-5)은 core facility 의뢰가 필요할 수 있음 — 예산/일정에 따라 ChIP-qPCR로 축소 가능.
