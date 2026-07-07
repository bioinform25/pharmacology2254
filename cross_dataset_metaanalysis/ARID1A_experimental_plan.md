# ARID1A 실험 검증 계획

## 배경 가설
ARID1A는 SMARCA4(BRG1)와 함께 SWI/SNF(BAF) 크로마틴 리모델링 복합체를 구성한다. SMARCA4는 이미 간 HSC-근섬유아세포 전환의 직접적 드라이버로 검증되어 있다(Li et al., Cell Death Dis 2023 — Lrat-Cre;Smarca4^fl/fl 마우스에서 HSC-특이적 결손이 여러 CCl4/BDL 모델에서 섬유화를 완화). ARID1A는 간암(간세포암) 맥락에서는 연구됐지만 **간섬유화/HSC 맥락에서는 문헌이 전혀 없음**을 확인했다. 저희 3개 독립 데이터셋(인체 MASLD, 마우스 CCl4, 마우스 TAA) 메타분석에서 ARID1A가 세 축 모두 유의 + 방향 일치로 나온 것이, 같은 복합체의 파트너가 함께 움직인다는 가설과 일치한다.

**핵심 가설**: ARID1A는 SMARCA4와 함께(같은 BAF 복합체 안에서) HSC 활성화에 필요하며, 독성물질(CCl4/TAA) 및 대사성 자극 양쪽에서 공통으로 유도된다.

## Phase 1 — In vitro 발현/기능 확인 (약 4-6주, 우선순위 1)

1. **발현 프로파일링**: LX-2(인체 HSC 세포주, 표준적으로 가장 많이 쓰임) 또는 1차 마우스 HSC를 quiescent → TGF-β1 자극 활성화 시간경과로 배양, ARID1A mRNA/단백질(qPCR, Western)이 활성화와 함께 증가하는지 확인 — 저희 bioinformatics 예측(진행에 따라 증가) 검증.
2. **Loss-of-function**: siRNA 또는 CRISPRi로 ARID1A knockdown → TGF-β1 자극 유무 각각에서 섬유화 마커(ACTA2/α-SMA, COL1A1, TIMP1, SERPINE1) mRNA+단백질, 증식(CCK-8/BrdU), 이동능(transwell) 측정.
   - *SMARCA4 2023 논문과 동일한 readout 패널을 그대로 사용 — 두 결과를 나란히 비교할 수 있어 "같은 복합체, 같은 표현형"이라는 주장이 훨씬 설득력 있어짐.*
3. **Gain-of-function**: ARID1A 과발현만으로 quiescent LX-2에서 활성화 마커가 유도되는지 확인 (충분조건 테스트).

## Phase 2 — 복합체 수준 기전 검증 (약 2-3개월)

4. **Epistasis 실험**: ARID1A 단독 KD, SMARCA4 단독 KD, 이중 KD를 비교 — 같은 복합체 안에서 작동한다면 이중 KD가 단독 KD보다 크게 강하지 않아야 함(non-additive). 이게 "우연히 같이 뜬 유전자"와 "진짜 같은 기능 단위" 가설을 구분하는 핵심 실험.
5. **크로마틴 접근성**: ARID1A KD 전후로 SMARCA4의 알려진 타겟 유전자좌(IGFBP5 등, 2023 논문에서 확인된 좌위) ATAC-seq 또는 ChIP-qPCR(H3K27ac, BRG1 occupancy) — ARID1A가 SMARCA4-의존적 크로마틴 개방에 필요한지 직접 검증.

## Phase 3 — 랩실 고유 독성물질 축과 연결 (Phase 1과 병행 가능)

6. TGF-β1 대신(또는 추가로) **CCl4/TAA 활성대사물 또는 처리된 간세포 조건배지**로 LX-2/1차 HSC 자극 → ARID1A 유도 여부 확인. 랩실이 이미 다루는 PHMG-P/CMIT-MIT를 간세포-HSC 공동배양에 적용해서 "폐 섬유화 유발 물질이 간 HSC의 ARID1A 축도 건드리는가"까지 확장하면 랩실 정체성과 직접 연결됨.

## Phase 4 — In vivo 확인 (약 6개월 이상, 자원 허용 시)

7. **HSC-특이적 조건부 KO**: Lrat-Cre(이미 공개/확보 가능, Friedman lab 유래) × Arid1a-flox(Jackson Labs에서 확보 가능) 교배 — 새로 만들 필요 없이 기존 마우스 라인 조합으로 가능. CCl4 또는 TAA로 섬유화 유도(저희가 이미 bioinformatics로 검증한 두 모델 그대로 재사용) 후 Sirius red, hydroxyproline, α-SMA IHC로 섬유화 정량, WT littermate와 비교.
8. **(장기) 약물 재창출 검증**: FHD-286류 BAF 복합체 억제제(임상 보류 이력이 있어 직접 사용은 제약이 있을 수 있음 — 접근 가능한 academic tool compound 대안 확인 필요)로 CCl4/TAA 모델에서 섬유화 완화 여부 테스트 — 랩실이 원래 하던 "독성 완화 물질 발굴" 논문 구조와 정확히 맞아떨어짐.

## 실현 가능성 메모
- LX-2 세포주, siRNA/CRISPRi 시약은 표준적으로 구하기 쉬움 — Phase 1은 일반적인 약리학 랩 셋업으로 바로 시작 가능.
- Lrat-Cre, Arid1a-flox 마우스 모두 기존에 공개된 라인이라 새로 만들 필요 없음 — Phase 4의 진입장벽이 생각보다 낮음.
- ATAC-seq(Phase 2-5)은 core facility 의뢰가 필요할 수 있음 — 예산/일정에 따라 ChIP-qPCR로 축소 가능.
