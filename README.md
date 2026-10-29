# 규준 공간 기반 ADHD 동공 반응 분석

`poster_final_analysis.Rmd`는 45번 참가자를 **규준 공간 구축부터 모든 분석에서 제외**하고 동공 반응의 MFPCA, 규준 공간 투영, 집단 비교, 행동 상관 및 시각화를 실행합니다. Rmd는 절 제목과 실행 가능한 코드블록으로 구성되어 있으며 주요 표와 그림을 직접 출력합니다.

## 실행

필요 패키지: `refund`, `data.table`, `lme4`, `ggplot2`, `patchwork`, `knitr`, `rmarkdown`. HTML 출력에는 Pandoc이 필요합니다(RStudio 포함).

```r
install.packages(c("refund", "data.table", "lme4", "ggplot2", "patchwork", "knitr", "rmarkdown"))
rmarkdown::render("poster_final_analysis.Rmd")
```

저장소 루트에서 실행하거나 RStudio에서 Rmd를 열어 Knit합니다. 수동 chunk 실행 시 작업 디렉터리는 Rmd가 있는 폴더로 설정합니다. 입력 경로는 Rmd 기준 상대경로이며, 별도 저장소에서는 다음과 같이 변경할 수 있습니다.

```r
rmarkdown::render(
  "poster_final_analysis.Rmd",
  params = list(
    data_file = "/absolute/path/pupil_curves_full.rds",
    timing_file = "/absolute/path/event_timing_trials.csv",
    output_dir = "outputs_reference_exclude45"
  )
)
```

## 입력

| RDS 원소 | 구조 |
|---|---|
| `Y_original` | 전처리·시행 내 표준화된 시행 × 시간 행렬 |
| `argvals` | 초 단위 시간축, 행렬 열에 대응 |
| `curve_meta` | 행렬 행과 같은 순서의 시행 메타데이터 |

메타데이터 필수 열은 `subject_id`, `session_id`, `group`, `coverage`, `load`, `performance`, `reaction_time`입니다. 집단 코드는 `Ctrl`, `off-ADHD`, `on-ADHD`, 반응시간 단위는 ms입니다. 같은 참가자의 On/Off는 같은 `subject_id`를 사용합니다.

기본 입력은 `../pupil_R_export/pupil_curves_full.rds`입니다. 추가 표준화·보간·결측 대치를 수행하지 않습니다. 이벤트 표 `../pupil_poster/results/event_timing_trials.csv`는 그림의 시간 구간 표시용이며, 없으면 명목 구간을 사용합니다.

## 분석 구성

1. **Helper functions**: 외부 중심화 FACE, 개인 평균, 사다리꼴 가중치, 관측 격자 투영 및 PVE 계산.
2. **곡선 QC**: `exclude_subject_ids: ["45"]`를 전체 분석에 적용. 원래 격자 coverage ≥ .50, load 1/2. 분석 격자 100 Hz(10개 중 1개).
3. **MFPCA**: 45번을 제외한 대조군 21명으로 규준평균·고유값·고유함수를 새로 추정. 개인 동일가중 평균, `weight="subj"`, `twoway=FALSE`, `p=3`, `m=2`, 35 knots. 최대 수준별 10개 중 CVE 80%에 처음 도달하는 성분 수 선택.
4. **투영·통계**: BPC는 세션 평균−규준평균의 전체 격자 내적. WPC는 시행−실제 세션 평균 잔차의 관측 격자 WLS(ridge 1e-8). 세션별 WPC 표본분산은 n−1 분모 사용.
5. **시각화**: 규준 평균·CVE·BPC·WPC 4패널, 성분별 고유함수, BPC2·WPC1·WPC2 평균 ±1SD, 상자 그림, 행동 상관 히트맵.

45번은 규준평균 계산, 공분산·고유함수 추정, 투영, 집단 비교, 행동 상관에서 모두 제외됩니다. QC 기록에는 제외 여부가 남습니다. 새 규준 기저로 모든 점수와 통계를 다시 계산하며 이전 모형이나 점수를 불러오지 않습니다.

집단 비교는 `test_value ~ group + (1 | subject_id)` ML LMM의 집단효과 LRT입니다. BPC는 점수, WPC는 자연로그 분산을 검정하며 BH 보정은 수준별로 적용합니다. 사후비교는 같은 LMM의 Wald z와 성분 내 세 쌍 Holm 보정입니다. 전체효과 q<.05인 성분에만 브라켓을 표시합니다. 개별 대비 CI는 보정 전 95% Wald 구간입니다. On−Off 대응 t검정을 별도 민감도 분석으로 출력합니다.

행동지표는 동공 QC 이전의 유효 반응 정답률과 정답·200–1500 ms 반응시간의 표준편차(최소 10개)입니다. 참가자×집단 BPC 평균 및 WPC 원척도 분산과 집단별 Pearson 상관을 계산하며 전체 상관검정에 BH 보정을 적용합니다. 제외 참가자는 행동지표에서도 제거합니다.

## 그림 디자인

기존 `poster_session_centered_plots.Rmd`, `session_centered_lmm_boxplots.R`, `correlation_heatmap_wide.Rmd`의 스타일을 적용합니다.

- BPC: `#17675F`, `#4F71BF`, `#77B6A8`, `#926BB0` 순서.
- WPC: `#963F35`, `#D18B29`, `#C07C62`, `#A44D88`, `#78639E` 순서.
- 집단: Control `#69A66F`, Off `#D87373`, On `#7199C6`.
- 함수형 그림: 0–8000 ms, 상단 이벤트 레이블·회색 자극 음영, 개별 패널 6.4 × 4 inch.
- ±1SD: 평균 검정 실선, 해당 성분색의 +SD 실선·−SD 점선. 고유값·고유함수 쌍으로 계산.
- 상자 그림: BPC/WPC 한 행 배치, 수준별 공통 y축, WPC 분산 로그축.
- 히트맵: 행동지표별 두 행 패널, Pearson r −1~1의 파랑–흰색–빨강, 상관계수 숫자 표시.

## 출력

현재 결과는 `outputs_reference_exclude45/`에 저장됩니다. 기존 `outputs/`는 이전 실행 결과이며 이번 분석에서 사용하지 않습니다.

| 출력 | 내용 |
|---|---|
| `reference_participants.csv`, `reference.rds` | 규준 학습 대상 및 새 모형·기저·평균 |
| `qc_flow.csv`, `curve_qc.csv`, `session_qc.csv`, `projection_qc.csv` | 단계별 QC |
| `sample_summary.csv`, `pve.csv` | 표본 수·시행 수, 고유값·PVE·CVE |
| `between_scores.csv`, `within_scores.csv`, `within_session_variance.csv` | 점수·분산 |
| `lmm_omnibus.csv`, `lmm_posthoc.csv`, `paired_sensitivity.csv` | 집단 및 대응 분석 |
| `correlations.csv`, `behavior_participants.csv` | 행동 상관 및 분석 자료 |
| `01_normative_four_panels.png`, `normative_*.png`, `eigen_*.png` | 규준 공간 합본·개별 그림 |
| `02_mean_plus_minus_sd.png`, `sd_*_panel.png`, `sd_curves.csv` | ±1SD 그림·곡선 |
| `03_boxplots.png`, `boxplots_*_LMM.png` | 상자 그림 |
| `04_correlation_heatmap.png` | 상관 히트맵 |
| `analysis_parameters.rds`, `input_manifest.csv`, `sessionInfo.txt` | 설정·입력 MD5·실행환경 |

그림은 300 dpi PNG입니다. HTML 보고서는 `poster_final_analysis.html`로 생성됩니다. 이 저장소에는 Rmd와 README를 제공합니다. 원자료·개인별 결과·생성 HTML은 포함하지 않습니다. 로컬 실행 결과를 추가 커밋할 경우 `outputs*/`와 HTML을 제외하세요.

## 실행 결과

대조군 규준 학습 표본은 21명·3,102시행이며, 비교 표본은 Control 21 / Off 28 / On 17의 총 66세션입니다. CVE 80% 기준으로 BPC 4개(83.32%), WPC 5개(81.46%)가 선택되었습니다.

| 성분 | LRT χ²(2) | BH q |
|---|---:|---:|
| BPC2 | 12.581049 | .007415 |
| WPC1 | 16.638938 | .000609 |
| WPC2 | 28.492631 | .00000325 |
| WPC5 | 7.297865 | .043365 |

WPC5는 singular fit입니다. 고유함수의 부호는 임의이며, 현재 재적합된 BPC2는 과제 초반 양수·Probe 구간 음수 방향입니다. 점수와 상관 부호는 해당 고유함수의 방향과 함께 해석합니다.

2026-10-29 전체 Rmd 실행과 HTML 생성을 완료했습니다. 규준 학습 대상·투영점수·행동 분석에서 45번 제외, LMM 9개·사후비교 27개·상관 54개, 그림 출력을 검증했습니다. 실행환경의 전체 버전은 `outputs_reference_exclude45/sessionInfo.txt`에 기록됩니다.

## References

Aston-Jones, G., & Cohen, J. D. (2005). An integrative theory of locus coeruleus-norepinephrine function: Adaptive gain and optimal performance. *Annual Review of Neuroscience, 28*(1), 403–450.

Cui, E., Li, R., Crainiceanu, C. M., & Xiao, L. (2023). Fast multilevel functional principal component analysis. *Journal of Computational and Graphical Statistics, 32*(2), 366–377.

Di, C. Z., Crainiceanu, C. M., Caffo, B. S., & Punjabi, N. M. (2009). Multilevel functional principal component analysis. *The Annals of Applied Statistics, 3*(1), 458.

Wainstein, G., Rojas-Líbano, D., Crossley, N. A., Carrasco, X., Aboitiz, F., & Ossandón, T. (2017). Pupil size tracks attentional performance in attention-deficit/hyperactivity disorder. *Scientific Reports, 7*(1), 8228.

Yadav, S., Aadya, & Shukla, J. (2026). Trial-level pupil dynamics reflect attention-related variability and medication-associated modulation in ADHD. *Scientific Reports*.
