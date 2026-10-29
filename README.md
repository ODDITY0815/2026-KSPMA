# 규준 공간 기반 ADHD 동공 반응 분석


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

## 분석 구성

1. **Helper functions**: 외부 중심화 FACE, 개인 평균, 사다리꼴 가중치, 관측 격자 투영 및 PVE 계산.
2. **곡선 QC**: `exclude_subject_ids: ["45"]`를 전체 분석에 적용. 원래 격자 coverage ≥ .50, load 1/2. 분석 격자 100 Hz(10개 중 1개).
3. **MFPCA**: 45번을 제외한 대조군 21명으로 규준평균·고유값·고유함수를 새로 추정. 개인 동일가중 평균, `weight="subj"`, `twoway=FALSE`, `p=3`, `m=2`, 35 knots. 최대 수준별 10개 중 CVE 80%에 처음 도달하는 성분 수 선택.
4. **투영·통계**: BPC는 세션 평균−규준평균의 전체 격자 내적. WPC는 시행−실제 세션 평균 잔차의 관측 격자 WLS(ridge 1e-8). 세션별 WPC 표본분산은 n−1 분모 사용.
5. **시각화**: 규준 평균·CVE·BPC·WPC 4패널, 성분별 고유함수, BPC2·WPC1·WPC2 평균 ±1SD, 상자 그림, 행동 상관 히트맵.

45번은 규준평균 계산, 공분산·고유함수 추정, 투영, 집단 비교, 행동 상관에서 모두 제외됩니다. 

집단 비교는 ML LMM의 집단효과 LRT입니다. BPC는 점수, WPC는 자연로그 분산을 검정하며 BH 보정은 수준별로 적용합니다. 사후비교는 같은 LMM의 Wald z와 성분 내 세 쌍 Holm 보정입니다. 전체효과 q<.05인 성분에만 브라켓을 표시합니다. 개별 대비 CI는 보정 전 95% Wald 구간입니다. On−Off 대응 t검정을 별도 민감도 분석으로 출력합니다.


$$
Y_{ij} = \beta_0 + \beta_1 D_{\mathrm{On},ij} + \beta_2 D_{\mathrm{Off},ij} + u_i + \epsilon_{ij}
$$


- $\beta_0$: 정상 대조군(Control)의 평균
- $\beta_1$: ADHD-On과 Control의 평균 차이
- $\beta_2$: ADHD-Off와 Control의 평균 차이
- $u_i$: 참가자별 무작위 절편
- $\epsilon_{ij}$: 잔차

행동지표는 동공 QC 이전의 유효 반응 정답률과 정답·200–1500 ms 반응시간의 표준편차(최소 10개)입니다. 참가자×집단 BPC 평균 및 WPC 원척도 분산과 집단별 Pearson 상관을 계산하며 전체 상관검정에 BH 보정을 적용합니다. 


## 출력


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



## 실행 결과

대조군 규준 학습 표본은 21명·3,102시행이며, 비교 표본은 Control 21 / Off 28 / On 17의 총 66세션입니다. CVE 80% 기준으로 BPC 4개(83.32%), WPC 5개(81.46%)가 선택되었습니다.

| 성분 | LRT χ²(2) | BH q |
|---|---:|---:|
| BPC2 | 12.581049 | .007415 |
| WPC1 | 16.638938 | .000609 |
| WPC2 | 28.492631 | .00000325 |
| WPC5 | 7.297865 | .043365 |

WPC5는 singular fit입니다. 고유함수의 부호는 임의이며, 현재 재적합된 BPC2는 과제 초반 양수·Probe 구간 음수 방향입니다. 점수와 상관 부호는 해당 고유함수의 방향과 함께 해석합니다.


## References

Aston-Jones, G., & Cohen, J. D. (2005). An integrative theory of locus coeruleus-norepinephrine function: Adaptive gain and optimal performance. *Annual Review of Neuroscience, 28*(1), 403–450.

Cui, E., Li, R., Crainiceanu, C. M., & Xiao, L. (2023). Fast multilevel functional principal component analysis. *Journal of Computational and Graphical Statistics, 32*(2), 366–377.

Di, C. Z., Crainiceanu, C. M., Caffo, B. S., & Punjabi, N. M. (2009). Multilevel functional principal component analysis. *The Annals of Applied Statistics, 3*(1), 458.

Wainstein, G., Rojas-Líbano, D., Crossley, N. A., Carrasco, X., Aboitiz, F., & Ossandón, T. (2017). Pupil size tracks attentional performance in attention-deficit/hyperactivity disorder. *Scientific Reports, 7*(1), 8228.

Yadav, S., Aadya, & Shukla, J. (2026). Trial-level pupil dynamics reflect attention-related variability and medication-associated modulation in ADHD. *Scientific Reports*.
