# 논문용 ggplot2 테마/팔레트/저장 규격을 한 곳에 고정.
# 모든 데이터셋 스크립트는 이 파일만 source()해서 시각화 스타일을 통일한다.

library(ggplot2)

theme_pub = function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title = element_text(face = "bold", size = base_size + 1),
      legend.position = "right",
      axis.title = element_text(face = "bold")
    )
}

# Okabe-Ito colorblind-safe qualitative palette, 질환 진행단계(4단계) 순서에 맞춰 배정.
# 랩실 스크립트에서 쓰던 group 이름(Normal/Steatosis_Only/Early_NASH_Fibrosis/Advanced_Fibrosis)과
# 시계열 timepoint 이름 둘 다에 재사용할 수 있도록 이름만 바꿔 끼우면 됨.
pub_palette_progression = c(
  "Normal"               = "#0072B2",
  "Steatosis_Only"       = "#009E73",
  "Early_NASH_Fibrosis"  = "#E69F00",
  "Advanced_Fibrosis"    = "#D55E00"
)

# Heatmap용 diverging color ramp (z-score 발현/활성도 표시용). red-blue 축은 색맹 친화적인
# RColorBrewer "RdBu"를 반전(양수=red, 음수=blue)해서 사용.
pub_heatmap_colors = function(n = 100) {
  rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(n))
}

# 모든 최종 figure는 이 함수로 저장 -> 크기/해상도/포맷(PDF, 벡터)이 스크립트마다 달라지는 것을 방지.
ggsave_pub = function(filename, plot = ggplot2::last_plot(), width = 7, height = 5) {
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height,
                   units = "in", device = "pdf")
}
