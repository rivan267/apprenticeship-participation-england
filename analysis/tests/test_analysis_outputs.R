# ONE test: the analysis pipeline produces its core outputs with internally
# consistent numbers. RED until analysis/03_analysis.R exists and runs clean.
test_analysis_outputs <- function() {
  out <- "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England/analysis/out"
  t1p <- file.path(out, "tab1_starts_age_level.csv")
  s2p <- file.path(out, "tab2_segmented.csv")
  stopifnot(file.exists(t1p))
  t1 <- read.csv(t1p, check.names = FALSE)
  stopifnot(all(2014:2024 %in% t1$year))                  # full decade present
  cellcols <- setdiff(names(t1), c("year", "year_label"))
  # 3 ages x 4 levels: Intermediate/Advanced/Higher plus a residual "Unknown"
  # level that the DfE file carries (80 starts in 2014/15, zero after)
  stopifnot(length(cellcols) == 12)
  stopifnot(file.exists(s2p))
  s2 <- read.csv(s2p)
  stopifnot(nrow(s2) == 7)                                 # Total + 3 ages + 3 levels
  cat("PASS: analysis outputs exist, cover 2014/15-2024/25, 12 cells, 7 segmented series\n")
}
test_analysis_outputs()

# Next behaviour: RQ2 deprivation gradient table covering the full decade.
test_rq2_gradient <- function() {
  out <- "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England/analysis/out"
  rq2 <- read.csv(file.path(out, "tab_rq2_gradient.csv"))
  stopifnot(all(2014:2024 %in% rq2$year))
  stopifnot(all(c("Q1", "Q5", "ratio_Q1_Q5", "conc_index") %in% names(rq2)))
  stopifnot(all(rq2$Q1 > 0) && all(rq2$Q5 > 0))
  cat("PASS: RQ2 gradient table covers 2014/15-2024/25 with Q1, Q5, ratio, CI\n")
}
test_rq2_gradient()

# Next behaviour: RQ3 regional starts counts for all nine English regions, full decade.
test_rq3_regions <- function() {
  out <- "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England/analysis/out"
  reg <- read.csv(file.path(out, "tab_region_counts.csv"))
  stopifnot(length(unique(reg$region)) == 9)
  stopifnot(all(2014:2024 %in% reg$year))
  stopifnot(all(reg$starts > 0, na.rm = TRUE))
  cat("PASS: RQ3 regional counts cover 9 regions across 2014/15-2024/25\n")
}
test_rq3_regions()

# Next behaviour: the five planned manuscript figures exist as 300 dpi PNGs.
test_figures <- function() {
  out <- "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England/analysis/out"
  figs <- file.path(out, sprintf("fig%d_%s.png", 1:5,
          c("total_starts", "age_composition", "level", "deprivation", "region_slope")))
  stopifnot(all(file.exists(figs)))
  stopifnot(all(file.size(figs) > 20000))   # non-trivial rendered images
  cat("PASS: five manuscript figures rendered\n")
}
test_figures()

# Next behaviour: RQ1 cell decomposition of the decade change sums to the total change.
test_decomposition <- function() {
  out <- "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England/analysis/out"
  dd <- read.csv(file.path(out, "tab_decomp_change.csv"))
  stopifnot(all(c("age_summary", "apps_level", "s0", "s1", "change", "pct_of_total_change") %in% names(dd)))
  t1 <- read.csv(file.path(out, "tab1_starts_age_level.csv"), check.names = FALSE)
  cellcols <- setdiff(names(t1), c("year", "year_label"))
  tot_chg <- sum(t1[t1$year == 2024, cellcols], na.rm = TRUE) - sum(t1[t1$year == 2014, cellcols], na.rm = TRUE)
  stopifnot(abs(sum(dd$change) - tot_chg) <= 100)   # cells reproduce the decade change
  cat("PASS: RQ1 decomposition sums to the decade change\n")
}
test_decomposition()
