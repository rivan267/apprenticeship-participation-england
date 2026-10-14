# =====================================================================
# Who still starts an apprenticeship?  England, 2014/15 to 2024/25
# Minimal implementation for tests/test_analysis_outputs.R:
#   RED evidence (Rscript, 16 Jul 2026):
#     Error in test_analysis_outputs() : file.exists(t1p) is not TRUE
# This step implements ONLY what that one test requires:
#   A) national age x level table (tab1_starts_age_level.csv)
#   D) segmented trend regressions (tab2_segmented.csv)
# Data: DfE Explore Education Statistics underlying CSVs (OGL v3.0).
# =====================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("RPCE_BASE", unset = "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England")
DATA <- file.path(BASE, "data")
OUT  <- file.path(BASE, "analysis", "out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(x))))
yr  <- function(tp) as.integer(substr(gsub('"', "", as.character(tp)), 1, 4))
ay  <- function(y) sprintf("%d/%02d", y, (y + 1) %% 100)
log_chk <- file(file.path(OUT, "checks_overlap.txt"), "w")
chk <- function(...) { cat(..., "\n"); cat(..., "\n", file = log_chk) }

# ---- A) national age x level series (single consistent DfE file) ----
hs <- fread(file.path(DATA, "cur/data/app-historical-summary-to-2425.csv"), colClasses = "character")
hs[, year := yr(time_period)]
hs[, starts := num(starts)]
stopifnot(all(hs$geographic_level == "National"))
al <- hs[ssa_tier_1 == "Total", .(starts = sum(starts, na.rm = TRUE)), by = .(year, age_summary, apps_level)]
chk("age_summary values:", paste(sort(unique(al$age_summary)), collapse = " | "))
chk("apps_level values:",  paste(sort(unique(al$apps_level)),  collapse = " | "))

dec <- al[year >= 2014 & year <= 2024]
tot <- dec[age_summary == "Total" & apps_level == "Total"][order(year)]
chk("\nTotal starts 2014/15-2024/25:")
chk(paste(capture.output(print(tot[, .(year = ay(year), starts)])), collapse = "\n"))

t1 <- dcast(dec[age_summary != "Total" & apps_level != "Total"],
            year ~ age_summary + apps_level, value.var = "starts")
t1[, year_label := ay(year)]
fwrite(t1, file.path(OUT, "tab1_starts_age_level.csv"))

# integrity: cells must sum to published total each year
cell_sum <- dec[age_summary != "Total" & apps_level != "Total", .(s = sum(starts)), by = year]
mrg <- merge(cell_sum, tot[, .(year, starts)], by = "year")
chk("\nCell-sum vs published total by year (diff):")
chk(paste(capture.output(print(mrg[, .(year = ay(year), cells = s, total = starts, diff = s - starts)])), collapse = "\n"))

# ---- RQ1 decomposition of the decade change into age x level cells ----
# Drives test_decomposition (RED: cannot open '.../tab_decomp_change.csv')
cells <- dec[age_summary != "Total" & apps_level != "Total"]
d0 <- cells[year == 2014, .(age_summary, apps_level, s0 = starts)]
dl <- cells[year == 2024, .(age_summary, apps_level, s1 = starts)]
ddc <- merge(d0, dl, by = c("age_summary", "apps_level"), all = TRUE)
ddc[is.na(s0), s0 := 0]; ddc[is.na(s1), s1 := 0]
ddc[, change := s1 - s0]
tot_chg <- tot[year == 2024, starts] - tot[year == 2014, starts]
ddc[, pct_of_total_change := round(100 * change / tot_chg, 1)]
chk("\nRQ1 decomposition: cells sum", sum(ddc$change), "vs total change", tot_chg)
fwrite(ddc[order(change)], file.path(OUT, "tab_decomp_change.csv"))

# ---- D) segmented trends, knots 2017/18 and 2020/21 ----
seg_fit <- function(df, label) {
  df <- as.data.table(df)[order(year)]
  df[, t := year - 2014]
  # knots: slope changes apply FROM 2017/18 (t=3) and FROM 2020/21 (t=6),
  # i.e. the first post-levy and first pandemic years open their segments
  m <- lm(starts ~ t + pmax(t - 2, 0) + pmax(t - 5, 0), data = df)
  cf <- summary(m)$coefficients
  data.table(series = label,
             slope_pre_levy    = cf[2, 1], se_pre   = cf[2, 2],
             slope_change_levy = cf[3, 1], se_levy  = cf[3, 2],
             slope_change_covid = cf[4, 1], se_covid = cf[4, 2],
             slope_post_levy   = cf[2, 1] + cf[3, 1],
             slope_post_covid  = cf[2, 1] + cf[3, 1] + cf[4, 1],
             r2 = summary(m)$r.squared)
}
segs <- rbindlist(list(
  seg_fit(tot[, .(year, starts)], "Total"),
  seg_fit(dec[age_summary == "Under 19" & apps_level == "Total", .(year, starts)], "Under 19"),
  seg_fit(dec[age_summary == "19-24" & apps_level == "Total", .(year, starts)], "19-24"),
  seg_fit(dec[age_summary == "25+" & apps_level == "Total", .(year, starts)], "25+"),
  seg_fit(dec[apps_level == "Intermediate Apprenticeship" & age_summary == "Total", .(year, starts)], "Intermediate"),
  seg_fit(dec[apps_level == "Advanced Apprenticeship" & age_summary == "Total", .(year, starts)], "Advanced"),
  seg_fit(dec[apps_level == "Higher Apprenticeship" & age_summary == "Total", .(year, starts)], "Higher")))
num_cols <- names(segs)[sapply(segs, is.numeric)]
segs[, (num_cols) := lapply(.SD, function(x) round(x, 1)), .SDcols = num_cols]
chk("\nRQ4 segmented regression (annual starts, knots 2017/18 & 2020/21):")
chk(paste(capture.output(print(segs)), collapse = "\n"))
fwrite(segs, file.path(OUT, "tab2_segmented.csv"))

# ---- B) RQ2 deprivation 2014/15-2024/25 (three releases spliced) ----
# Drives tests/test_analysis_outputs.R::test_rq2_gradient
# (RED: cannot open file '.../out/tab_rq2_gradient.csv')
std_imd <- function(x) {
  x <- trimws(gsub('"', "", x))
  out <- rep(NA_character_, length(x))
  out[grepl("^(1|One)",   x, ignore.case = TRUE)] <- "Q1 most deprived"
  out[grepl("^(2|Two)",   x, ignore.case = TRUE)] <- "Q2"
  out[grepl("^(3|Three)", x, ignore.case = TRUE)] <- "Q3"
  out[grepl("^(4|Four)",  x, ignore.case = TRUE)] <- "Q4"
  out[grepl("^(5|Five)",  x, ignore.case = TRUE)] <- "Q5 least deprived"
  out[is.na(out)] <- x[is.na(out)]
  out
}
imd_levels <- c("Q1 most deprived", "Q2", "Q3", "Q4", "Q5 least deprived")

d1 <- fread(file.path(DATA, "hist1920/data/app-learner-deprivation-starts-1920-q4.csv"), colClasses = "character")
d1 <- d1[geographic_level == "National" & ssa_t1_desc == "Total" & age_youth_adult == "Total"]
d1[, `:=`(year = yr(time_period), starts = num(starts), imd = std_imd(learner_home_depriv),
          level = apps_level, age = age_summary, src = "rel1920")]

d2 <- fread(file.path(DATA, "hist/data/app-learner-deprivation-starts-202223-q4.csv"), colClasses = "character")
d2 <- d2[geographic_level == "National" & ssa_t1_desc == "Total" & age_youth_adult == "Total"]
d2[, `:=`(year = yr(time_period), starts = num(starts), imd = std_imd(learner_home_depriv),
          level = apps_level_detailed, age = age_summary, src = "rel2223")]

d3 <- fread(file.path(DATA, "cur/data/apps_11_learner_deprivation_202526_10.csv"), colClasses = "character")
d3 <- d3[geographic_level == "National" & ssa_tier_1 == "Total" & age_youth_adult == "Total"]
d3[, `:=`(year = yr(time_period), starts = num(starts), imd = std_imd(learner_home_depriv),
          level = apps_level, age = age_summary, src = "cur")]
d3 <- d3[year <= 2024]

dep <- rbindlist(list(
  d1[, .(year, imd, level, age, starts, src)],
  d2[, .(year, imd, level, age, starts, src)],
  d3[, .(year, imd, level, age, starts, src)]))
chk("\nIMD labels after standardisation:", paste(sort(unique(dep$imd)), collapse = " | "))
chk("dep age labels:", paste(sort(unique(dep$age)), collapse = " | "))

ovl <- dcast(dep[imd %in% imd_levels & age == "Total" & level == "Total"],
             year + imd ~ src, value.var = "starts", fun.aggregate = sum)
chk("\nDeprivation overlap (totals by quintile, by source), 2017-2022:")
chk(paste(capture.output(print(ovl[year %in% 2017:2022])), collapse = "\n"))

dep[, prio := fifelse(src == "cur", 3L, fifelse(src == "rel2223", 2L, 1L))]
setorder(dep, year, imd, level, age, prio)
depu <- dep[imd %in% imd_levels, .SD[.N], by = .(year, imd, level, age)]

g <- depu[age == "Total" & level == "Total"][order(year, match(imd, imd_levels))]
g[, share := 100 * starts / sum(starts), by = year]
rat <- dcast(g, year ~ imd, value.var = "starts")
setnames(rat, c("Q1 most deprived", "Q5 least deprived"), c("Q1", "Q5"))
rat[, ratio_Q1_Q5 := Q1 / Q5]
ci_fun <- function(s) { p <- s / sum(s); r <- seq(0.1, 0.9, 0.2); 2 * sum(p * r) - 1 }
ci <- g[, .(conc_index = ci_fun(starts)), by = year]
rq2 <- merge(rat[, .(year, Q1, Q5, ratio_Q1_Q5)], ci, by = "year")
chk("\nRQ2 gradient (Q1, Q5, ratio, concentration index):")
chk(paste(capture.output(print(rq2[, .(year = ay(year), Q1, Q5, ratio = round(ratio_Q1_Q5, 3), CI = round(conc_index, 4))])), collapse = "\n"))
fwrite(g,   file.path(OUT, "dep_totals_by_quintile.csv"))
fwrite(rq2, file.path(OUT, "tab_rq2_gradient.csv"))

wg_age <- depu[level == "Total" & age != "Total"]
wg_age[, share := 100 * starts / sum(starts), by = .(year, age)]
fwrite(wg_age, file.path(OUT, "dep_within_age.csv"))
wg_lv <- depu[age == "Total" & level != "Total"]
wg_lv[, share := 100 * starts / sum(starts), by = .(year, level)]
fwrite(wg_lv, file.path(OUT, "dep_within_level.csv"))

# ---- C) RQ3 regions: counts 2014/15-2024/25, rates where published ----
# Drives tests/test_analysis_outputs.R::test_rq3_regions
# (RED: cannot open file '.../out/tab_region_counts.csv')
r1 <- fread(file.path(DATA, "hist/data/app-starts-since-202223-q4.csv"), colClasses = "character")
r1 <- r1[geographic_level == "Regional"]
r1[, `:=`(year = yr(time_period), starts = num(starts), region = gsub('"', "", region_name), src = "since")]

r2 <- fread(file.path(DATA, "cur/data/app-geography-population-202526-q3.csv"), colClasses = "character")
r2 <- r2[geographic_level == "Regional" & apps_level == "Total" & age_summary == "Total"]
r2[, `:=`(year = yr(time_period), starts = num(starts), region = gsub('"', "", region_name),
          rate = num(starts_rate_per_100000_population), src = "curgeo")]
r2 <- r2[year <= 2024]

reg <- rbindlist(list(r1[, .(year, region, starts, src)],
                      r2[, .(year, region, starts, src)]), fill = TRUE)
chk("\nRegion labels:", paste(sort(unique(reg$region)), collapse = " | "))
ovr <- dcast(reg[year %in% 2020:2022], year + region ~ src, value.var = "starts")
chk("\nRegional overlap check 2020/21-2022/23 (since vs curgeo), first rows:")
chk(paste(capture.output(print(head(ovr[order(year, region)], 10))), collapse = "\n"))

reg[, prio := fifelse(src == "curgeo", 2L, 1L)]
setorder(reg, year, region, prio)
regu <- reg[!region %in% c("", "Outside of England and unknown"),
            .SD[.N], by = .(year, region)][year >= 2014 & year <= 2024]
fwrite(regu[order(region, year)], file.path(OUT, "tab_region_counts.csv"))

h223 <- fread(file.path(DATA, "hist/data/apps-geography-population-202223-q4.csv"), colClasses = "character")
h223 <- h223[geographic_level == "Regional" & apps_level == "Total" & age_group == "Total"]
h223[, `:=`(year = yr(time_period), region = gsub('"', "", region_name), starts = num(starts))]
rc <- grep("starts_rate", names(h223), value = TRUE)
if (length(rc) == 1) h223[, rate := num(get(rc))] else h223[, rate := NA_real_]
rates_full <- rbindlist(list(
  h223[year <= 2019, .(year, region, rate, src = "rel2223")],
  r2[, .(year, region, rate, src = "curgeo")]), fill = TRUE)
fwrite(rates_full[order(region, year)], file.path(OUT, "tab_region_rates_per100k.csv"))

# ---- E) Five manuscript figures (drives test_figures; RED:
#         'all(file.exists(figs)) is not TRUE') ----
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442")
thm <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        legend.position = "bottom")

f1d <- al[age_summary == "Total" & apps_level == "Total" & year >= 2009 & year <= 2024]
f1 <- ggplot(f1d, aes(year, starts / 1000)) +
  geom_vline(xintercept = c(2017, 2020), linetype = "dashed", colour = "grey55") +
  geom_line(colour = okabe[1], linewidth = 0.9) + geom_point(colour = okabe[1], size = 1.6) +
  annotate("text", x = 2017.1, y = max(f1d$starts / 1000) * 0.99, label = "Levy (April 2017)", hjust = 0, size = 3, colour = "grey30") +
  annotate("text", x = 2020.1, y = max(f1d$starts / 1000) * 0.90, label = "COVID-19 (March 2020)", hjust = 0, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq(2009, 2024, 2), labels = ay(seq(2009, 2024, 2))) +
  labs(title = "Apprenticeship starts in England, 2009/10 to 2024/25",
       subtitle = "Thousands. Dashed lines mark the apprenticeship levy and the COVID-19 shock.",
       x = NULL, y = "Starts (thousands)",
       caption = "Source: DfE, Apprenticeships historical summary (OGL v3.0).") + thm
ggsave(file.path(OUT, "fig1_total_starts.png"), f1, width = 7.2, height = 4.2, dpi = 300)

compA <- dec[apps_level == "Total" & age_summary != "Total"]
compA[, share := 100 * starts / sum(starts), by = year]
fwrite(compA, file.path(OUT, "comp_age.csv"))
f2 <- ggplot(compA, aes(year, share, fill = factor(age_summary, c("Under 19", "19-24", "25+")))) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = okabe[c(1, 3, 2)], name = NULL) +
  scale_x_continuous(breaks = 2014:2024, labels = ay(2014:2024)) +
  labs(title = "Age composition of apprenticeship starts, 2014/15 to 2024/25",
       subtitle = "Per cent of all starts.", x = NULL, y = "Share of starts (%)",
       caption = "Source: DfE (OGL v3.0).") + thm +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "fig2_age_composition.png"), f2, width = 7.2, height = 4.4, dpi = 300)

compL <- dec[age_summary == "Total" & apps_level != "Total"]
compL[, share := 100 * starts / sum(starts), by = year]
fwrite(compL, file.path(OUT, "comp_level.csv"))
lv_order <- c("Intermediate Apprenticeship", "Advanced Apprenticeship", "Higher Apprenticeship")
f3 <- ggplot(compL[apps_level %in% lv_order],
             aes(year, starts / 1000, colour = factor(apps_level, lv_order))) +
  geom_vline(xintercept = c(2017, 2020), linetype = "dashed", colour = "grey75") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_manual(values = okabe[c(4, 1, 3)], name = NULL,
                      labels = c("Intermediate (level 2)", "Advanced (level 3)", "Higher (levels 4+)")) +
  scale_x_continuous(breaks = 2014:2024, labels = ay(2014:2024)) +
  labs(title = "Apprenticeship starts by level, 2014/15 to 2024/25",
       subtitle = "Thousands. Dashed lines: levy (2017) and COVID-19 (2020).",
       x = NULL, y = "Starts (thousands)", caption = "Source: DfE (OGL v3.0).") + thm +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "fig3_level.png"), f3, width = 7.2, height = 4.4, dpi = 300)

f4 <- ggplot(g, aes(year, share, colour = factor(imd, imd_levels))) +
  geom_vline(xintercept = c(2017, 2020), linetype = "dashed", colour = "grey75") +
  geom_line(linewidth = 0.85) + geom_point(size = 1.4) +
  scale_colour_manual(values = okabe[c(4, 2, 7, 6, 1)], name = "IMD quintile") +
  scale_x_continuous(breaks = 2014:2024, labels = ay(2014:2024)) +
  labs(title = "Apprenticeship starts by learner home deprivation, 2014/15 to 2024/25",
       subtitle = "Share of starts with known IMD quintile (%).", x = NULL, y = "Share (%)",
       caption = "Source: DfE (OGL v3.0). Quintiles of the Index of Multiple Deprivation.") +
  thm + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "fig4_deprivation.png"), f4, width = 7.2, height = 4.6, dpi = 300)

sl <- regu[year %in% c(2014, 2024)]
sl[, yl := ay(year)]
f5 <- ggplot(sl, aes(yl, starts / 1000, group = region)) +
  geom_line(colour = "grey60") +
  geom_point(colour = okabe[1], size = 2) +
  geom_text(data = sl[year == 2024], aes(label = region), hjust = -0.06, size = 2.9) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.55))) +
  labs(title = "Apprenticeship starts by English region: 2014/15 versus 2024/25",
       subtitle = "Thousands of starts.", x = NULL, y = "Starts (thousands)",
       caption = "Source: DfE (OGL v3.0).") + thm
ggsave(file.path(OUT, "fig5_region_slope.png"), f5, width = 7.2, height = 5.2, dpi = 300)

close(log_chk)
cat("\nALL PARTS DONE. Outputs in", OUT, "\n")
