# =====================================================================
# Who still starts an apprenticeship?  England, 2014/15 to 2024/25
# Revision 1 analysis (October 2026), responding to the RPCE reviewer.
# Run 04a_build_population_denominators.R first. Works with Rscript or
# Source in RStudio (set RPCE_BASE or open the file from the project).
#
# New in this script (all outputs go to analysis/out_r1/):
#   R1-A  Deprivation ratio (Q1:Q5) and concentration index WITHIN LEVEL
#         and within AGE x LEVEL cells, 2014/15 to 2024/25, plus a
#         fixed-composition counterfactual.
#   R1-B  Starts per 1,000 population by IMD quintile (denominators:
#         ONS mid-2020 LSOA estimates mapped to IMD 2019 quintiles;
#         sensitivity with mid-2015 IMD File 7 denominators).
#   R1-C  Segmented regression extended back to 2002/03, with
#         Newey-West (HAC) standard errors and 95% confidence intervals;
#         sensitivity models (two-knot model, 2014/15 window, excluding
#         2020/21, HAC lag 1).
#   R1-D  Regional rates per 100,000 population aged 16-64 for the full
#         decade using ONS/Nomis mid-year estimates; dispersion measures.
#   R1-E  Composition of starts by sex, ethnicity and learning
#         difficulty/disability (LLDD), 2014/15 to 2024/25.
#   R1-F  Achievement rates by deprivation, age and level, 2024/25.
#   R1-G  Timing evidence: monthly starts around May 2017 and the share
#         of starts on standards (vs frameworks) by year.
#   R1-H  Journal-style figures (black, white and grey; Arial).
#
# Data: DfE Explore Education Statistics underlying CSVs (OGL v3.0);
#       ONS mid-year population estimates (Nomis API and ONS LSOA file);
#       MHCLG English Indices of Deprivation 2019, File 7.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sandwich)
  library(lmtest)
})
options(bitmapType = "quartz")   # macOS: lets png() use the system Arial font

BASE <- Sys.getenv("RPCE_BASE", unset = "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England")
if (!dir.exists(BASE) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
DATA <- file.path(BASE, "data")
EXT  <- file.path(DATA, "ext")
OUT  <- file.path(BASE, "analysis", "out_r1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(x))))
# DfE symbols: 'low' = fewer than 5 (treated as 0 here), 'z' = not
# applicable (programme did not exist; treated as 0 for series use)
num0 <- function(x) { v <- num(x); v[is.na(v) & trimws(x) %in% c("low", "z", "~", "c", "x")] <- 0; v }
yr  <- function(tp) as.integer(substr(gsub('"', "", as.character(tp)), 1, 4))
ay  <- function(y) sprintf("%d/%02d", y, (y + 1) %% 100)
logf <- file(file.path(OUT, "checks_r1.txt"), "w")
chk <- function(...) { cat(..., "\n"); cat(..., "\n", file = logf) }
rnd <- function(dt, d = 2) dt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, d) else x)]

imd_levels <- c("Q1", "Q2", "Q3", "Q4", "Q5")
std_imd <- function(x) {
  x <- trimws(gsub('"', "", x))
  out <- rep(NA_character_, length(x))
  out[grepl("^(1|One)",   x, ignore.case = TRUE)] <- "Q1"
  out[grepl("^(2|Two)",   x, ignore.case = TRUE)] <- "Q2"
  out[grepl("^(3|Three)", x, ignore.case = TRUE)] <- "Q3"
  out[grepl("^(4|Four)",  x, ignore.case = TRUE)] <- "Q4"
  out[grepl("^(5|Five)",  x, ignore.case = TRUE)] <- "Q5"
  out[is.na(out)] <- x[is.na(out)]
  out
}
std_level <- function(x) {
  x <- trimws(x)
  out <- rep(NA_character_, length(x))
  out[x %in% c("Intermediate Apprenticeship", "Intermediate", "Level 2", "2")] <- "Intermediate"
  out[x %in% c("Advanced Apprenticeship", "Advanced", "Level 3", "3")] <- "Advanced"
  out[x %in% c("Higher Apprenticeship", "Higher", "Level 4", "Level 5", "Level 6", "Level 7", "4", "5", "6", "7")] <- "Higher"
  out[x == "Total"] <- "Total"
  out[is.na(out)] <- "Unknown"
  out
}
conc_index <- function(q) 2 * sum(q * seq(0.1, 0.9, 0.2)) / sum(q) - 1   # q = starts in Q1..Q5

# =====================================================================
# 0. National historical series 2002/03 to 2024/25 (age x level)
# =====================================================================
hs <- fread(file.path(DATA, "cur/data/app-historical-summary-to-2425.csv"), colClasses = "character")
hs <- hs[ssa_tier_1 == "Total"]
hs[, `:=`(year = yr(time_period), starts = num0(starts), level = std_level(apps_level), age = age_summary)]
nat <- hs[, .(starts = sum(starts)), by = .(year, age, level)]
chk("Historical series: years", min(nat$year), "to", max(nat$year))
ser <- function(a, l) nat[age == a & level == l][order(year), .(year, starts)]
series_list <- list(
  "Total"        = ser("Total", "Total"),
  "Under 19"     = ser("Under 19", "Total"),
  "19-24"        = ser("19-24", "Total"),
  "25+"          = ser("25+", "Total"),
  "Intermediate" = ser("Total", "Intermediate"),
  "Advanced"     = ser("Total", "Advanced"),
  "Higher"       = ser("Total", "Higher"))
long_nat <- rbindlist(lapply(names(series_list), function(n) cbind(series = n, series_list[[n]])))
long_nat[, year_label := ay(year)]
fwrite(dcast(long_nat, year + year_label ~ series, value.var = "starts"), file.path(OUT, "tab_national_series_2002_2024.csv"))

# =====================================================================
# R1-C. Segmented regressions with HAC standard errors
# =====================================================================
# Two specifications.
# (i) seg_ls: interrupted time series with a LEVEL change and a SLOPE change
#     at each interruption (Wagner et al. 2002): for an interruption in year
#     k, step = 1[year >= k] and ramp = max(year - k, 0). This is the main
#     specification in the revision: with the series extended to 2002/03
#     the 2010/11 expansion and the 2017/18 levy year are step changes.
# (ii) seg_hac: the continuous piecewise-linear (slope-change only) model of
#     the submitted version, kept as a sensitivity analysis; the hinge for a
#     segment opening in year k is pmax(year - (k - 1), 0).
seg_ls <- function(df, label, knots, window, lag = 2, drop_year = NULL) {
  df <- as.data.table(df)[year >= window[1] & year <= window[2]][order(year)]
  if (!is.null(drop_year)) df <- df[!year %in% drop_year]
  X <- data.frame(y = df$starts, t = df$year - window[1])
  for (i in seq_along(knots)) { X[[paste0("L", i)]] <- as.integer(df$year >= knots[i]); X[[paste0("S", i)]] <- pmax(df$year - knots[i], 0) }
  m <- lm(y ~ ., data = X)
  V_hac <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  b <- coef(m); se_h <- sqrt(diag(V_hac)); se_o <- sqrt(diag(vcov(m)))
  dfres <- df.residual(m); tq <- qt(0.975, dfres)
  dw <- tryCatch(unname(dwtest(m)$statistic), error = function(e) NA_real_)
  nk <- length(knots); nb <- length(b)
  slope_idx <- c(2, 2 + 2 * seq_len(nk))            # t, S1, S2, ...
  level_idx <- 1 + 2 * seq_len(nk)                  # L1, L2, ...
  out <- data.table(series = label, n_obs = nrow(df), df_resid = dfres, dw = round(dw, 2), r2 = round(summary(m)$r.squared, 3))
  add <- function(nm, est, se) { out[[nm]] <<- est; out[[paste0(nm, "_se_hac")]] <<- se; out[[paste0(nm, "_lo95")]] <<- est - tq * se; out[[paste0(nm, "_hi95")]] <<- est + tq * se }
  for (i in seq_len(nk)) { add(paste0("level_change", i), b[level_idx[i]], se_h[level_idx[i]]); add(paste0("slope_change", i), b[slope_idx[i + 1]], se_h[slope_idx[i + 1]]) }
  out[["level_change_se_ols_levy"]] <- if (nk >= 2) se_o[level_idx[2]] else se_o[level_idx[1]]
  # implied segment slopes: cumulative sums of t, S1, S2, ...
  for (j in seq_len(nk + 1)) {
    Lm <- rep(0, nb); Lm[slope_idx[1:j]] <- 1
    add(paste0("slope_seg", j), sum(Lm * b), sqrt(as.numeric(t(Lm) %*% V_hac %*% Lm)))
  }
  attr(out, "fitted") <- data.table(year = df$year, fitted = fitted(m), observed = df$starts)
  out
}
run_ls <- function(knots, window, model_name, drop_year = NULL, lag = 2) {
  res <- rbindlist(lapply(names(series_list), function(n) seg_ls(series_list[[n]], n, knots, window, lag = lag, drop_year = drop_year)), fill = TRUE)
  res[, model := model_name]
  res
}
seg_hac <- function(df, label, knots, window, lag = 2, drop_year = NULL) {
  df <- as.data.table(df)[year >= window[1] & year <= window[2]][order(year)]
  if (!is.null(drop_year)) df <- df[!year %in% drop_year]
  X <- data.frame(y = df$starts, t = df$year - window[1])
  for (i in seq_along(knots)) X[[paste0("h", i)]] <- pmax(df$year - (knots[i] - 1), 0)
  m <- lm(y ~ ., data = X)
  V_hac <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  b <- coef(m); se_h <- sqrt(diag(V_hac)); se_o <- sqrt(diag(vcov(m)))
  dfres <- df.residual(m); tq <- qt(0.975, dfres)
  dw <- tryCatch(unname(dwtest(m)$statistic), error = function(e) NA_real_)
  nb <- length(b)
  # implied segment slopes = cumulative sums of the slope coefficients
  L <- matrix(0, nrow = nb - 1, ncol = nb); for (j in 1:(nb - 1)) L[j, 2:(j + 1)] <- 1
  seg_est <- as.vector(L %*% b); seg_se <- sqrt(diag(L %*% V_hac %*% t(L)))
  out <- data.table(series = label, n_obs = nrow(df), df_resid = dfres, dw = round(dw, 2), r2 = round(summary(m)$r.squared, 3))
  coef_names <- c("slope_base", paste0("change_knot", seq_along(knots)))
  for (j in 2:nb) {
    nm <- coef_names[j - 1]
    out[[nm]] <- b[j]; out[[paste0(nm, "_se_hac")]] <- se_h[j]; out[[paste0(nm, "_se_ols")]] <- se_o[j]
    out[[paste0(nm, "_lo95")]] <- b[j] - tq * se_h[j]; out[[paste0(nm, "_hi95")]] <- b[j] + tq * se_h[j]
  }
  for (j in seq_along(seg_est)) {
    nm <- paste0("slope_seg", j)
    out[[nm]] <- seg_est[j]; out[[paste0(nm, "_se_hac")]] <- seg_se[j]
    out[[paste0(nm, "_lo95")]] <- seg_est[j] - tq * seg_se[j]; out[[paste0(nm, "_hi95")]] <- seg_est[j] + tq * seg_se[j]
  }
  attr(out, "fitted") <- data.table(year = df$year, fitted = fitted(m), observed = df$starts)
  out
}
run_models <- function(knots, window, model_name, drop_year = NULL, lag = 2) {
  res <- rbindlist(lapply(names(series_list), function(n)
    seg_hac(series_list[[n]], n, knots, window, lag = lag, drop_year = drop_year)), fill = TRUE)
  res[, model := model_name]
  res
}
# ---- Main model M: level + slope changes at 2010/11, 2017/18, 2020/21; 2002/03-2024/25; Newey-West lag 2 ----
mM  <- run_ls(c(2010, 2017, 2020), c(2002, 2024), "M: level+slope changes 2010/11, 2017/18, 2020/21; 2002/03-2024/25; NW lag 2 (main)")
mM1 <- run_ls(c(2017, 2020), c(2010, 2024), "M1: level+slope changes 2017/18, 2020/21; window 2010/11-2024/25")
mM2 <- run_ls(c(2010, 2017, 2020), c(2002, 2024), "M2: as M, excluding 2020/21", drop_year = 2020)
mM3 <- run_ls(c(2010, 2017, 2020), c(2002, 2024), "M3: as M, Newey-West lag 1", lag = 1)
mM4 <- run_ls(c(2010, 2017, 2020), c(2002, 2024), "M4: as M, Newey-West lag 3", lag = 3)
# ---- Sensitivity: continuous slope-only models (submitted specification) ----
mA <- run_models(c(2010, 2017, 2020), c(2002, 2024), "A: slope-only, knots 2010/11, 2017/18, 2020/21; 2002/03-2024/25")
mB <- run_models(c(2017, 2020), c(2002, 2024), "B: slope-only, knots 2017/18, 2020/21; 2002/03-2024/25")
mC <- run_models(c(2017, 2020), c(2014, 2024), "C: slope-only, knots 2017/18, 2020/21; 2014/15-2024/25 (submitted model, now with NW SEs)")
fwrite(rbindlist(list(mM, mM1, mM2, mM3, mM4), fill = TRUE), file.path(OUT, "tab_segmented_levelslope_all_models.csv"))
fwrite(rbindlist(list(mA, mB, mC), fill = TRUE), file.path(OUT, "tab_segmented_slopeonly_models.csv"))

# Main manuscript table (Model M): pre-levy slope (2010/11-2016/17 segment), level change at the levy (2017/18),
# slope change at the levy, level change at the pandemic (2020/21), post-pandemic slope; HAC 95% CIs
mainM <- mM[, .(series, n_obs, df_resid, dw, r2,
                slope_2002_2009 = slope_seg1, level_2010 = level_change1, level_2010_lo = level_change1_lo95, level_2010_hi = level_change1_hi95,
                pre_levy_slope = slope_seg2, pre_lo = slope_seg2_lo95, pre_hi = slope_seg2_hi95,
                level_levy = level_change2, level_levy_lo = level_change2_lo95, level_levy_hi = level_change2_hi95, level_levy_se_hac = level_change2_se_hac, level_levy_se_ols = level_change_se_ols_levy,
                slope_levy = slope_change2, slope_levy_lo = slope_change2_lo95, slope_levy_hi = slope_change2_hi95,
                post_levy_slope = slope_seg3, post_levy_lo = slope_seg3_lo95, post_levy_hi = slope_seg3_hi95,
                level_covid = level_change3, level_covid_lo = level_change3_lo95, level_covid_hi = level_change3_hi95,
                slope_covid = slope_change3, slope_covid_lo = slope_change3_lo95, slope_covid_hi = slope_change3_hi95,
                post_covid_slope = slope_seg4, post_covid_lo = slope_seg4_lo95, post_covid_hi = slope_seg4_hi95)]
fwrite(rnd(mainM, 1), file.path(OUT, "tab_segmented_main_modelM.csv"))
chk("\nModel M (main): level + slope changes, Newey-West lag 2:")
chk(paste(capture.output(print(rnd(mainM[, .(series, pre_levy_slope, pre_lo, pre_hi, level_levy, level_levy_lo, level_levy_hi, slope_levy, slope_levy_lo, slope_levy_hi, post_levy_slope)], 0))), collapse = "\n"))
chk(paste(capture.output(print(rnd(mainM[, .(series, level_covid, level_covid_lo, level_covid_hi, slope_covid, slope_covid_lo, slope_covid_hi, post_covid_slope, post_covid_lo, post_covid_hi, level_2010, slope_2002_2009, dw, r2, level_levy_se_ols, level_levy_se_hac)], 0))), collapse = "\n"))
pickM <- function(m, pre, lev, slp, post) m[, .(model, series, pre_levy_slope = get(pre), level_levy = get(lev), level_levy_lo = get(paste0(lev, "_lo95")), level_levy_hi = get(paste0(lev, "_hi95")),
                                                slope_levy = get(slp), slope_levy_lo = get(paste0(slp, "_lo95")), slope_levy_hi = get(paste0(slp, "_hi95")),
                                                post_covid_slope = get(post), post_covid_lo = get(paste0(post, "_lo95")), post_covid_hi = get(paste0(post, "_hi95")))]
sensM <- rbindlist(list(pickM(mM, "slope_seg2", "level_change2", "slope_change2", "slope_seg4"),
                        pickM(mM1, "slope_seg1", "level_change1", "slope_change1", "slope_seg3"),
                        pickM(mM2, "slope_seg2", "level_change2", "slope_change2", "slope_seg4"),
                        pickM(mM3, "slope_seg2", "level_change2", "slope_change2", "slope_seg4"),
                        pickM(mM4, "slope_seg2", "level_change2", "slope_change2", "slope_seg4")))
pick <- function(m, pre, lev, cov, post) m[, .(model, series, pre_levy_slope = get(pre), slope_levy = get(lev), slope_levy_lo = get(paste0(lev, "_lo95")), slope_levy_hi = get(paste0(lev, "_hi95")),
                                               post_covid_slope = get(post), post_covid_lo = get(paste0(post, "_lo95")), post_covid_hi = get(paste0(post, "_hi95")))]
sensS <- rbindlist(list(pick(mA, "slope_seg2", "change_knot2", "change_knot3", "slope_seg4"),
                        pick(mB, "slope_seg1", "change_knot1", "change_knot2", "slope_seg3"),
                        pick(mC, "slope_seg1", "change_knot1", "change_knot2", "slope_seg3")))
sens <- rbindlist(list(sensM, sensS), fill = TRUE)
fwrite(rnd(sens, 0), file.path(OUT, "tab_segmented_sensitivity.csv"))
chk("\nSensitivity across models (level change at levy where estimated; slope change at levy; post-COVID slope):")
chk(paste(capture.output(print(rnd(sens[, .(model = sub(":.*", "", model), series, pre_levy_slope, level_levy, level_levy_lo, level_levy_hi, slope_levy, slope_levy_lo, slope_levy_hi, post_covid_slope, post_covid_lo, post_covid_hi)], 0))), collapse = "\n"))
fitM_total <- attr(seg_ls(series_list[["Total"]], "Total", c(2010, 2017, 2020), c(2002, 2024)), "fitted")
fwrite(fitM_total, file.path(OUT, "fitted_modelM_total.csv"))
fitA_total <- fitM_total

# =====================================================================
# 1. Deprivation files spliced: level x age x IMD, 2014/15 to 2024/25
# =====================================================================
d1 <- fread(file.path(DATA, "hist1920/data/app-learner-deprivation-starts-1920-q4.csv"), colClasses = "character")
d1 <- d1[geographic_level == "National" & ssa_t1_desc == "Total" & age_youth_adult == "Total"]
d1[, `:=`(year = yr(time_period), starts = num0(starts), imd = std_imd(learner_home_depriv),
          level = std_level(apps_level), age = age_summary, src = "rel1920")]
# In the 2022/23 and 2025/26 releases the age columns are hierarchical: age_summary rows sit under
# their age_youth_adult parent ('Under 19' or '19+'), so select the parent-level Total row and the
# age_summary rows (each age_summary value occurs once).
age_rows <- function(dt) dt[(age_youth_adult == "Total" & age_summary == "Total") | age_summary != "Total"]
d2 <- fread(file.path(DATA, "hist/data/app-learner-deprivation-starts-202223-q4.csv"), colClasses = "character")
d2 <- age_rows(d2[geographic_level == "National" & ssa_t1_desc == "Total"])
d2[, `:=`(year = yr(time_period), starts = num0(starts), imd = std_imd(learner_home_depriv),
          level = std_level(apps_level_detailed), age = age_summary, src = "rel2223")]
d2 <- d2[, .(starts = sum(starts)), by = .(year, imd, level, age, src)]   # levels 4-7 -> Higher
d3 <- fread(file.path(DATA, "cur/data/apps_11_learner_deprivation_202526_10.csv"), colClasses = "character")
d3 <- age_rows(d3[geographic_level == "National" & ssa_tier_1 == "Total" & apps_level_detailed == "Total"])
d3[, `:=`(year = yr(time_period), starts = num0(starts), imd = std_imd(learner_home_depriv),
          level = std_level(apps_level), age = age_summary, src = "cur")]
d3 <- d3[year <= 2024]
stopifnot(!anyDuplicated(d3[, .(year, imd, level, age)]), !anyDuplicated(d2[, .(year, imd, level, age)]))
dep <- rbindlist(list(d1[, .(year, imd, level, age, starts, src)], d2[, .(year, imd, level, age, starts, src)], d3[, .(year, imd, level, age, starts, src)]))
dep <- dep[imd %in% c(imd_levels, "Total", "Unknown") & level != "Unknown"]
ovl <- dcast(dep[imd %in% imd_levels & age == "Total" & level != "Total"], year + level + imd ~ src, value.var = "starts", fun.aggregate = sum)
chk("\nDeprivation x level overlap: max |rel2223 - rel1920| 2017-2019 =", max(abs(ovl[year %in% 2017:2019, rel2223 - rel1920]), na.rm = TRUE),
    "; max |cur - rel2223| 2020-2022 =", max(abs(ovl[year %in% 2020:2022, cur - rel2223]), na.rm = TRUE))
dep[, prio := fifelse(src == "cur", 3L, fifelse(src == "rel2223", 2L, 1L))]
setorder(dep, year, imd, level, age, prio)
depu <- dep[, .SD[.N], by = .(year, imd, level, age)]
depk <- depu[imd %in% imd_levels]
chk_lv <- merge(depk[level != "Total", .(s = sum(starts)), by = .(year, age)], depk[level == "Total", .(tot = sum(starts)), by = .(year, age)], by = c("year", "age"))
chk("Level cells vs level total (known quintile): max abs diff =", max(abs(chk_lv$s - chk_lv$tot)))
unk <- depu[age == "Total" & level == "Total", .(unknown_pct = 100 * sum(starts[imd == "Unknown"]) / sum(starts[imd != "Total"])), by = year][order(year)]
chk("Share of starts with unknown IMD quintile by year (%):", paste(ay(unk$year), round(unk$unknown_pct, 2), collapse = "; "))

# ---- R1-A: ratios and concentration indices within groups ----
grp <- depk[, .(starts = sum(starts)), by = .(year, imd, level, age)]
grp_w <- dcast(grp, year + level + age ~ imd, value.var = "starts")
grp_w[, `:=`(ratio_Q1_Q5 = Q1 / Q5, known_starts = Q1 + Q2 + Q3 + Q4 + Q5, year_label = ay(year))]
grp_w[, conc_index := apply(.SD, 1, conc_index), .SDcols = imd_levels]
fwrite(grp_w[order(level, age, year)], file.path(OUT, "tab_dep_within_groups_all_years.csv"))
sel_years <- c(2014, 2016, 2017, 2019, 2024)
within_tab <- grp_w[year %in% sel_years & !(level == "Total" & age == "Total")]
within_tab[, group := fifelse(level == "Total", paste("Age", age), fifelse(age == "Total", paste("Level", level), paste0(age, " x ", level)))]
wt_ratio <- dcast(within_tab, group + level + age ~ year_label, value.var = "ratio_Q1_Q5")
wt_ci    <- dcast(within_tab, group + level + age ~ year_label, value.var = "conc_index")
wt_n     <- dcast(within_tab, group + level + age ~ year_label, value.var = "known_starts")
fwrite(wt_ratio, file.path(OUT, "tab_within_ratio_selected_years.csv"))
fwrite(wt_ci, file.path(OUT, "tab_within_ci_selected_years.csv"))
fwrite(wt_n, file.path(OUT, "tab_within_n_selected_years.csv"))
chk("\nR1-A Q1:Q5 ratio within groups:"); chk(paste(capture.output(print(rnd(wt_ratio, 2))), collapse = "\n"))
chk("\nR1-A concentration index within groups:"); chk(paste(capture.output(print(rnd(wt_ci, 3))), collapse = "\n"))
chk("\nR1-A known-quintile starts within groups:"); chk(paste(capture.output(print(wt_n)), collapse = "\n"))
# Counterfactual: hold the 2014/15 age x level mix fixed, apply each year's observed within-cell quintile shares
cells <- grp[level != "Total" & age != "Total"]
cells[, share_q := starts / sum(starts), by = .(year, level, age)]
w14 <- cells[year == 2014, .(w = sum(starts)), by = .(level, age)][, w := w / sum(w)]
cf <- merge(cells, w14, by = c("level", "age"))[, .(cf = sum(share_q * w)), by = .(year, imd)]
cf <- dcast(cf, year ~ imd, value.var = "cf")
cf[, `:=`(cf_ratio_fixed_2014_mix = Q1 / Q5, cf_ci_fixed_2014_mix = apply(.SD, 1, conc_index)), .SDcols = imd_levels]
cf_out <- merge(grp_w[level == "Total" & age == "Total", .(year, year_label, observed_ratio = ratio_Q1_Q5, observed_ci = conc_index)],
                cf[, .(year, cf_ratio_fixed_2014_mix, cf_ci_fixed_2014_mix)], by = "year")
fwrite(cf_out, file.path(OUT, "tab_dep_counterfactual_fixed_mix.csv"))
chk("\nCounterfactual (fixed 2014/15 age x level mix) vs observed:"); chk(paste(capture.output(print(rnd(cf_out, 3))), collapse = "\n"))

# =====================================================================
# R1-B. Starts per 1,000 population by IMD quintile
# =====================================================================
pop <- fread(file.path(EXT, "imd2019_quintile_population.csv"))
pop[, imd := paste0("Q", imd_quintile)]
setnames(pop, grep("^imd2015_Population aged 16-59", names(pop), value = TRUE), "pop_16_59_mid2015")
setnames(pop, grep("^imd2015_Total population", names(pop), value = TRUE), "pop_all_mid2015")
pop_long <- melt(pop[, .(imd, pop_16_18, pop_19_24, pop_16_24, pop_25_64, pop_16_64, pop_16_59_mid2015, pop_all_mid2015)],
                 id.vars = "imd", variable.name = "denom", value.name = "pop")
chk("\nPopulation shares by quintile (mid-2020 unless stated):")
chk(paste(capture.output(print(pop[, .(imd, share_16_18 = round(100 * pop_16_18 / sum(pop_16_18), 1), share_19_24 = round(100 * pop_19_24 / sum(pop_19_24), 1),
      share_16_24 = round(100 * pop_16_24 / sum(pop_16_24), 1), share_25_64 = round(100 * pop_25_64 / sum(pop_25_64), 1), share_16_64 = round(100 * pop_16_64 / sum(pop_16_64), 1),
      share_16_59_mid2015 = round(100 * pop_16_59_mid2015 / sum(pop_16_59_mid2015), 1))])), collapse = "\n"))
age_tot <- grp[level == "Total"]
u25 <- age_tot[age %in% c("Under 19", "19-24"), .(starts = sum(starts), level = "Total", age = "Under 25"), by = .(year, imd)]
age_tot2 <- rbindlist(list(age_tot, u25), use.names = TRUE)
rate_def <- data.table(age = c("Under 19", "19-24", "25+", "Total", "Under 25"), denom = c("pop_16_18", "pop_19_24", "pop_25_64", "pop_16_64", "pop_16_24"))
rates <- merge(merge(age_tot2, rate_def, by = "age"), pop_long, by = c("imd", "denom"))
rates[, `:=`(rate_per_1000 = 1000 * starts / pop, year_label = ay(year))]
rate_w <- dcast(rates, age + denom + year + year_label ~ imd, value.var = "rate_per_1000")
rate_w[, rate_ratio_Q1_Q5 := Q1 / Q5]
fwrite(rate_w[order(age, year)], file.path(OUT, "tab_rates_per_1000_by_quintile.csv"))
chk("\nR1-B starts per 1,000 population by quintile (all years):")
chk(paste(capture.output(print(rnd(rate_w[, .(age, year_label, Q1, Q2, Q3, Q4, Q5, ratio = rate_ratio_Q1_Q5)][order(age, year_label)], 2))), collapse = "\n"))
sens_rate <- merge(age_tot[age == "Total"], pop_long[denom == "pop_16_59_mid2015"], by = "imd")[, rate := 1000 * starts / pop]
sens_w <- dcast(sens_rate, year ~ imd, value.var = "rate")[, .(year, year_label = ay(year), ratio_mid2015_16_59 = Q1 / Q5)]
sens_w <- merge(sens_w, rate_w[age == "Total", .(year, ratio_mid2020_16_64 = rate_ratio_Q1_Q5)], by = "year")
fwrite(sens_w, file.path(OUT, "tab_rates_sensitivity_denominator_vintage.csv"))
chk("\nRate ratio sensitivity to denominator vintage (all ages):"); chk(paste(capture.output(print(rnd(sens_w, 3))), collapse = "\n"))

# =====================================================================
# R1-D. Regional rates per 100,000 population aged 16-64, 2014/15-2024/25
# =====================================================================
r1 <- fread(file.path(DATA, "hist/data/app-starts-since-202223-q4.csv"), colClasses = "character")
r1 <- r1[geographic_level == "Regional"]
r1[, `:=`(year = yr(time_period), starts = num0(starts), region = gsub('"', "", region_name), src = "since")]
r2 <- fread(file.path(DATA, "cur/data/app-geography-population-202526-q3.csv"), colClasses = "character")
r2 <- r2[geographic_level == "Regional" & apps_level == "Total" & age_summary == "Total"]
r2[, `:=`(year = yr(time_period), starts = num0(starts), region = gsub('"', "", region_name),
          dfe_rate = num(starts_rate_per_100000_population), dfe_pop = num(population_estimate), src = "curgeo")]
r2 <- r2[year <= 2024]
h223 <- fread(file.path(DATA, "hist/data/apps-geography-population-202223-q4.csv"), colClasses = "character")
h223 <- h223[geographic_level == "Regional" & apps_level == "Total" & age_group == "Total"]
h223[, `:=`(year = yr(time_period), region = gsub('"', "", region_name), starts = num0(starts),
            dfe_rate = num(starts_rate_per_100000_population), dfe_pop = num(population_estimate))]
reg <- rbindlist(list(r1[, .(year, region, starts, src)], r2[, .(year, region, starts, src)]), fill = TRUE)
reg[, prio := fifelse(src == "curgeo", 2L, 1L)]
setorder(reg, year, region, prio)
regu <- reg[!region %in% c("", "Outside of England and unknown"), .SD[.N], by = .(year, region)][year >= 2009 & year <= 2024]
nom <- fread(file.path(EXT, "nomis_region_pop_2009_2024.csv"))
nom[, region := fifelse(GEOGRAPHY_NAME == "East", "East of England", GEOGRAPHY_NAME)]
nom16 <- nom[C_AGE_NAME == "Aged 16 to 64", .(pop_year = as.integer(DATE_NAME), region, pop_16_64 = as.numeric(OBS_VALUE))]
dfe_rates <- rbindlist(list(h223[year <= 2019, .(year, region, dfe_rate, dfe_pop)], r2[, .(year, region, dfe_rate, dfe_pop)]))
chk_conv <- merge(dfe_rates, nom16[, .(year = pop_year, region, pop_start = pop_16_64)], by = c("year", "region"))
chk_conv <- merge(chk_conv, nom16[, .(year = pop_year - 1, region, pop_end = pop_16_64)], by = c("year", "region"))
chk("\nDfE denominator vs Nomis 16-64: mean abs % diff, start-year mid-year =", round(100 * mean(abs(chk_conv$dfe_pop - chk_conv$pop_start) / chk_conv$dfe_pop), 2),
    "%; end-year mid-year =", round(100 * mean(abs(chk_conv$dfe_pop - chk_conv$pop_end) / chk_conv$dfe_pop), 2), "%")
rr <- merge(regu, nom16[, .(year = pop_year, region, pop_16_64)], by = c("year", "region"))
rr[, rate_per_100k := 1e5 * starts / pop_16_64]
rr <- merge(rr, dfe_rates[, .(year, region, dfe_rate)], by = c("year", "region"), all.x = TRUE)
rr[, year_label := ay(year)]
chk("Own rate vs DfE published rate (2017/18+): correlation =", round(cor(rr$rate_per_100k, rr$dfe_rate, use = "complete.obs"), 4),
    "; mean abs diff =", round(mean(abs(rr$rate_per_100k - rr$dfe_rate), na.rm = TRUE), 1), "per 100,000")
fwrite(rr[order(region, year)], file.path(OUT, "tab_region_rates_full_decade.csv"))
disp <- rr[year >= 2014, .(max_min_ratio = max(rate_per_100k) / min(rate_per_100k), cv = sd(rate_per_100k) / mean(rate_per_100k),
                           england_rate = 1e5 * sum(starts) / sum(pop_16_64),
                           max_region = region[which.max(rate_per_100k)], min_region = region[which.min(rate_per_100k)]), by = year][order(year)]
disp[, year_label := ay(year)]
fwrite(disp, file.path(OUT, "tab_region_dispersion_full_decade.csv"))
chk("\nR1-D regional dispersion of rates per 100,000 aged 16-64:"); chk(paste(capture.output(print(rnd(disp, 3))), collapse = "\n"))
rr_sel <- dcast(rr[year %in% c(2014, 2016, 2019, 2024)], region ~ year_label, value.var = "rate_per_100k")
rr_sel <- merge(rr_sel, rr[year %in% c(2014, 2024), .(change_pct = 100 * (rate_per_100k[year == 2024] / rate_per_100k[year == 2014] - 1)), by = region], by = "region")
fwrite(rr_sel[order(change_pct)], file.path(OUT, "tab_region_rates_selected.csv"))
chk(paste(capture.output(print(rnd(rr_sel[order(change_pct)], 1))), collapse = "\n"))

# =====================================================================
# R1-E. Sex, ethnicity and LLDD composition, 2014/15 to 2024/25
# =====================================================================
l1 <- fread(file.path(DATA, "hist1920/data/app-learner-detailed-1920-q4-20201125.csv"), colClasses = "character")
l1 <- l1[age_youth_adult == "Total" & ethnicty_group == "Total" & ethnicity_detailed == "Total"]
l1[, `:=`(year = yr(time_period), starts = num0(starts), level = std_level(level), age = age_summary, src = "rel1920", sex = gender,
          eth = fifelse(bame == "BAME", "Ethnic minorities (excluding white minorities)", fifelse(bame == "Not BAME", "White", fifelse(bame == "Total", "Total", "Unknown"))),
          lldd = fifelse(lldd == "LLDD - yes", "With LLDD", fifelse(lldd == "LLDD - no", "No LLDD", fifelse(lldd == "Total", "Total", "Unknown"))))]
l2 <- fread(file.path(DATA, "hist/data/app-learner-detailed-202223-q4.csv"), colClasses = "character")
l2 <- age_rows(l2[ethnicity_major == "Total" & ethnicity_minor == "Total"])
l2[, `:=`(year = yr(time_period), starts = num0(starts), level = std_level(level), age = age_summary, src = "rel2223",
          eth = fifelse(minority_ethnic == "Ethnic minorities (excluding white minorities)", "Ethnic minorities (excluding white minorities)",
                        fifelse(minority_ethnic == "White", "White", fifelse(minority_ethnic == "Total", "Total", "Unknown"))),
          lldd = fifelse(lldd == "LLDD - yes", "With LLDD", fifelse(lldd == "LLDD - no", "No LLDD", fifelse(lldd == "Total", "Total", "Unknown"))))]
l3 <- fread(file.path(DATA, "cur/data/apps_09_learner_detailed_202526_10.csv"), colClasses = "character")
l3 <- age_rows(l3[ethnicity_major == "Total" & ethnicity_minor == "Total"])
l3[, `:=`(year = yr(time_period), starts = num0(starts), level = std_level(apps_level), age = age_summary, src = "cur",
          eth = fifelse(minority_ethnic == "Ethnic minorities (excluding white minorities)", "Ethnic minorities (excluding white minorities)",
                        fifelse(minority_ethnic == "White", "White", fifelse(minority_ethnic == "Total", "Total", "Unknown"))),
          lldd = fifelse(learning_difficulty_disability == "With learning difficulty / disability", "With LLDD",
                         fifelse(learning_difficulty_disability == "No learning difficulty / disability", "No LLDD",
                                 fifelse(learning_difficulty_disability == "Total", "Total", "Unknown"))))]
l3 <- l3[year <= 2024]
lc <- rbindlist(list(l1[, .(year, age, level, sex, eth, lldd, starts, src)], l2[, .(year, age, level, sex, eth, lldd, starts, src)], l3[, .(year, age, level, sex, eth, lldd, starts, src)]))
lc <- lc[level != "Unknown"]
lc[, prio := fifelse(src == "cur", 3L, fifelse(src == "rel2223", 2L, 1L))]
ovl_sex  <- dcast(lc[age == "Total" & level == "Total" & eth == "Total" & lldd == "Total" & sex %in% c("Female", "Male")], year + sex ~ src, value.var = "starts", fun.aggregate = sum)
ovl_eth  <- dcast(lc[age == "Total" & level == "Total" & sex == "Total" & lldd == "Total" & eth != "Total"], year + eth ~ src, value.var = "starts", fun.aggregate = sum)
ovl_lldd <- dcast(lc[age == "Total" & level == "Total" & sex == "Total" & eth == "Total" & lldd != "Total"], year + lldd ~ src, value.var = "starts", fun.aggregate = sum)
chk("\nR1-E overlap checks (max abs diff between releases in overlapping years):")
for (nm in c("ovl_sex", "ovl_eth", "ovl_lldd")) {
  o <- get(nm)
  if (all(c("rel1920", "rel2223") %in% names(o))) chk(" ", nm, "rel1920 vs rel2223 (2017-2019):", max(abs(o[year %in% 2017:2019, rel2223 - rel1920]), na.rm = TRUE))
  if (all(c("rel2223", "cur") %in% names(o))) chk(" ", nm, "rel2223 vs cur (2018-2022):", max(abs(o[year %in% 2018:2022, cur - rel2223]), na.rm = TRUE))
}
setorder(lc, year, age, level, sex, eth, lldd, prio)
lcu <- lc[, .SD[.N], by = .(year, age, level, sex, eth, lldd)]
sex_t <- lcu[age == "Total" & level == "Total" & eth == "Total" & lldd == "Total" & sex %in% c("Female", "Male")][, share := 100 * starts / sum(starts), by = year]
eth_t <- lcu[age == "Total" & level == "Total" & sex == "Total" & lldd == "Total" & eth %in% c("White", "Ethnic minorities (excluding white minorities)")][, share := 100 * starts / sum(starts), by = year]
eth_unk <- lcu[age == "Total" & level == "Total" & sex == "Total" & lldd == "Total" & eth != "Total"][, .(eth_unknown_pct = 100 * sum(starts[eth == "Unknown"]) / sum(starts)), by = year]
lldd_t <- lcu[age == "Total" & level == "Total" & sex == "Total" & eth == "Total" & lldd %in% c("With LLDD", "No LLDD")][, share := 100 * starts / sum(starts), by = year]
lldd_unk <- lcu[age == "Total" & level == "Total" & sex == "Total" & eth == "Total" & lldd != "Total"][, .(lldd_unknown_pct = 100 * sum(starts[lldd == "Unknown"]) / sum(starts)), by = year]
comp <- Reduce(function(a, b) merge(a, b, by = "year"), list(
  sex_t[sex == "Female", .(year, female_pct = share, female_starts = starts)],
  eth_t[eth != "White", .(year, ethnic_minority_pct = share, ethnic_minority_starts = starts)], eth_unk,
  lldd_t[lldd == "With LLDD", .(year, lldd_pct = share, lldd_starts = starts)], lldd_unk))
comp[, year_label := ay(year)]
fwrite(comp, file.path(OUT, "tab_sex_ethnicity_lldd_by_year.csv"))
chk("\nR1-E composition by year:"); chk(paste(capture.output(print(rnd(comp[, .(year_label, female_pct, ethnic_minority_pct, eth_unknown_pct, lldd_pct, lldd_unknown_pct)], 1))), collapse = "\n"))
bl <- function(var, cats, keep) {
  x <- lcu[year %in% c(2014, 2016, 2019, 2024) & get(var) %in% cats]
  x <- x[if (var == "sex") eth == "Total" & lldd == "Total" else if (var == "eth") sex == "Total" & lldd == "Total" else sex == "Total" & eth == "Total"]
  x <- x[(level != "Total" & age == "Total") | (level == "Total" & age != "Total") | (level == "Total" & age == "Total")]
  x[, share := 100 * starts / sum(starts), by = .(year, level, age)]
  x[get(var) == keep, .(year_label = ay(year), level, age, share = round(share, 1), starts)]
}
sub_tab <- rbindlist(list(cbind(characteristic = "Female", bl("sex", c("Female", "Male"), "Female")),
                          cbind(characteristic = "Ethnic minorities (excl. white minorities)", bl("eth", c("White", "Ethnic minorities (excluding white minorities)"), "Ethnic minorities (excluding white minorities)")),
                          cbind(characteristic = "With LLDD", bl("lldd", c("With LLDD", "No LLDD"), "With LLDD"))))
fwrite(sub_tab, file.path(OUT, "tab_sex_ethnicity_lldd_by_level_age.csv"))
chk("\nR1-E shares by level and age (selected years):"); chk(paste(capture.output(print(dcast(sub_tab, characteristic + level + age ~ year_label, value.var = "share"))), collapse = "\n"))

# =====================================================================
# R1-F. Achievement rates 2024/25 by deprivation, age and level (NARTS)
# =====================================================================
na3 <- fread(file.path(DATA, "cur/data/na03_app_narts_deprivation_202526_14.csv"), colClasses = "character")
na3 <- na3[time_period == "202425" & std_fwk_flag == "Total" & stem == "Total" & ssa_tier_1 == "Total" & detailed_level == "Total" & age_youth_adult == "Total"]
na3[, `:=`(imd = std_imd(learner_home_depriv), ach = num(achievement_rate), ret = num(retention_rate), leavers = num(leavers))]
ach <- na3[imd %in% c(imd_levels, "Total"), .(level, age_group, imd, leavers, achievement_rate = ach, retention_rate = ret)]
fwrite(ach, file.path(OUT, "tab_achievement_rates_2024_25.csv"))
chk("\nR1-F achievement rates 2024/25 by IMD:"); chk(paste(capture.output(print(dcast(ach, level + age_group ~ imd, value.var = "achievement_rate"))), collapse = "\n"))
d3a <- fread(file.path(DATA, "cur/data/apps_11_learner_deprivation_202526_10.csv"), colClasses = "character")
d3a <- d3a[geographic_level == "National" & ssa_tier_1 == "Total" & age_youth_adult == "Total" & age_summary == "Total" & apps_level_detailed == "Total" & apps_level == "Total" & time_period == "202425"]
d3a[, `:=`(imd = std_imd(learner_home_depriv), achievements = num0(achievements), starts = num0(starts))]
d3a <- merge(d3a[imd %in% imd_levels, .(imd, starts, achievements)], pop[, .(imd, pop_16_64)], by = "imd")
d3a[, `:=`(starts_per_1000 = 1000 * starts / pop_16_64, ach_per_1000 = 1000 * achievements / pop_16_64)]
fwrite(d3a, file.path(OUT, "tab_starts_achievements_per_1000_2024_25.csv"))
chk("\nStarts and achievements per 1,000 aged 16-64 by quintile, 2024/25:"); chk(paste(capture.output(print(rnd(d3a, 2))), collapse = "\n"))

# =====================================================================
# R1-G. Timing evidence: monthly starts and standards share
# =====================================================================
u10 <- fread(file.path(DATA, "cur/supporting-files/apps_u10_underlying_monthly_starts_202526_10.csv"), colClasses = "character")
u10 <- u10[grepl("Full Year", academic_year)]
u10[, `:=`(year = as.integer(substr(academic_year, 1, 4)), starts = num0(starts), m = as.integer(start_month_number))]
std_share <- u10[, .(starts = sum(starts), standards = sum(starts[std_fwk_flag == "Standard"]), levy_funded = sum(starts[funding_type == "Supported by ASA levy funds"])), by = year][order(year)]
std_share[, `:=`(standards_pct = 100 * standards / starts, levy_funded_pct = 100 * levy_funded / starts, year_label = ay(year))]
fwrite(std_share, file.path(OUT, "tab_standards_levy_share_by_year.csv"))
chk("\nR1-G standards and levy-funded share of starts by year:"); chk(paste(capture.output(print(rnd(std_share[, .(year_label, starts, standards_pct, levy_funded_pct)], 1))), collapse = "\n"))
mon <- u10[year %in% 2015:2017, .(starts = sum(starts), framework = sum(starts[std_fwk_flag == "Framework"]), standard = sum(starts[std_fwk_flag == "Standard"]),
                                  under19 = sum(starts[age_summary == "Under 19"]), a1924 = sum(starts[age_summary == "19-24"]), a25 = sum(starts[age_summary == "25+"]),
                                  intermediate = sum(starts[apps_level == "Intermediate Apprenticeship"])), by = .(year, m, start_month)][order(year, m)]
fwrite(mon, file.path(OUT, "tab_monthly_starts_2015_2017.csv"))
mj <- merge(mon[m >= 10, .(may_jul = sum(starts), may_jul_framework = sum(framework), may_jul_under19 = sum(under19), may_jul_25plus = sum(a25)), by = year],
            mon[m <= 9, .(aug_apr = sum(starts), aug_apr_framework = sum(framework)), by = year], by = "year")[, year_label := ay(year)]
fwrite(mj, file.path(OUT, "tab_may_july_vs_aug_april.csv"))
chk("\nMay-July vs August-April starts, 2015/16 to 2017/18:"); chk(paste(capture.output(print(mj)), collapse = "\n"))
chk("Monthly 2016/17:"); chk(paste(capture.output(print(mon[year == 2016, .(start_month, starts, framework, standard, under19, a25, intermediate)])), collapse = "\n"))

# =====================================================================
# R1-H. Figures: journal line art (black/white/grey, Arial)
# =====================================================================
thm <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "black"), axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9),
        legend.position = "bottom", legend.title = element_blank(), legend.key.width = unit(1.1, "cm"), legend.text = element_text(size = 8),
        strip.background = element_blank(), strip.text = element_text(size = 8.5, hjust = 0), plot.margin = margin(4, 8, 4, 4))
f1d <- series_list[["Total"]]
f1 <- ggplot(f1d, aes(year, starts / 1000)) +
  geom_vline(xintercept = c(2016.5, 2019.5), linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_line(data = fitA_total, aes(year, fitted / 1000), colour = "grey55", linewidth = 0.5, linetype = "longdash") +
  geom_line(colour = "black", linewidth = 0.6) + geom_point(colour = "black", size = 1.3) +
  annotate("text", x = 2016.4, y = 560, label = "Levy\n(April 2017)", hjust = 1, size = 2.6, family = "Arial", colour = "grey25") +
  annotate("text", x = 2019.6, y = 560, label = "COVID-19\n(March 2020)", hjust = 0, size = 2.6, family = "Arial", colour = "grey25") +
  scale_x_continuous(breaks = seq(2002, 2024, 2), labels = ay(seq(2002, 2024, 2)), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 580), breaks = seq(0, 500, 100), expand = expansion(mult = 0)) +
  labs(x = NULL, y = "Starts (thousands)") + thm + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "fig1_total_starts_2002_2024.png"), f1, width = 6.5, height = 3.6, dpi = 600, bg = "white")

f2d <- long_nat[series %in% c("Intermediate", "Advanced", "Higher") & year >= 2014]
f2d[, series := factor(series, c("Intermediate", "Advanced", "Higher"), c("Intermediate (level 2)", "Advanced (level 3)", "Higher (levels 4 and above)"))]
f2 <- ggplot(f2d, aes(year, starts / 1000, linetype = series, shape = series)) +
  geom_vline(xintercept = c(2016.5, 2019.5), linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_line(colour = "black", linewidth = 0.55) + geom_point(colour = "black", size = 1.5) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash")) + scale_shape_manual(values = c(16, 15, 17)) +
  scale_x_continuous(breaks = 2014:2024, labels = ay(2014:2024)) + scale_y_continuous(limits = c(0, 320), expand = expansion(mult = 0)) +
  labs(x = NULL, y = "Starts (thousands)") + thm + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "fig2_level.png"), f2, width = 6.5, height = 3.6, dpi = 600, bg = "white")

sh <- grp_w[level == "Total" & age == "Total"]
sh_long <- melt(sh[, .(year, Q1, Q2, Q3, Q4, Q5)], id.vars = "year", variable.name = "imd", value.name = "starts")
sh_long[, share := 100 * starts / sum(starts), by = year][, panel := "(a) Share of starts with a known IMD quintile (%)"]
rt <- melt(rate_w[age == "Under 25", .(year, Q1, Q2, Q3, Q4, Q5)], id.vars = "year", variable.name = "imd", value.name = "value")
rt[, panel := "(b) Under-25 starts per 1,000 residents aged 16 to 24"]
f3d <- rbindlist(list(sh_long[, .(year, imd, value = share, panel)], rt), use.names = TRUE)
f3d[, imd := factor(imd, imd_levels, c("Q1 (most deprived)", "Q2", "Q3", "Q4", "Q5 (least deprived)"))]
f3 <- ggplot(f3d, aes(year, value, linetype = imd, shape = imd, colour = imd)) +
  geom_vline(xintercept = c(2016.5, 2019.5), linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.4) + facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash", "solid")) + scale_shape_manual(values = c(16, 15, 17, 18, 1)) +
  scale_colour_manual(values = c("black", "grey30", "grey45", "grey60", "black")) +
  scale_x_continuous(breaks = seq(2014, 2024, 2), labels = ay(seq(2014, 2024, 2))) +
  labs(x = NULL, y = NULL) + thm + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.width = unit(0.9, "cm"))
ggsave(file.path(OUT, "fig3_deprivation_two_panel.png"), f3, width = 6.5, height = 3.8, dpi = 600, bg = "white")

sl <- regu[year %in% c(2014, 2024)][, yl := ay(year)]
# spread the 2024/25 labels so that none overlap (minimum gap of 3.2 thousand on the y axis)
lab <- sl[year == 2024][order(starts)]
lab[, y_lab := starts / 1000]
for (i in 2:nrow(lab)) if (lab$y_lab[i] - lab$y_lab[i - 1] < 3.2) lab$y_lab[i] <- lab$y_lab[i - 1] + 3.2
lab[, x_end := 2.08]
f4 <- ggplot(sl, aes(yl, starts / 1000, group = region)) +
  geom_line(colour = "grey40", linewidth = 0.5) + geom_point(colour = "black", size = 1.8) +
  geom_segment(data = lab, aes(x = 2.02, xend = x_end, y = starts / 1000, yend = y_lab), colour = "grey60", linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = lab, aes(x = x_end + 0.02, y = y_lab, label = region), hjust = 0, size = 2.6, family = "Arial", inherit.aes = FALSE) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.9))) + scale_y_continuous(limits = c(0, 85), expand = expansion(mult = 0)) +
  labs(x = NULL, y = "Starts (thousands)") + thm
ggsave(file.path(OUT, "fig4_region_slope.png"), f4, width = 6.5, height = 4.2, dpi = 600, bg = "white")

close(logf)
cat("\nREVISION ANALYSIS DONE. Outputs in", OUT, "\n")
