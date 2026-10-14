# Validation tests for the revision-1 outputs (run after 04a and 04): Rscript analysis/tests/test_revision_outputs.R
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("RPCE_BASE", unset = "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England")
OUT <- file.path(BASE, "analysis", "out_r1")
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS " else "FAIL ", msg, "\n"); if (!isTRUE(cond)) stop(msg) }
need <- c("tab_segmented_main_modelM.csv", "tab_segmented_sensitivity.csv", "tab_within_ratio_selected_years.csv", "tab_dep_counterfactual_fixed_mix.csv",
          "tab_rates_per_1000_by_quintile.csv", "tab_region_dispersion_full_decade.csv", "tab_sex_ethnicity_lldd_by_year.csv",
          "tab_achievement_rates_2024_25.csv", "tab_standards_levy_share_by_year.csv", "tab_may_july_vs_aug_april.csv",
          "fig1_total_starts_2002_2024.png", "fig2_level.png", "fig3_deprivation_two_panel.png", "fig4_region_slope.png")
ok(all(file.exists(file.path(OUT, need))), "all revision outputs exist")
m <- fread(file.path(OUT, "tab_segmented_main_modelM.csv"))
ok(nrow(m) == 7 && all(m$n_obs == 23), "Model M: 7 series on 23 observations")
ok(abs(m[series == "Total", level_levy] - (-111321)) < 100, "Model M: total level change at levy about -111,300")
ok(m[series == "Under 19", post_covid_lo] < 0 && m[series == "Under 19", post_covid_hi] > 0, "Model M: under-19 post-pandemic slope CI includes zero")
r <- fread(file.path(OUT, "tab_rates_per_1000_by_quintile.csv"))
ok(abs(r[age == "Under 25" & year == 2024, rate_ratio_Q1_Q5] - 0.71) < 0.01, "Rates: under-25 Q1:Q5 rate ratio 0.71 in 2024/25")
ok(abs(r[age == "Total" & year == 2014, rate_ratio_Q1_Q5] - 1.66) < 0.01, "Rates: all-age Q1:Q5 rate ratio 1.66 in 2014/15")
cf <- fread(file.path(OUT, "tab_dep_counterfactual_fixed_mix.csv"))
ok(abs(cf[year == 2024, cf_ratio_fixed_2014_mix] - 1.41) < 0.01 && abs(cf[year == 2024, observed_ratio] - 1.04) < 0.01, "Counterfactual: 1.41 vs observed 1.04 in 2024/25")
w <- fread(file.path(OUT, "tab_within_ratio_selected_years.csv"))
ok(abs(w[group == "Level Intermediate", `2024/25`] - 1.44) < 0.01, "Within level: intermediate ratio 1.44 in 2024/25")
d <- fread(file.path(OUT, "tab_region_dispersion_full_decade.csv"))
ok(abs(d[year == 2014, max_min_ratio] - 2.74) < 0.01 && abs(d[year == 2024, max_min_ratio] - 1.76) < 0.01, "Regions: max/min ratio 2.74 -> 1.76")
s <- fread(file.path(OUT, "tab_standards_levy_share_by_year.csv"))
ok(abs(s[year == 2016, standards_pct] - 5.0) < 0.1 && abs(s[year == 2017, standards_pct] - 43.6) < 0.1, "Standards share 5.0% (2016/17) and 43.6% (2017/18)")
cat("ALL REVISION TESTS PASSED\n")
