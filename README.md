# Who still starts an apprenticeship? Data and code

Replication package for the article *Who still starts an apprenticeship? Age, level, deprivation and the changing composition of apprenticeship participation in England, 2014/15 to 2024/25* (Research in Post-Compulsory Education, revision 1, October 2026).

Everything in the article can be reproduced from this package with R (>= 4.5) and the packages `data.table`, `ggplot2`, `sandwich`, `lmtest` and `readxl`.

## How to run

```r
Sys.setenv(RPCE_BASE = "/path/to/this/folder")      # or open the scripts from RStudio inside the folder
source("analysis/04a_build_population_denominators.R")   # population denominators by IMD quintile
source("analysis/04_revision_analysis.R")                # all revision-1 analyses, tables and figures
source("analysis/03_analysis.R")                         # the original (submitted) analysis, kept for reference
source("analysis/tests/test_revision_outputs.R")         # validation tests (all should PASS)
```

Outputs are written to `analysis/out_r1/` (revision 1) and `analysis/out/` (original). `checks_r1.txt` in `analysis/out_r1/` logs every validation check (splicing overlaps, cell sums, denominator comparisons).

## Data (all open licences; downloaded 16 July 2026 unless stated)

| Folder | Source | Files used |
|---|---|---|
| `data/cur/` | DfE, *Apprenticeships*, academic year 2025/26 (August 2025 to April 2026) release, Explore Education Statistics, Open Government Licence v3.0 | historical summary 2002/03 to 2024/25; learner deprivation (apps_11); learner detailed (apps_09); geography and population; national achievement rates by deprivation (na03); monthly starts (apps_24) and underlying monthly starts (apps_u10) |
| `data/hist/` | DfE, *Apprenticeships and traineeships*, 2022/23 release | learner deprivation starts; learner detailed; starts since 2009/10 by region; geography and population |
| `data/hist1920/` | DfE, *Apprenticeships and traineeships*, 2019/20 release | learner deprivation starts; learner detailed |
| `data/ext/` | ONS mid-2020 LSOA population estimates by single year of age (SAPE23DT2); ONS mid-year population estimates by region 2009 to 2024 (Nomis dataset NM_2002_1, downloaded 14 October 2026); MHCLG English Indices of Deprivation 2019, File 7 | denominators by IMD quintile and by region |

The `data-guidance.txt` files are the DfE's own metadata for each release. DfE counts are rounded to the nearest 10 (nearest 100 before 2014/15); `low` means fewer than 5 and `z` not applicable.

## Scripts

- `analysis/04a_build_population_denominators.R`: matches the ONS LSOA estimates to IMD 2019 deciles and writes `data/ext/imd2019_quintile_population.csv`.
- `analysis/04_revision_analysis.R`: splices the three DfE releases; within-level and age-by-level deprivation ratios and concentration indices; fixed-composition counterfactual; starts per 1,000 residents by quintile; interrupted time series models (level and slope changes in 2010/11, 2017/18 and 2020/21; Newey-West standard errors) and sensitivity models; regional rates and dispersion for 2014/15 to 2024/25; sex, ethnicity and LLDD composition; achievement rates; monthly timing and standards share; figures.
- `analysis/03_analysis.R`: the original analysis submitted in July 2026 (two-knot slope-only regressions on 2014/15 to 2024/25).
- `analysis/tests/`: automated checks on the outputs.

`supplementary_tables_R1.xlsx` collects the supplementary tables S1 to S9 cited in the article.

## Licence

Code: MIT License (see LICENSE). Data: Open Government Licence v3.0 (DfE, ONS, MHCLG); the ONS and MHCLG files are redistributed unchanged.
