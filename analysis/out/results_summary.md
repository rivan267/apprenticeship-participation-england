# Results summary: Who still starts an apprenticeship? (England, 2014/15 to 2024/25)

Analysis run: 16 July 2026. All inputs are DfE Explore Education Statistics underlying CSVs (Open Government Licence v3.0), downloaded 16 July 2026. Full pipeline: `analysis/03_analysis.R`; integrity tests: `analysis/tests/test_analysis_outputs.R` (5 tests, all passing). DfE counts are rounded to the nearest 10.

## Data sources and validation

| Component | File | Coverage used |
|---|---|---|
| Age x level (national) | `cur/data/app-historical-summary-to-2425.csv` (single consistent DfE historical summary) | 2009/10 to 2024/25 (analysis window 2014/15 to 2024/25) |
| Deprivation (IMD quintiles) | 2019-20 release + 2022-23 release + 2025/26 release, spliced | 2014/15 to 2024/25 |
| Regions | `app-starts-since-202223-q4.csv` + `app-geography-population-202526-q3.csv` | 2014/15 to 2024/25 (rates 2017/18+) |

Validation results: (1) age x level cells reproduce published national totals in every year (max difference 10, i.e. rounding); (2) deprivation splice: overlapping years agree exactly between releases (2017/18 to 2019/20 across releases 2019-20 vs 2022-23; 2020/21 to 2022/23 across 2022-23 vs 2025/26); (3) regional splice: all 30 overlapping year x region cells agree exactly (max difference 0). See `checks_overlap.txt`.

## RQ1. Volume and age x level composition

Total starts: 499,890 (2014/15) -> 375,760 (2017/18, first levy year) -> 322,530 (2019/20) -> 321,440 (2020/21, COVID floor) -> 353,500 (2024/25). Net decade change: -146,390 (-29.3%).

Composition shift (share of all starts):

| Dimension | 2014/15 | 2024/25 |
|---|---|---|
| Intermediate (level 2) | 59.7% | 18.6% |
| Advanced (level 3) | 36.4% | 41.6% |
| Higher (levels 4+) | 4.0% | 39.8% |
| Under 19 | 25.2% | 21.2% |
| 19-24 | 32.0% | 27.5% |
| 25+ | 42.8% | 51.3% |

Cell decomposition of the net change (-146,390), `tab_decomp_change.csv`:
- Intermediate cells account for -232,600 (159% of the net decline): 25+ intermediate -101,450; 19-24 intermediate -74,920; under-19 intermediate -56,230.
- Higher apprenticeships offset +120,960: 25+ higher +86,750; 19-24 higher +28,320; under-19 higher only +5,890.
- In 2024/25, 25+ starts at levels 4+ were 101,160; under-19 starts at levels 4+ were 7,030.

## RQ2. Deprivation gradient (learner home IMD quintile, known-quintile starts)

| Year | Q1 most deprived | Q5 least deprived | Ratio Q1:Q5 | Concentration index |
|---|---|---|---|---|
| 2014/15 | 128,290 | 69,950 | 1.83 | -0.110 |
| 2016/17 | 127,770 | 69,180 | 1.85 | -0.114 |
| 2017/18 | 83,810 | 61,410 | 1.36 | -0.054 |
| 2019/20 | 66,300 | 58,040 | 1.14 | -0.022 |
| 2024/25 | 68,720 | 66,240 | 1.04 | -0.003 |

The pro-deprived gradient collapsed, and the collapse is timed exactly at the levy year: the ratio fell from 1.85 (2016/17) to 1.36 (2017/18) in a single year, then drifted to near-parity. By 2024/25 the most deprived quintile no longer supplies the largest share of starts. Within-group tables (`dep_within_age.csv`, `dep_within_level.csv`) allow the manuscript to show this is not purely a level-composition artefact.

## RQ3. Regions

Change in starts, 2014/15 to 2024/25 (`tab_region_counts.csv`):

| Region | 2014/15 | 2024/25 | Change |
|---|---|---|---|
| North East | 35,220 | 18,580 | -47.2% |
| Yorkshire and The Humber | 62,550 | 38,490 | -38.5% |
| West Midlands | 61,240 | 38,720 | -36.8% |
| North West | 79,310 | 50,580 | -36.2% |
| East Midlands | 48,060 | 30,900 | -35.7% |
| South West | 51,480 | 39,790 | -22.7% |
| East of England | 45,790 | 37,010 | -19.2% |
| South East | 65,030 | 55,260 | -15.0% |
| London | 45,550 | 40,380 | -11.4% |

A clear north-south gradient: the regions that historically used apprenticeships most intensively lost the most. Rates per 100,000 population available 2017/18+ in `tab_region_rates_per100k.csv`.

## RQ4. Segmented trends (knots: 2017/18 opens the post-levy segment; 2020/21 opens the post-COVID segment)

Annual slopes in starts per year (`tab2_segmented.csv`; 11 annual points, descriptive):

| Series | Pre-levy slope | Post-levy slope | Post-COVID slope |
|---|---|---|---|
| Total | -15,634 | -51,825 | +6,079 |
| Under 19 | -1,196 | -17,106 | +511 |
| 19-24 | -11,833 | -13,653 | +245 |
| 25+ | -2,606 | -21,068 | +5,324 |
| Intermediate | -29,833 | -51,563 | -4,828 |
| Advanced | +6,210 | -16,957 | +769 |
| Higher | +8,027 | +16,690 | +10,139 |

Reading: before the levy the total was roughly stable (slope not significant); the decline concentrates in the post-levy window 2017/18 to 2019/20; the post-COVID recovery (+6,079 per year overall) is carried almost entirely by 25+ and higher-level starts, while under-19 (+511 per year) and 19-24 (+245 per year) starts are essentially flat at their post-levy floor.

## Headline answer to the title question

Who still starts an apprenticeship in England? Increasingly: adults aged 25+, at levels 4 and above, from mid- and less-deprived neighbourhoods, in the South. The under-19 intermediate apprentice from a deprived area, the archetype of the pre-levy system, has largely disappeared: under-19 intermediate starts fell from 85,600 (2014/15) to 29,370 (2024/25).

## Caveats to carry into the manuscript

1. Aggregate published data; the segmented regressions describe timing, not causal effects; the 2017/18 break also bundles frameworks-to-standards migration, minimum duration and co-investment.
2. IMD vintage changes across releases (IMD 2015 vs IMD 2019); verify the vintage note in each release's metadata during drafting; shares use known-quintile starts only.
3. 2020/21 is a pandemic year; robustness excluding it is planned in the manuscript.
4. DfE rounding to nearest 10; "Unknown" level category (80 starts, 2014/15 only) retained in totals.
