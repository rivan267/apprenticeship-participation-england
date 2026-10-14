# =====================================================================
# Revision 1, step 0: population denominators by IMD 2019 quintile
#
# Builds data/ext/imd2019_quintile_population.csv from
#   (1) ONS mid-2020 LSOA population estimates by single year of age
#       (SAPE23DT2, 2011 LSOA boundaries), data/ext/ons_lsoa_mid2020_syoa.xlsx
#   (2) MHCLG English Indices of Deprivation 2019, File 7 (all ranks,
#       deciles, scores and mid-2015 population denominators),
#       data/ext/imd2019_file7.csv
# The DfE assigns each apprentice an IMD 2019 quintile from the home
# postcode; quintiles are fifths of the 32,844 English LSOAs ranked by
# IMD score (decile 1-2 = quintile 1, most deprived).
# Run this before 04_revision_analysis.R (Rscript or Source in RStudio).
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE <- Sys.getenv("RPCE_BASE", unset = "/Users/gestiardirrgmail.com/hairur new topik/RPCE_Apprenticeship_England")
if (!dir.exists(BASE) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
EXT <- file.path(BASE, "data", "ext")

# ---- (1) ONS LSOA single year of age, persons ----
xl <- file.path(EXT, "ons_lsoa_mid2020_syoa.xlsx")
ons <- as.data.table(read_excel(xl, sheet = "Mid-2020 Persons", skip = 4))
setnames(ons, make.names(names(ons)))
code_col <- grep("^LSOA.Code", names(ons), value = TRUE)[1]
age_cols <- grep("^X[0-9]+$|^X90\\.$", names(ons), value = TRUE)
age_num  <- as.integer(gsub("[^0-9]", "", age_cols))
stopifnot(length(age_cols) == 91, all(sort(age_num) == 0:90))
ons <- ons[grepl("^E", get(code_col))]                     # England LSOAs only
band <- function(a, b) rowSums(ons[, age_cols[age_num >= a & age_num <= b], with = FALSE])
lsoa <- data.table(lsoa11 = ons[[code_col]], all_ages = ons[[grep("^All.Ages", names(ons))[1]]],
                   pop_16_18 = band(16, 18), pop_19_24 = band(19, 24), pop_16_24 = band(16, 24),
                   pop_25_64 = band(25, 64), pop_16_64 = band(16, 64), pop_0_15 = band(0, 15))
stopifnot(nrow(lsoa) == 32844)
fwrite(lsoa, file.path(EXT, "ons_lsoa_mid2020_agebands.csv"))

# ---- (2) IMD 2019 File 7: decile and mid-2015 denominators ----
imd <- fread(file.path(EXT, "imd2019_file7.csv"))
dec_col <- grep("^Index of Multiple Deprivation \\(IMD\\) Decile", names(imd), value = TRUE)[1]
pop_cols <- grep("mid 2015|Working age population", names(imd), value = TRUE)
imd <- imd[, c("LSOA code (2011)", dec_col, pop_cols), with = FALSE]
setnames(imd, c("lsoa11", "imd_decile", paste0("imd2015_", pop_cols)))
imd[, imd_quintile := (imd_decile + 1L) %/% 2L]
stopifnot(nrow(imd) == 32844)

# ---- (3) merge and aggregate ----
m <- merge(lsoa, imd, by = "lsoa11")
stopifnot(nrow(m) == 32844)
q <- m[, c(list(n_lsoa = .N), lapply(.SD, sum)), by = imd_quintile,
       .SDcols = c("all_ages", "pop_16_18", "pop_19_24", "pop_16_24", "pop_25_64", "pop_16_64", grep("^imd2015_", names(m), value = TRUE))][order(imd_quintile)]
fwrite(q, file.path(EXT, "imd2019_quintile_population.csv"))
cat("Quintile population shares (mid-2020):\n")
print(q[, .(imd_quintile, n_lsoa, share_16_18 = round(100 * pop_16_18 / sum(pop_16_18), 1), share_19_24 = round(100 * pop_19_24 / sum(pop_19_24), 1),
            share_16_24 = round(100 * pop_16_24 / sum(pop_16_24), 1), share_25_64 = round(100 * pop_25_64 / sum(pop_25_64), 1),
            share_16_64 = round(100 * pop_16_64 / sum(pop_16_64), 1))])
cat("Written:", file.path(EXT, "imd2019_quintile_population.csv"), "\n")
