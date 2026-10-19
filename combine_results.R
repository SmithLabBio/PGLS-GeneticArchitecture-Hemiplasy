########################################################################
# combine_results.R
# Combines all per-combo CSVs in sim_output/ into one compressed CSV per
# Nspec x X_arch in results_combined/, so figure scripts can load only
# the slice they need.
#
# Safe to re-run: it deletes and rebuilds results_combined/ from whatever
# CSVs currently exist, and reports any combos that are missing or incomplete.
########################################################################

library(data.table)

Reps <- 1000   # must match run_task.R

# Output files are named Lambda<>_X<>_Y<>_Beta<>_N<>.csv. Unfinished files end
# in .csv.tmp, so they are not matched.
files <- list.files("sim_output", pattern = "^Lambda.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "result files\n")
if (length(files) == 0) stop("No result CSVs found in sim_output/")

# Read Y_arch as text in every file. Otherwise fread guesses integer for
# files without "same", so its type would depend on which files exist.
combined <- rbindlist(lapply(files, fread, colClasses = list(character = "Y_arch")))
cat("Combined", nrow(combined), "rows,", ncol(combined), "columns\n")

# ---- checks ---------------------------------------------------------------
params <- c("Nspec", "Lambda", "X_arch", "Y_arch", "Beta")

# 1. columns the figure scripts expect
expected <- c(params, "rep", "rep_seed",
              "species_p", "x_p", "y_p", "random_p", "random500_p", "unrelated_p",
              "species_sigma", "x_sigma", "y_sigma", "random_sigma",
              "random500_sigma", "unrelated_sigma",
              "RF_X", "RF_Y", "RF_Random", "RF_Random500", "RF_Unrelated",
              "species_height", "height_x", "height_y")
missing_cols <- setdiff(expected, names(combined))
if (length(missing_cols) > 0) stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))

# 2. every combo has exactly Reps replicates, and none is duplicated
counts <- combined[, .N, by = params]
short  <- counts[N != Reps]
if (nrow(short) > 0) {
  warning(nrow(short), " combos do not have ", Reps, " rows:")
  print(short)
}
dups <- combined[, .N, by = c(params, "rep")][N > 1]
if (nrow(dups) > 0) stop(nrow(dups), " duplicated (combo, rep) rows -- check for stale files in sim_output/")

# 3. combos in the grid with no results yet (same grid as run_task.R)
grid <- CJ(Lambda = 10^(-4:2), X_arch = c(1, 2, 5, 10, 100),
           Y_arch = c("same", "1", "2", "5", "10", "100"),
           Beta = c(0, 0.25, 0.5, 0.75, 1.0), Nspec = c(10, 100))
# round Lambda so 10^-4 from the grid matches 1e-04 read back from CSV
grid[, Lambda := signif(Lambda, 6)]
combined[, Lambda := signif(Lambda, 6)]
counts <- combined[, .N, by = params]   # recount after rounding Lambda
not_run <- grid[!counts, on = params]
cat(nrow(grid) - nrow(not_run), "of", nrow(grid), "combos present\n")
if (nrow(not_run) > 0) {
  warning(nrow(not_run), " grid combos have no results file yet")
  fwrite(not_run, "missing_combos.csv")
  cat("Listed in missing_combos.csv\n")
}

# 4. any failed fits (NA/NaN) in the p-values
na_counts <- sapply(combined[, .SD, .SDcols = patterns("_p$")], function(v) sum(!is.finite(v)))
if (any(na_counts > 0)) { warning("Non-finite p-values:"); print(na_counts[na_counts > 0]) }

# ---- write ----------------------------------------------------------------
# One compressed CSV per Nspec x X_arch. Deleted and rebuilt on every run.
out_dir <- "results_combined"
unlink(out_dir, recursive = TRUE)
dir.create(out_dir)
for (g in split(combined, by = c("Nspec", "X_arch"))) {
  fwrite(g, file.path(out_dir, sprintf("Nspec%s_X%s.csv.gz", g$Nspec[1], g$X_arch[1])))
}
cat("Wrote", length(list.files(out_dir)), "files to", out_dir, "/\n")