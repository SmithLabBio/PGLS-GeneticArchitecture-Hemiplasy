########################################################################
# combine_results.R
# Combines the per-combo CSVs from one simulation run into one compressed
# CSV per Nspec x X_arch, so figure scripts can load only the slice they
# need, and reports any combos that are missing, incomplete or duplicated.
#
# Usage:
#   Rscript combine_results.R main      # sim_output/         -> results_combined/
#   Rscript combine_results.R weak      # sim_output_weak/    -> results_combined_weak/
#   Rscript combine_results.R overlap   # sim_output_overlap/ -> results_combined_overlap/
#
# Each run's folders, parameters and grid are set in `runs` below; the grid
# must match the one in that run's run_task script. Safe to re-run: the
# output folder is deleted and rebuilt from whatever CSVs currently exist.
########################################################################

library(data.table)

# ---- settings for each run --------------------------------------------------
# params: columns that identify a parameter combination
# grid:   every combination the run should contain (same as its run_task script)
runs <- list(
  main = list(
    in_dir  = "sim_output",
    out_dir = "results_combined",
    params  = c("Nspec", "Lambda", "X_arch", "Y_arch", "Beta"),
    grid    = CJ(Lambda = 10^(-4:2), X_arch = c(1, 2, 5, 10, 100),
                 Y_arch = c("same", "1", "2", "5", "10", "100"),
                 Beta = c(0, 0.25, 0.5, 0.75, 1.0), Nspec = c(10, 100))
  ),
  weak = list(    # Pagel's lambda < 1 on X (run_task_weak.R)
    in_dir  = "sim_output_weak",
    out_dir = "results_combined_weak",
    params  = c("Nspec", "Lambda", "X_arch", "Y_arch", "Beta", "pagel_lambda_x"),
    grid    = CJ(Lambda = 10^(-4:2), X_arch = 1, Y_arch = c("same", "1", "100"),
                 Beta = c(0, 0.25, 0.5, 0.75, 1.0), Nspec = c(10, 100),
                 pagel_lambda_x = 0.5)
  ),
  overlap = list( # X and Y share n_shared of K = 10 genes (run_task_overlap.R)
    in_dir  = "sim_output_overlap",
    out_dir = "results_combined_overlap",
    params  = c("Nspec", "Lambda", "X_arch", "n_shared", "Beta"),
    grid    = CJ(Lambda = 10^(-4:2), X_arch = 10, n_shared = c(0, 2, 5, 8, 10),
                 Beta = c(0, 0.5), Nspec = c(10, 100))
  )
)

Reps <- 1000   # replicates per combo (same in every run_task script)

# columns every run must have, besides its own params
common_cols <- c("rep", "rep_seed",
                 "species_p", "x_p", "y_p", "random_p", "random500_p", "unrelated_p",
                 "species_sigma", "x_sigma", "y_sigma", "random_sigma",
                 "random500_sigma", "unrelated_sigma",
                 "RF_X", "RF_Y", "RF_Random", "RF_Random500", "RF_Unrelated",
                 "species_height", "height_x", "height_y")

# ---- which run ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% names(runs)) {
  stop("Usage: Rscript combine_results.R <run>, where <run> is one of: ",
       paste(names(runs), collapse = ", "))
}
run <- args[1]
cfg <- runs[[run]]
cat("Combining run '", run, "' from ", cfg$in_dir, "/\n", sep = "")

# ---- read -----------------------------------------------------------------
# Output files start with "Lambda" and end in .csv; unfinished files end in
# .csv.tmp and sessionInfo files in .txt, so neither is matched.
files <- list.files(cfg$in_dir, pattern = "^Lambda.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "result files\n")
if (length(files) == 0) stop("No result CSVs found in ", cfg$in_dir, "/")

# Read Y_arch as text where it exists; otherwise fread guesses integer for
# files without "same", so its type would depend on which files exist.
read_one <- function(f) {
  has_y <- "Y_arch" %in% names(fread(f, nrows = 0))
  if (has_y) fread(f, colClasses = list(character = "Y_arch")) else fread(f)
}
tables <- lapply(files, read_one)

combined <- rbindlist(tables, use.names = TRUE)
cat("Combined", nrow(combined), "rows,", ncol(combined), "columns\n")

params <- cfg$params
grid   <- copy(cfg$grid)

# ---- checks ---------------------------------------------------------------
# 1. expected columns present
missing_cols <- setdiff(c(params, common_cols), names(combined))
if (length(missing_cols) > 0) stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))

# round Lambda so 10^-4 from the grid matches 1e-04 read back from CSV
grid[, Lambda := signif(Lambda, 6)]
combined[, Lambda := signif(Lambda, 6)]

# 2. every parameter takes only values in the grid (catches results from the
#    wrong run, or a parameter that never reached the simulation)
for (p in params) {
  bad <- setdiff(unique(combined[[p]]), unique(grid[[p]]))
  if (length(bad) > 0) stop("Unexpected values of ", p, ": ", paste(bad, collapse = ", "))
}

# 3. every combo has exactly Reps replicates, and none is duplicated
counts <- combined[, .N, by = params]
short  <- counts[N != Reps]
if (nrow(short) > 0) {
  warning(nrow(short), " combos do not have ", Reps, " rows:")
  print(short)
}
dups <- combined[, .N, by = c(params, "rep")][N > 1]
if (nrow(dups) > 0) stop(nrow(dups), " duplicated (combo, rep) rows -- check for stale files in ", cfg$in_dir, "/")

# 4. grid combos with no results yet
not_run <- grid[!counts, on = params]
cat(nrow(grid) - nrow(not_run), "of", nrow(grid), "combos present\n")
missing_file <- sprintf("missing_combos_%s.csv", run)
if (nrow(not_run) > 0) {
  warning(nrow(not_run), " grid combos have no results file yet")
  fwrite(not_run, missing_file)
  cat("Listed in", missing_file, "\n")
} else if (file.exists(missing_file)) {
  file.remove(missing_file)   # don't leave an out-of-date list behind
}

# 5. failed fits (NA/NaN) in the p-values
na_counts <- sapply(combined[, .SD, .SDcols = patterns("_p$")], function(v) sum(!is.finite(v)))
if (any(na_counts > 0)) { warning("Non-finite p-values:"); print(na_counts[na_counts > 0]) }

# ---- write ----------------------------------------------------------------
# One compressed CSV per Nspec x X_arch. Deleted and rebuilt on every run.
unlink(cfg$out_dir, recursive = TRUE)
dir.create(cfg$out_dir)
for (g in split(combined, by = c("Nspec", "X_arch"))) {
  fwrite(g, file.path(cfg$out_dir, sprintf("Nspec%s_X%s.csv.gz", g$Nspec[1], g$X_arch[1])))
}
cat("Wrote", length(list.files(cfg$out_dir)), "files to", cfg$out_dir, "/\n")