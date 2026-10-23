########################################################################
# summarize_results.R
# Reduces one run's combined results (one row per replicate) to small
# summary tables for the repository:
#
#   results/<prefix>summary_by_method.csv
#     One row per parameter combination x fitting method:
#     rejection rate at alpha = 0.05 (FPR when Beta = 0, power otherwise),
#     and the mean, bias, SD and RMSE of beta-hat,
#     mean SE, 95% CI coverage and mean sigma.
#
#   results/<prefix>summary_trees.csv
#     One row per parameter combination: mean RF distances and tree heights.
#
# Usage:
#   Rscript summarize_results.R main      # results_combined/         -> results/summary_*.csv
#   Rscript summarize_results.R weak      # results_combined_weak/    -> results/weak_summary_*.csv
#   Rscript summarize_results.R overlap   # results_combined_overlap/ -> results/overlap_summary_*.csv
#
# Run after `Rscript combine_results.R <run>`.
########################################################################

library(data.table)

# ---- settings for each run --------------------------------------------------
# params: columns that identify a parameter combination, in the order rows
#         are sorted (same parameters as in combine_results.R)
runs <- list(
  main    = list(in_dir = "results_combined",         prefix = "",
                 params = c("Nspec", "X_arch", "Y_arch", "Beta", "Lambda")),
  weak    = list(in_dir = "results_combined_weak",    prefix = "weak_",
                 params = c("Nspec", "X_arch", "Y_arch", "pagel_lambda_x", "Beta", "Lambda")),
  overlap = list(in_dir = "results_combined_overlap", prefix = "overlap_",
                 params = c("Nspec", "X_arch", "n_shared", "Beta", "Lambda"))
)

alpha   <- 0.05
methods <- c("species", "x", "y", "random", "random500", "unrelated")
y_arch_order <- c("same", "1", "2", "5", "10", "100")   # natural, not alphabetical

# ---- which run ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% names(runs)) {
  stop("Usage: Rscript summarize_results.R <run>, where <run> is one of: ",
       paste(names(runs), collapse = ", "))
}
run    <- args[1]
cfg    <- runs[[run]]          # this run's settings (in_dir, prefix, params)
params <- cfg$params

# ---- read -----------------------------------------------------------------
files <- list.files(cfg$in_dir, pattern = "^Nspec.*\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("No files in ", cfg$in_dir, "/ -- run `Rscript combine_results.R ", run, "` first")

# Read Y_arch as text where it exists (the overlap run has no Y_arch column)
read_one <- function(f) {
  if ("Y_arch" %in% names(fread(f, nrows = 0))) fread(f, colClasses = list(character = "Y_arch")) else fread(f)
}
d <- rbindlist(lapply(files, read_one))
cat("Read", nrow(d), "replicates from", cfg$in_dir, "/\n")

missing <- setdiff(params, names(d))
if (length(missing) > 0) stop("Missing parameter columns: ", paste(missing, collapse = ", "))

# ---- per method -----------------------------------------------------------
# Reshape to one row per replicate x method, with generic column names.
long <- rbindlist(lapply(methods, function(m) {
  cbind(d[, params, with = FALSE],
        data.table(method = m,
                   beta  = d[[paste0(m, "_beta")]],
                   se    = d[[paste0(m, "_se")]],
                   p     = d[[paste0(m, "_p")]],
                   sigma = d[[paste0(m, "_sigma")]]))
}))

by_method <- long[, {
  ok   <- is.finite(p)
  rate <- mean(p[ok] < alpha)
  tcrit <- qt(0.975, df = Nspec[1] - 2)   # same df as fit_pgls(): N - 2 coefficients
  list(
    n            = .N,
    n_valid      = sum(ok),
    reject_rate  = rate,
    mean_beta    = mean(beta[ok]),
    bias_beta    = mean(beta[ok]) - Beta[1],
    sd_beta      = sd(beta[ok]),
    rmse_beta    = sqrt(mean((beta[ok] - Beta[1])^2)),
    mean_se      = mean(se[ok]),
    coverage_95  = mean(abs(beta[ok] - Beta[1]) <= tcrit * se[ok]),
    mean_sigma   = mean(sigma[ok])
  )
}, by = c(params, "method")]

# ---- per combination: trees -----------------------------------------------
tree_cols <- intersect(c("RF_X", "RF_Y", "RF_Random", "RF_Random500", "RF_Unrelated",
                         "species_height", "height_x", "height_y", "random_arch"), names(d))
trees <- d[, c(list(n = .N), lapply(.SD, function(v) mean(v, na.rm = TRUE))),
           by = params, .SDcols = tree_cols]
setnames(trees, tree_cols, paste0("mean_", tree_cols))

# ---- write ------------------------------------------------------------------
# Sort rows by params in the order given (Y_arch and method in their natural
# order, not alphabetical) and round to 5 significant digits.
tidy <- function(dt) {
  dt <- copy(dt)
  sort_by <- params
  if ("Y_arch" %in% params) {
    dt[, y_ord := match(Y_arch, y_arch_order)]
    sort_by[sort_by == "Y_arch"] <- "y_ord"
  }
  if ("method" %in% names(dt)) {
    dt[, m_ord := match(method, methods)]
    sort_by <- c(sort_by, "m_ord")
  }
  setorderv(dt, sort_by)
  helpers <- intersect(c("y_ord", "m_ord"), names(dt))   # temporary sort columns
  if (length(helpers) > 0) dt[, (helpers) := NULL]
  setcolorder(dt, params)   # parameter columns first, in the order above
  num <- setdiff(names(dt)[sapply(dt, is.double)], params)
  dt[, (num) := lapply(.SD, signif, 5), .SDcols = num]
  dt
}
by_method <- tidy(by_method)
trees     <- tidy(trees)

dir.create("results", showWarnings = FALSE)
out_method <- sprintf("results/%ssummary_by_method.csv", cfg$prefix)
out_trees  <- sprintf("results/%ssummary_trees.csv", cfg$prefix)
fwrite(by_method, out_method)
fwrite(trees,     out_trees)
cat("Wrote ", out_method, " (", nrow(by_method), " rows, ",
    round(file.size(out_method) / 1e6, 2), " MB)\n", sep = "")
cat("Wrote ", out_trees, " (", nrow(trees), " rows, ",
    round(file.size(out_trees) / 1e6, 2), " MB)\n", sep = "")
