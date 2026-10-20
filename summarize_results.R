########################################################################
# summarize_results.R
# Reduces results_combined/ (one row per replicate) to small summary
# tables for the repository:
#
#   results/summary_by_method.csv
#     One row per parameter combination x fitting method:
#     rejection rate at alpha = 0.05 (FPR when Beta = 0, power otherwise),
#     and the mean, bias, SD and RMSE of beta-hat,
#     mean SE, 95% CI coverage and mean sigma.
#
#   results/summary_trees.csv
#     One row per parameter combination: mean RF distances and tree heights.
#
# Run after combine_results.R.
########################################################################

library(data.table)

alpha   <- 0.05
methods <- c("species", "x", "y", "random", "random500", "unrelated")
params  <- c("Nspec", "Lambda", "X_arch", "Y_arch", "Beta")

files <- list.files("results_combined", pattern = "^Nspec.*\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("No files in results_combined/ -- run combine_results.R first")
d <- rbindlist(lapply(files, fread, colClasses = list(character = "Y_arch")))
cat("Read", nrow(d), "replicates\n")

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
# Sort rows in a readable order (Y_arch and method in their natural order,
# not alphabetical) and round to 5 significant digits to keep files small.
tidy <- function(dt) {
  dt <- copy(dt)
  dt[, y_ord := match(Y_arch, c("same", "1", "2", "5", "10", "100"))]
  sort_by <- c("Nspec", "X_arch", "y_ord", "Beta", "Lambda")
  if ("method" %in% names(dt)) { dt[, m_ord := match(method, methods)]; sort_by <- c(sort_by, "m_ord") }
  setorderv(dt, sort_by)
  dt[, intersect(c("y_ord", "m_ord"), names(dt)) := NULL]
  num <- setdiff(names(dt)[sapply(dt, is.double)], params)
  dt[, (num) := lapply(.SD, signif, 5), .SDcols = num]
  dt
}
by_method <- tidy(by_method)
trees     <- tidy(trees)

dir.create("results", showWarnings = FALSE)
fwrite(by_method, "results/summary_by_method.csv")
fwrite(trees,     "results/summary_trees.csv")
cat("Wrote results/summary_by_method.csv (", nrow(by_method), " rows, ",
    round(file.size("results/summary_by_method.csv") / 1e6, 2), " MB)\n", sep = "")
cat("Wrote results/summary_trees.csv (", nrow(trees), " rows, ",
    round(file.size("results/summary_trees.csv") / 1e6, 2), " MB)\n", sep = "")