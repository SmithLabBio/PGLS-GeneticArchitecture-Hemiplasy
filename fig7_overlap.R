########################################################################
# fig7_overlap.R
# Partial gene overlap: X and Y each have K = 10 gene trees, of which
# n_shared are shared. Same layout as fig1/fig2:
#   rows    = number of species
#   columns = genes shared by X and Y (0 = separate genes ... 10 = identical)
#   x       = speciation rate (Lambda), log scale
#   colour  = tree(s) used for fitting
#
# Figures:
#   figures/fig7_overlap_fpr.png          FPR (Beta = 0)
#   figures/fig7_overlap_tpr_beta0.5.png  power (Beta = 0.5)
#
# The end columns repeat main-run conditions (0 shared = X 10 genes /
# Y 10 genes; 10 shared = X 10 genes / Y same as X), so they should agree
# with fig1/fig2 for X_arch = 10 to within simulation noise.
#
# Needs results_combined_overlap/ (Rscript combine_results.R overlap).
########################################################################

source("plot_helpers.R")

Reps <- 1000   # replicates per combo (must match run_task_overlap.R)
K    <- 10     # genes per trait in this run

plot_overlap <- function(beta_value, save = TRUE) {
  p_cols <- scenario_p_cols   # all six methods, names/colours from plot_helpers.R
  is_fpr <- beta_value == 0

  df <- open_results(path = "results_combined_overlap") %>%
    filter(Beta == beta_value) %>%
    select(all_of(c("Nspec", "Lambda", "n_shared", unname(p_cols))))
  if (nrow(df) == 0) stop("No rows for Beta = ", beta_value)

  summ <- summarize_rejection_rate(df, c("Nspec", "Lambda", "n_shared"), p_cols) %>%
    mutate(scenario = factor(scenario, levels = names(p_cols)),
           # column order: none shared -> all shared
           n_shared = factor(n_shared, levels = sort(unique(n_shared))))

  short <- summ %>% filter(n_valid != Reps)
  if (nrow(short) > 0) {
    warning("Beta=", beta_value, ": ", nrow(short), " points use fewer than ", Reps, " replicates")
  }

  g <- ggplot(summ, aes(x = Lambda, y = rejection_rate, color = scenario, linetype = scenario)) +
    { if (is_fpr) geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey50") } +
    geom_line(linewidth = 1) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = scenario_colors) +
    scale_linetype_manual(values = scenario_linetypes) +
    facet_grid(Nspec ~ n_shared, labeller = labeller(
      Nspec    = function(x) paste(x, "species"),
      n_shared = function(x) sprintf("%s of %s genes shared", x, K))) +
    guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    labs(
      x = expression(paste("Speciation rate (", lambda, ")")),
      y = if (is_fpr) "False positive rate" else bquote("True positive rate (" * beta == .(beta_value) * ")"),
      color = NULL, linetype = NULL
    ) +
    base_theme

  if (save) {
    fname <- if (is_fpr) "figures/fig7_overlap_fpr.png" else sprintf("figures/fig7_overlap_tpr_beta%s.png", beta_value)
    ggsave(fname, g, width = 15, height = 7, dpi = 300)
    cat("Saved", fname, "\n")
  }
  g
}

plot_overlap(beta_value = 0)     # FPR
plot_overlap(beta_value = 0.5)   # power