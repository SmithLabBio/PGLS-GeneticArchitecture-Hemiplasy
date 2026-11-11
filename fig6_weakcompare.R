########################################################################
# fig_weakcompare.R
# Effect of weaker phylogenetic signal in X (Pagel's lambda = 0.5) on
# false positive rate and power, compared with the main run (lambda = 1).
#
# X is always a single gene (X_arch = 1). Both figures share one layout:
#   rows    = number of species
#   columns = Y architecture (same as X, 1 gene, 100 genes)
#   x       = speciation rate (Lambda), log scale
#   colour  = tree(s) used for fitting (a subset of methods, below)
#   linetype= Pagel's lambda on X: solid = 1 (main run), dashed = 0.5
#
# Figures:
#   figures/fig6_weak_fpr.png           FPR (Beta = 0)
#   figures/fig6_weak_tpr_beta<b>.png   power at one Beta (default 0.5)
#
# Needs results_combined/ (main run) and results_combined_weak/ (from
# combine_results_weak.R).
########################################################################

source("plot_helpers.R")

Reps <- 1000   # replicates per combo (must match run_task.R / run_task_weak.R)

# Fitting methods to show (display names from plot_helpers.R). Dropped:
# "Random architecture" and "500 gene trees, true species tree".
weak_methods <- c("Y's true tree(s)", "X's true tree(s)")
weak_y_arch  <- c("same", "1", "100")   # Y architectures in the weak run

# --- load: main run (lambda = 1) + weak run (lambda = 0.5), X_arch = 1 only ---
load_weak_comparison <- function() {
  p_cols <- scenario_p_cols[weak_methods]
  keep   <- c("Nspec", "Lambda", "Y_arch", "Beta", unname(p_cols))

  main <- open_results(x_arch = 1) %>%
    filter(Y_arch %in% weak_y_arch) %>%
    select(all_of(keep)) %>%
    mutate(pagel_lambda_x = 1)

  weak <- open_results(x_arch = 1, path = "results_combined_weak") %>%
    select(all_of(c(keep, "pagel_lambda_x")))

  bind_rows(main, weak)
}

plot_weak <- function(d, beta_value, save = TRUE) {
  p_cols <- scenario_p_cols[weak_methods]
  is_fpr <- beta_value == 0

  df <- d %>% filter(Beta == beta_value)
  if (nrow(df) == 0) stop("No rows for Beta = ", beta_value)
  df <- order_y_arch(df)

  summ <- summarize_rejection_rate(df, c("Nspec", "Lambda", "Y_arch", "pagel_lambda_x"), p_cols)
  summ <- summ %>%
    mutate(scenario = factor(scenario, levels = weak_methods),
           pagel = factor(pagel_lambda_x, levels = c(1, 0.5)))

  short <- summ %>% filter(n_valid != Reps)
  if (nrow(short) > 0) {
    warning("Beta=", beta_value, ": ", nrow(short), " points use fewer than ", Reps, " replicates")
  }

  g <- ggplot(summ, aes(x = Lambda, y = rejection_rate,
                        color = scenario, linetype = pagel,
                        group = interaction(scenario, pagel))) +
    { if (is_fpr) geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey50") } +
    geom_line(linewidth = 1) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = scenario_colors[weak_methods]) +
    scale_linetype_manual(values = c("1" = "solid", "0.5" = "dashed"),
                          name = expression(paste("Pagel's ", lambda, " on X"))) +
    facet_grid(Nspec ~ Y_arch, labeller = labeller(
      Nspec  = function(x) paste(x, "species"),
      Y_arch = function(x) paste("Y:", x))) +
    guides(color = guide_legend(nrow = 2, order = 1),
           linetype = guide_legend(order = 2, override.aes = list(color = "grey30"))) +
    labs(
      x = expression(paste("Speciation rate (", lambda, ")")),
      y = if (is_fpr) "False positive rate" else bquote("True positive rate (" * beta == .(beta_value) * ")"),
      color = NULL
    ) +
    base_theme

  if (save) {
    fname <- if (is_fpr) "figures/fig6_weak_fpr.png" else sprintf("figures/fig6_weak_tpr_beta%s.png", beta_value)
    ggsave(fname, g, width = 12, height = 7, dpi = 300)
    cat("Saved", fname, "\n")
  }
  g
}

d <- load_weak_comparison()
plot_weak(d, beta_value = 0)     # FPR
plot_weak(d, beta_value = 0.5)   # power