########################################################################
# fig1_fpr_grid.R
# False positive rate (Beta = 0) vs Lambda. One figure per X_arch.
# Rows = number of species; columns = Y_arch (matched -> increasing
# mismatch); lines = tree(s) used for fitting.
########################################################################

source("plot_helpers.R")

Reps <- 1000   # replicates per combo (must match run_task.R)

plot_fpr_grid <- function(x_arch_value, save = TRUE) {

  df <- open_results(x_arch = x_arch_value) %>%
    filter(Beta == 0) %>%
    select(Nspec, Lambda, Y_arch, species_p, x_p, y_p, random_p, random500_p, unrelated_p)

  if (nrow(df) == 0) {
    warning("No rows for X_arch=", x_arch_value, " -- skipping")
    return(invisible(NULL))
  }

  df <- order_y_arch(df)

  p_cols <- scenario_p_cols   # defined in plot_helpers.R

  summ <- summarize_rejection_rate(df, c("Nspec", "Lambda", "Y_arch"), p_cols)
  summ$scenario <- factor(summ$scenario, levels = names(p_cols))

  # flag points not based on the full set of replicates
  short <- summ %>% filter(n_valid != Reps)
  if (nrow(short) > 0) {
    warning("X_arch=", x_arch_value, ": ", nrow(short), " points use fewer than ", Reps, " replicates")
  }

  g <- ggplot(summ, aes(x = Lambda, y = rejection_rate, color = scenario, linetype = scenario)) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey50") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = scenario_colors) +
    scale_linetype_manual(values = scenario_linetypes) +
    facet_grid(Nspec ~ Y_arch, labeller = labeller(
      Nspec  = function(x) paste(x, "species"),
      Y_arch = function(x) paste("Y:", x))) +
    guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    labs(
      x = expression(paste("Speciation rate (", lambda, ")")),
      y = "False positive rate (FPR)",
      color = NULL, linetype = NULL
    ) +
    base_theme

  if (save) {
    fname <- sprintf("figures/fig1_fpr_Xarch%s.png", x_arch_value)
    ggsave(fname, g, width = 14, height = 7, dpi = 300)
    cat("Saved", fname, "\n")
  }
  g
}

# Create all plots
for (xa in c(1, 2, 5, 10, 100)) plot_fpr_grid(x_arch_value = xa)