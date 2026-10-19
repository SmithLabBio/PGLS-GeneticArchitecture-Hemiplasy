########################################################################
# fig3_rf_qc.R
# QC/diagnostic: normalized Robinson-Foulds distance between gene trees
# and the species tree vs Lambda (mean +/- SD across replicates).
# One figure: columns = number of species.
#
# Uses RF_X at a single X_arch. RF_X is the MEAN RF over X's gene trees,
# so its mean does not depend on X_arch, but its SD does: with X_arch = 1
# the band shows variation among individual gene trees; with X_arch = 100
# it shows variation in a 100-tree average (about 10x narrower).
#
# Reference line: RF to gene trees from an unrelated species tree (the
# level expected with no shared history).
#
# Also prints a cross-check that RF_X, RF_Y, RF_Random and RF_Random500
# have the same mean at each Lambda x Nspec, since they are all gene trees
# simulated on the same species tree.
########################################################################

source("plot_helpers.R")

plot_rf_qc <- function(x_arch_value = 100, save = TRUE) {

  df <- open_results(x_arch = x_arch_value) %>%
      select(Nspec, Lambda, RF_X, RF_Y, RF_Random, RF_Random500, RF_Unrelated)

  if (nrow(df) == 0) {
    warning("No rows for X_arch=", x_arch_value, " -- skipping")
    return(invisible(NULL))
  }

  summ <- df %>%
    group_by(Nspec, Lambda) %>%
    summarize(n = n(),
      mean_rf = mean(RF_X, na.rm = TRUE),
              sd_rf = sd(RF_X, na.rm = TRUE),
              mean_rf_unrelated = mean(RF_Unrelated, na.rm = TRUE),
              .groups = "drop")


  g <- ggplot(summ, aes(x = Lambda, y = mean_rf)) +
    geom_ribbon(aes(ymin = pmax(mean_rf - sd_rf, 0), ymax = pmin(mean_rf + sd_rf, 1)),
                fill = "#0072B2", alpha = 0.15) +
    geom_line(color = "#0072B2", linewidth = 0.8) +
    geom_line(aes(y = mean_rf_unrelated), color = "grey40", linetype = "dashed") +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(~ Nspec, labeller = labeller(Nspec = function(x) paste(x, "species"))) +
    labs(
      x = expression(paste("Speciation rate (", lambda, ")")),
      y = expression("Normalized RF distance (mean" %+-% "SD)")
    ) +
    base_theme

  if (save) {
    fname <- sprintf("figures/fig3_rf_qc_Xarch%s.png", x_arch_value)
    ggsave(fname, g, width = 10, height = 6, dpi = 300)
    cat("Saved", fname, "\n")
  }
  g
}

plot_rf_qc()