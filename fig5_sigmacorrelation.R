########################################################################
# fig5_sigmacorrelation.R
# Shows that tree height (which varies by ~8 orders of magnitude across
# Lambda, because trees are not rescaled) does not affect inference.
# Uses the correctly specified fit only: PGLS on Y's true gene tree(s).
#
#   (a) Residual variance sigma^2 vs mean height of Y's gene trees.
#       Y's residual evolves by Brownian motion (rate 1) on those trees, so
#       sigma^2 should equal height_y * (N - 2) / N (ML divides by N):
#       tree height is absorbed entirely into sigma^2.
#   (b) False positive rate (Beta = 0) vs the same heights. If height
#       mattered for the p-values, this would drift away from 0.05.
#
# One point per parameter combination (mean over replicates); color =number of species.
########################################################################

source("plot_helpers.R")
library(patchwork)

alpha <- 0.05

plot_tree_height <- function(save = TRUE) {

  df <- open_results() %>% filter(Beta == 0)

  summ <- df %>%
    group_by(Nspec, Lambda, X_arch, Y_arch) %>%
    summarize(n = sum(is.finite(y_p)),
              height_y = mean(height_y),
              sigma2 = mean(y_sigma^2, na.rm = TRUE),
              fpr = mean(y_p < alpha, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Nspec = factor(Nspec))

  species_colors <- c("10" = "#0072B2", "100" = "#D55E00")

  # (a) expected sigma^2 for a correct fit: height_y * (N - 2) / N
  hrange <- range(summ$height_y)
  ref <- expand.grid(Nspec = levels(summ$Nspec), height_y = hrange) %>%
    mutate(sigma2 = height_y * (as.numeric(as.character(Nspec)) - 2) /
                    as.numeric(as.character(Nspec)))

  pa <- ggplot(summ, aes(x = height_y, y = sigma2, color = Nspec)) +
    geom_line(data = ref, aes(group = Nspec), linetype = "dashed",
              show.legend = FALSE) +   # point-only key, so the two legends merge
    geom_point(size = 1.6, alpha = 0.7) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = species_colors, name = "Species") +
    labs(x = "Mean height of Y's gene trees",
         y = expression(paste("Residual variance (", sigma^2, ")")),
         title = "(a)") +
    base_theme

  # (b) false positive rate, with the 95% range expected from simulation
  #     noise alone around 0.05 (binomial, using the typical replicate count)
  n_typ <- median(summ$n)
  band  <- 1.96 * sqrt(alpha * (1 - alpha) / n_typ)

  pb <- ggplot(summ, aes(x = height_y, y = fpr, color = Nspec)) +
    annotate("rect", xmin = hrange[1], xmax = hrange[2],
             ymin = max(0, alpha - band), ymax = alpha + band, fill = "grey85") +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "grey40") +
    geom_point(size = 1.6, alpha = 0.7) +
    scale_x_log10() +
    # zoom with coord_cartesian (not scale limits) so nothing is dropped
    coord_cartesian(ylim = c(0, max(0.15, max(summ$fpr), alpha + band))) +
    scale_color_manual(values = species_colors, name = "Species") +
    labs(x = "Mean height of Y's gene trees",
         y = "False positive rate",
         title = "(b)") +
    base_theme

  g <- pa + pb + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0))

  if (save) {
    ggsave("figures/fig5_sigma_correlation.png", g, width = 12, height = 6, dpi = 300)
    cat("Saved figures/fig5_sigma_correlation.png\n")
  }
  g
}

plot_tree_height()