########################################################################
# plot_helpers.R
# Shared setup, sourced by every fig*.R script. Keeps colors, theme,
# and the "proportion of reps with p < 0.05" summary logic consistent
# across figures.
########################################################################

library(data.table)
library(dplyr)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)

# --- load results (written by combine_results.R) ---
# results_combined/ holds one file per Nspec x X_arch (Nspec<n>_X<x>.csv.gz).
# Pass x_arch and/or nspec to read only those files; leave NULL to read all.
open_results <- function(x_arch = NULL, nspec = NULL, path = "results_combined") {
  files <- list.files(path, pattern = "^Nspec.*\\.csv\\.gz$", full.names = TRUE)
  keep  <- rep(TRUE, length(files))
  if (!is.null(x_arch)) keep <- keep & grepl(paste0("_X(", paste(x_arch, collapse = "|"), ")\\.csv\\.gz$"), files)
  if (!is.null(nspec))  keep <- keep & grepl(paste0("^Nspec(", paste(nspec, collapse = "|"), ")_"), basename(files))
  if (!any(keep)) stop("No results files in ", path, " for x_arch = ", toString(x_arch), ", nspec = ", toString(nspec))
  # Y_arch as text in every file (see combine_results.R)
  rbindlist(lapply(files[keep], fread, colClasses = list(character = "Y_arch")))
}

# --- fitting scenarios: display name -> p-value column ---
# Defined once here so every figure uses the same names, and the names
# always match the color/linetype vectors below (matched by exact name).
scenario_p_cols <- c(
  "Y's true tree(s)"                    = "y_p",
  "Species tree"                        = "species_p",
  "X's true tree(s)"                    = "x_p",
  "Random architecture"                 = "random_p",
  "500 gene trees, true species tree"   = "random500_p",
  "500 gene trees, random species tree" = "unrelated_p"
)
# --- consistent scenario colors/linetypes across all figures ---
scenario_colors <- setNames(
  c("#000000",   # black
    "#D55E00",   # vermillion
    "#0072B2",   # blue
    "#999999",   # grey
    "#009E73",   # bluish green
    "#E69F00"),  # orange
  names(scenario_p_cols))
scenario_linetypes <- setNames(
  c("solid", "solid", "solid", "dashed", "dotted", "solid"),
  names(scenario_p_cols))

# --- Y_arch display order: matched first, then increasing mismatch ---
y_arch_levels <- c("same", "1", "2", "5", "10", "100")
y_arch_labels <- c("Same as X", "1 gene", "2 genes", "5 genes", "10 genes", "100 genes")

order_y_arch <- function(df) {
  df$Y_arch <- factor(df$Y_arch, levels = y_arch_levels, labels = y_arch_labels)
  df
}

# --- collapse per-replicate p-values into rejection rate (FPR or TPR) ---
# `p_cols` = named vector: display name -> column name, e.g.
# c("Species tree" = "species_p", "X's true tree(s)" = "x_p", ...)
# n = replicates in the group; n_valid = replicates with a usable p-value.

summarize_rejection_rate <- function(df, group_vars, p_cols, alpha = 0.05) {
  out_list <- lapply(names(p_cols), function(scenario_name) {
    col <- p_cols[[scenario_name]]
    df %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(
        n = n(),
        n_valid = sum(is.finite(.data[[col]])),
        rejection_rate = mean(.data[[col]] < alpha, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(scenario = scenario_name)
  })
  bind_rows(out_list)
}

base_theme <- theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 17)
  )
