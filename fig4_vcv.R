########################################################################
# fig_vcv.R
# Heatmaps of the saved variance-covariance matrices for one replicate:
# species tree, X's and Y's gene trees, random architecture, 500 gene
# trees from the species tree, and 500 gene trees from an unrelated
# species tree.
#
# Matrices are shown as CORRELATIONS (cov2cor), because that is what
# fit_pgls() uses; the overall scale of a VCV does not affect the fit.
#
# Species are ordered by clustering the SPECIES-TREE matrix, and the same
# order is used in every panel, so the species tree appears as nested
# blocks and gene-tree discordance shows up as departures from them.
#
# Reads sim_output/matrices/<combo>.rds (written by run_task.R for the
# first N_SAVE replicates) and the matching results CSV (for the number
# of genes drawn in the random architecture).
########################################################################

library(ggplot2)
source("plot_helpers.R")

fmt <- function(x) format(signif(x, 6), scientific = FALSE, trim = TRUE)  # same as run_task.R

# panel titles, in display order (names = list names in the .rds)
vcv_titles <- c(
  species   = "Species tree",
  X         = "X's true tree(s)",
  Y         = "Y's true tree(s)",
  random    = "Random architecture",
  random500 = "500 gene trees, true species tree",
  unrelated = "500 gene trees, random species tree"
)

# Load one replicate's VCVs and return them in long format for ggplot.
vcv_long <- function(Lambda, X_arch, Y_arch, Beta, Nspec, rep = 1, dir = "sim_output") {
  stem <- sprintf("Lambda%s_X%s_Y%s_Beta%s_N%s",
                  fmt(Lambda), fmt(X_arch), Y_arch, fmt(Beta), fmt(Nspec))
  rds <- file.path(dir, "matrices", paste0(stem, ".rds"))
  if (!file.exists(rds)) stop("No saved matrices: ", rds)
  saved <- readRDS(rds)
  if (rep > length(saved)) stop("Only ", length(saved), " replicates saved in ", rds)
  vcvs <- saved[[rep]]

  # number of genes in the random architecture, from the results CSV
  csv <- file.path(dir, paste0(stem, ".csv"))
  titles <- vcv_titles
  if (file.exists(csv)) {
    res <- read.csv(csv)
    k <- res$random_arch[res$rep == rep]
    if (length(k) == 1) titles["random"] <- sprintf("Random architecture (%s gene%s)", k, if (k == 1) "" else "s")
  }

  # species order: hierarchical clustering of the species-tree matrix
  C_sp <- cov2cor(vcvs$species)
  ord  <- rownames(C_sp)[hclust(as.dist(1 - C_sp), method = "average")$order]

  do.call(rbind, lapply(names(vcv_titles), function(nm) {
    C <- cov2cor(vcvs[[nm]])[ord, ord]
    data.frame(panel = titles[[nm]],
               row = factor(rep(ord, times = length(ord)), levels = rev(ord)),
               col = factor(rep(ord, each = length(ord)), levels = ord),
               cor = as.vector(C))
  })) |>
    transform(panel = factor(panel, levels = unname(titles)))
}

plot_vcv <- function(Lambda, X_arch, Y_arch, Beta, Nspec, rep = 1, save = TRUE) {
  d <- vcv_long(Lambda, X_arch, Y_arch, Beta, Nspec, rep)

  g <- ggplot(d, aes(x = col, y = row, fill = cor)) +
    geom_raster() +
    scale_fill_viridis_c(limits = c(0, 1), name = "Correlation") +
    facet_wrap(~ panel, nrow = 2, labeller = label_wrap_gen(width = 26)) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = sprintf("Lambda = %s, X = %s, Y = %s, %s species (replicate %s)",
                         fmt(Lambda), X_arch, Y_arch, Nspec, rep)) +
    base_theme +
    theme(axis.text = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right")

  if (save) {
    fname <- sprintf("figures/fig4_vcv_Lambda%s_X%s_Y%s_Beta%s_N%s_rep%s.png",
                     fmt(Lambda), X_arch, Y_arch, fmt(Beta), Nspec, rep)
    ggsave(fname, g, width = 12, height = 8, dpi = 300)
    cat("Saved", fname, "\n")
  }
  g
}

# Example: how discordance grows with Lambda for one setting.
# Beta doesn't affect the matrices, so any saved Beta works.
for (L in c(1e-4, 1, 100)) plot_vcv(Lambda = L, X_arch = 5, Y_arch = "5", Beta = 0, Nspec = 10)