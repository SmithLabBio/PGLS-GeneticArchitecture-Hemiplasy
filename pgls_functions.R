# Functions for simulating species trees, gene trees and traits, and fitting PGLS/PIC regressions. Main entry point: one_rep().

library(ape)
library(castor)
library(phangorn)
library(MASS)
library(nlme)
library(TreeSim)

# one_rep(): simulate one dataset and fit every method.
# Returns a 1-row data.frame of results. If return_matrices = TRUE, returns
# list(results = <that data.frame>, vcvs = <VCVs used>) instead.
one_rep <- function(Nspec, Lambda, X_arch, Y_arch, Beta, extinction="half", return_matrices=FALSE, pagel_lambda_x = 1, n_shared = NULL) {
  
  if (extinction == "half"){
    Mu <- Lambda / 2
  } else if (extinction == "none"){
    Mu <- 0
  } else{
    stop("extinction must be 'half' or 'none'")
  }

  # simulations
  SpeciesTree <- TreeSim::sim.bd.taxa(n = Nspec, numbsim = 1, lambda = Lambda, mu = Mu, complete = FALSE)
  SpeciesTree[[1]] <- drop_extinct(SpeciesTree[[1]], Nspec)

  GenesX <- gene_trees(X_arch, SpeciesTree)
  if (!is.null(n_shared)) {
    if (as.character(Y_arch) != as.character(X_arch)) stop("n_shared is only supported when X_arch and Y_arch are the same")
    GenesY <- overlap_gene_trees(GenesX, n_shared, SpeciesTree)
  }else if(Y_arch == "same"){
    GenesY <- GenesX
  }else{
    GenesY <- gene_trees(Y_arch, SpeciesTree)
  }
  GenesRandom <- gene_trees(NULL, SpeciesTree)
  GenesRandom500 <- gene_trees(500, SpeciesTree)
  
  # Actually random gene trees: 500 gene trees from an independent species tree with
  # the same number of tips and the same height, used only for fitting.
  # mrca = TRUE makes `age` the crown age (= tree height); the default treats it
  # as the stem age, giving a shorter tree.
  height <- function(tr) max(node.depth.edgelength(tr))
  UnrelatedTree <- TreeSim::sim.bd.taxa.age(n = Nspec, numbsim = 1, lambda = Lambda, mu = Mu,
                                            age = height(SpeciesTree[[1]]), mrca = TRUE)
  GenesUnrelated <- gene_trees(500, UnrelatedTree)
  # Both trees label tips t1..tN, but in different orders. Reorder the VCV by
  # name to the species tree's order, which the traits are in. Which unrelated
  # tip pairs with which species is arbitrary, which is what makes it unrelated.
  sp <- SpeciesTree[[1]]$tip.label
  if (!setequal(rownames(GenesUnrelated$vcv), sp)) stop("Unrelated tree tip labels differ from species tree")
  vcv_unrelated <- GenesUnrelated$vcv[sp, sp]
  # RF to the true species tree (gene_trees() measured it against UnrelatedTree)
  RF_unrelated <- mean(sapply(GenesUnrelated$trees, function(g)
    RF.dist(SpeciesTree[[1]], g$tree, normalize = TRUE)))

  # Gene Tree Heights
  X_height = mean(diag(GenesX$vcv))   # mean height of X's gene trees
  Y_height = mean(diag(GenesY$vcv))   # mean height of Y's gene trees

  # traits
  Traits <- traits(GenesX, GenesY, Beta, pagel_lambda_x)
  
  # PGLS 
  pgls_species <- fit_pgls(Traits$TraitY, Traits$TraitX, NULL, tree = SpeciesTree[[1]]) # on species tree
  if (X_arch == 1){
    pgls_x <- fit_pgls(Traits$TraitY, Traits$TraitX, NULL, tree=GenesX$trees[[1]]$tree)
  } else{
    pgls_x <- fit_pgls(Traits$TraitY, Traits$TraitX, GenesX$vcv, tree=NULL)
  }
  if (Y_arch ==  1 | (Y_arch == "same" & X_arch == 1)){
    pgls_y <- fit_pgls(Traits$TraitY, Traits$TraitX, NULL, tree=GenesY$trees[[1]]$tree)
  } else{
    pgls_y <- fit_pgls(Traits$TraitY, Traits$TraitX, GenesY$vcv, tree=NULL)
  }
  if (GenesRandom$architecture == 1){
    pgls_random <- fit_pgls(Traits$TraitY, Traits$TraitX, NULL, tree=GenesRandom$trees[[1]]$tree)
  }else{
    pgls_random <- fit_pgls(Traits$TraitY, Traits$TraitX, GenesRandom$vcv, tree=NULL)
  }
  pgls_random500 <- fit_pgls(Traits$TraitY, Traits$TraitX, GenesRandom500$vcv, tree=NULL)
  pgls_unrelated <- fit_pgls(Traits$TraitY, Traits$TraitX, vcv_unrelated, tree=NULL)

  results <- data.frame(RF_X= GenesX$rf_mean, RF_Y = GenesY$rf_mean, RF_Random = GenesRandom$rf_mean, RF_Random500 = GenesRandom500$rf_mean, RF_Unrelated = RF_unrelated,
                        species_height = height(SpeciesTree[[1]]),
                        species_beta = pgls_species$beta, species_se = pgls_species$se, species_t = pgls_species$t, species_p = pgls_species$p, species_sigma = pgls_species$sigma,
                        x_beta = pgls_x$beta, x_se = pgls_x$se, x_t = pgls_x$t, x_p = pgls_x$p, x_sigma = pgls_x$sigma,
                        y_beta = pgls_y$beta, y_se = pgls_y$se, y_t = pgls_y$t, y_p = pgls_y$p, y_sigma = pgls_y$sigma,
                        random_beta = pgls_random$beta, random_se = pgls_random$se, random_t = pgls_random$t, random_p = pgls_random$p, random_sigma = pgls_random$sigma,
                        random500_beta = pgls_random500$beta, random500_se = pgls_random500$se, random500_t = pgls_random500$t, random500_p = pgls_random500$p, random500_sigma = pgls_random500$sigma,
                        unrelated_beta = pgls_unrelated$beta, unrelated_se = pgls_unrelated$se, unrelated_t = pgls_unrelated$t, unrelated_p = pgls_unrelated$p,  unrelated_sigma = pgls_unrelated$sigma,
                        random_arch = GenesRandom$architecture, height_x = X_height, height_y = Y_height)
  if (!is.null(n_shared)) results$n_shared <- n_shared
  if (!return_matrices) return(results)
 
  # The VCV behind each fit (averaged over gene trees where there are several).
  # All are in species-tree tip order.
  vcvs <- list(
    species   = ape::vcv(SpeciesTree[[1]]),
    X         = GenesX$vcv,
    Y         = GenesY$vcv,
    random    = GenesRandom$vcv,
    random500 = GenesRandom500$vcv,
    unrelated = vcv_unrelated
  )
  list(results = results, vcvs = vcvs)
}

gene_trees <- function(architecture, SpeciesTree) {
  if (is.null(architecture)){
    architecture <- sample(c(1, 2, 5, 10, 100), size=1, replace=FALSE)
  }
  gtrees <- list()
  rfs <- c()
  for (g in seq_len(architecture)){
    gtree <- generate_gene_tree_msc(SpeciesTree[[1]], population_sizes = 1, generation_times = 1)
    if (!identical(gtree$tree$tip.label, SpeciesTree[[1]]$tip.label)) {
      stop("Gene tree ", g, " tip order differs from the species tree")
    }
    gtrees[[g]] <- gtree
    rf <- RF.dist(SpeciesTree[[1]], gtree$tree, normalize=TRUE)
    rfs <- c(rfs, rf)
  }
  vcv <- lapply(gtrees, function(x) ape::vcv(x$tree))
  vcv_avg <- Reduce('+', vcv) / length(vcv)
  return(list(vcv = vcv_avg , rf_mean = mean(rfs), trees = gtrees, architecture=architecture))
}

traits <- function(X, Y, Beta, pagel_lambda_x = 1){
  if (pagel_lambda_x != 1){
    TEMP_X <- pagel_lambda_x * X$vcv + (1 - pagel_lambda_x) * diag(diag(X$vcv))
  } else{
    TEMP_X <- X$vcv
  }
  TraitX <- mvrnorm(1, rep(0, nrow(TEMP_X)), TEMP_X)
  ExpectedTraitY <- TraitX * Beta
  ResidualTraitY <- mvrnorm(1, rep(0, nrow(Y$vcv)), Y$vcv)
  TraitY <- ExpectedTraitY + ResidualTraitY
  return(list(TraitX = TraitX, TraitY = TraitY))
}

fit_pgls <- function(y, x, V, tree = NULL) {
  d <- data.frame(y = y, x = x)
  
  if (!is.null(tree)) {
    corstr <- corBrownian(1, phy = tree)
  } else {
    cor_matrix <- cov2cor(V)
    lower_tri  <- cor_matrix[lower.tri(cor_matrix)]
    corstr <- corSymm(lower_tri, fixed = TRUE)
  }
  
  fit <- gls(y ~ x, data = d, correlation = corstr, method = "ML")
  
  cf  <- coef(fit)
  vc  <- vcov(fit)
  se  <- sqrt(diag(vc))
  tval <- cf / se
  df   <- fit$dims$N - fit$dims$p
  pval <- 2 * pt(-abs(tval), df = df)
  
  return(list(beta = unname(cf["x"]), se = unname(se["x"]),
       t = unname(tval["x"]), p = unname(pval["x"]), sigma=fit$sigma))
}

# TreeSim treats lineages that went extinct within 1e-6 time units of the
# present as still alive (fixed tolerance in its drop.extinct call). Trees get
# shorter as Lambda grows, so this happens more often at high Lambda. TreeSim
# conditions on exactly n living species, so keep the n tips that reach the
# present and drop the rest. Relabel t1..tN so labels match the unrelated tree.
drop_extinct <- function(tr, n) {
  extra <- Ntip(tr) - n
  if (extra > 0) {
    depth <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
    tr <- drop.tip(tr, order(depth)[seq_len(extra)])
    tr$tip.label <- paste0("t", seq_len(Ntip(tr)))
  }
  if (Ntip(tr) != n || !is.ultrametric(tr)) {
    stop("Species tree has ", Ntip(tr), " tips (wanted ", n, ") or is not ultrametric")
  }
  tr
}

# Y's gene trees for the partial-overlap runs: the first n_shared of X's
# gene trees plus (K - n_shared) newly simulated ones, K = number of X's
# trees. X's trees are exchangeable, so which ones are shared is arbitrary.
# Returns the same structure as gene_trees().
overlap_gene_trees <- function(GenesX, n_shared, SpeciesTree) {
  K <- length(GenesX$trees)
  if (n_shared < 0 || n_shared > K || n_shared != round(n_shared)) {
    stop("n_shared must be a whole number from 0 to ", K)
  }
  shared <- GenesX$trees[seq_len(n_shared)]
  new    <- if (n_shared < K) gene_trees(K - n_shared, SpeciesTree)$trees else list()
  trees  <- c(shared, new)
  vcv    <- lapply(trees, function(x) ape::vcv(x$tree))
  rfs    <- sapply(trees, function(x) RF.dist(SpeciesTree[[1]], x$tree, normalize = TRUE))
  list(vcv = Reduce('+', vcv) / length(vcv), rf_mean = mean(rfs),
       trees = trees, architecture = K)
}
