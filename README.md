# PGLS under gene-tree discordance

Simulation code for **The phylogenetic regression conundrum revisited: The
role of genetic architecture** (Megan L. Smith, Richard Adams, & Fábio K. Mendes).

We test how phylogenetic generalized least squares (PGLS) performs when traits evolve on gene trees that differ from the species tree, evaluating hte impacts of different genetic architectures for the two traits of interest. Traits are controlled by one or more loci whose gene trees are simulated under the multispecies coalescent. PGLS is then fitted using the species tree, the traits' true gene trees, and several misspecified alternatives. We report false positive rates (β = 0) and power (β ≠ 0). The design follows Adams et al., with the changes described in the manuscript.

---

## Repository contents

| File | Purpose |
|---|---|
| `pgls_functions.R` | Simulation and model fitting. Main entry point: `one_rep()` |
| `run_task.R` | Runs 1,000 replicates of one parameter combination (one SLURM array task) |
| `submit_block.sh` | Submits a range of combinations as SLURM array jobs |
| `combine_results.R` | Checks completeness and combines per-combination output |
| `plot_helpers.R` | Shared data loading, colors and theme for figures |
| `fig1_fpr_grid.R` | False positive rate vs λ |
| `fig2_tpr_grid.R` | True positive rate (power) vs λ |
| `fig3_rf_qc.R` | Gene-tree/species-tree discordance (Robinson–Foulds) vs λ |
| `fig_vcv.R` | Heatmaps of the covariance matrices used for fitting |
| `results/` | Summary tables (rejection rates per combination and method) |
| `figures/` | Figures from the paper |
| `versions.txt` | Exact R package versions used |

---

## Requirements

- R [4.4.x]
- R packages: ape, castor, phangorn, MASS, nlme, TreeSim, digest, data.table, dplyr, ggplot2. Exact versions are in `versions.txt`.

---

## Reproducing the results

Run all commands from the repository root; scripts use relative paths.

**1. Simulations and Inference** (2,100 parameter combinations × 1,000 replicates):

```bash
./submit_block.sh 1 2100
```

Each task writes `sim_output/<combination>.csv`, plus the covariance matrices for its first replicates in `sim_output/matrices/`. A task skips combinations that are already finished, so rerunning the same command fills in any that failed. To run one combination locally:

```bash
Rscript run_task.R 1
```

**2. Combine and check:**

```bash
Rscript combine_results.R
```

This writes `results_combined/` and lists any missing combinations in `missing_combos.csv`.

**3. Figures:**

```bash
Rscript fig1_fpr_grid.R
Rscript fig2_tpr_grid.R
Rscript fig3_rf_qc.R
Rscript fig_vcv.R
```

### Reproducibility

Each combination's seed is derived from its parameter values (an MD5 hash of its output filename), and replicate *r* uses `seed + r`. Any single replicate can be regenerated exactly:

```r
source("pgls_functions.R")
set.seed(rep_seed)   # from the rep_seed column of the results
one_rep(Nspec, Lambda, X_arch, Y_arch, Beta, return_matrices = TRUE)
```

---

## Parameters

| Parameter | Values | Meaning |
|---|---|---|
| `Lambda` | 10⁻⁴, 10⁻³, …, 10² | Speciation rate (extinction μ = λ/2). Sets tree depth in coalescent units, so higher λ means more incomplete lineage sorting |
| `X_arch` | 1, 2, 5, 10, 100 | Number of loci (gene trees) underlying trait X |
| `Y_arch` | `same`, 1, 2, 5, 10, 100 | Loci underlying trait Y. `same` means Y uses X's gene trees |
| `Beta` | 0, 0.25, 0.5, 0.75, 1 | True slope of Y on X |
| `Nspec` | 10, 100 | Number of extant species |

## Simulation and fitting

For each replicate:

1. Simulate a birth–death species tree with `Nspec` extant tips (`TreeSim::sim.bd.taxa`).
2. Simulate gene trees under the multispecies coalescent (`castor::generate_gene_tree_msc`; population size 1, generation time 1).
3. Simulate X by Brownian motion on the average covariance matrix of its gene trees. Y = βX + ε, where ε evolves the same way on Y's gene trees.
4. Fit Y ~ X by PGLS (`nlme::gls`, Brownian correlation structure) using each of:

| Method | Covariance used for fitting |
|---|---|
| `species` | Species tree |
| `x` | X's true gene tree(s) |
| `y` | Y's true gene tree(s) |
| `random` | Gene trees from the species tree, number of loci drawn at random from {1, 2, 5, 10, 100} |
| `random500` | 500 gene trees from the true species tree |
| `unrelated` | 500 gene trees from an independent species tree with the same tip count and height |

Where a method uses several gene trees, their covariance matrices are averaged.

## Output columns

One row per replicate:

| Column(s) | Description |
|---|---|
| `Nspec`, `Lambda`, `X_arch`, `Y_arch`, `Beta` | Parameter values |
| `rep`, `seed`, `rep_seed` | Replicate number, combination seed, replicate seed |
| `<method>_beta`, `_se`, `_t`, `_p` | Slope estimate, standard error, t statistic, two-sided p-value |
| `<method>_sigma` | Residual SD estimated by `gls` (ML) |
| `RF_X`, `RF_Y`, `RF_Random`, `RF_Random500`, `RF_Unrelated` | Mean normalized Robinson–Foulds distance of each gene-tree set to the species tree |
| `species_height` | Species tree height |
| `X_height`, `Y_height` | Mean gene-tree height for X and Y |
| `random_arch` | Number of loci drawn for the `random` method |

`<method>` is one of `species`, `x`, `y`, `random`, `random500`, `unrelated`.

---



## Contact

Megan Smith — ms4438@msstate.edu