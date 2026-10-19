# Usage:
#   Rscript run_task.R <task_id> [offset]
#   task_id : integer, usually $SLURM_ARRAY_TASK_ID
#   offset  : optional integer added to task_id. If omitted, the env var
#             COMBO_OFFSET is used (default 0). Needed because the grid has
#             2100 combos but SLURM caps array indices (often at 1000).
# Output:  sim_output/Lambda<>_X<>_Y<>_Beta<>_N<>.csv, one row per replicate
#          sim_output/matrices/<same name>.rds, VCVs for the first N_SAVE replicates
# Depends: pgls_functions.R (defines one_rep()), digest

library(digest)

source("pgls_functions.R")

args    <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2 && !is.na(suppressWarnings(as.integer(args[2])))) {
  offset <- as.integer(args[2])
} else {
  offset <- as.integer(Sys.getenv("COMBO_OFFSET", "0"))
  if (is.na(offset)) stop("COMBO_OFFSET env var is set but not a valid integer: '",
                           Sys.getenv("COMBO_OFFSET"), "'")
}

task_id <- as.integer(args[1]) + offset
Reps    <- 1000 # replicates per parameter combination
message("Running task_id = ", task_id, " (raw arg = ", args[1], ", offset = ", offset, ")")

dir.create("sim_output", showWarnings = FALSE)

# Lambda : speciation rate  (values: 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)
# X_arch : architecture for trait X (values: 1 = single gene, 2 = 2 genes, 5 = 5 genes, 10 = 10 genes, 100 = 100 genes)
# Y_arch : architecture for trait Y ("same" = uses the same gene trees as X)
# Beta   : true regression slope of Y on X
# Nspec  : number of tips (species)

grid <- expand.grid(
  Lambda = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
  X_arch = c(1, 2, 5, 10, 100),
  Y_arch = c("same", "1", "2", "5", "10", "100"),
  Beta   = c(0, 0.25, 0.5, 0.75, 1.0),
  Nspec  = c(10, 100),
  stringsAsFactors = FALSE
)

if (task_id > nrow(grid)) stop("task_id ", task_id, " exceeds grid size (", nrow(grid), ")")
if (task_id < 1) stop("task_id ", task_id, " is out of bounds (", nrow(grid), ")")
p <- grid[task_id, ]

# Build the output filename from the parameter values.
fmt <- function(x) format(signif(x, 6), scientific = FALSE, trim = TRUE)

outfile <- sprintf(
  "sim_output/Lambda%s_X%s_Y%s_Beta%s_N%s.csv",
  fmt(p$Lambda), fmt(p$X_arch),
  p$Y_arch, fmt(p$Beta), fmt(p$Nspec)
)
if (file.exists(outfile)) quit(save = "no", status = 0)

# Reproducible seed tied to the parameter combo
seed <- strtoi(substr(digest::digest(basename(outfile), algo = "md5"), 1, 7), 16L)

N_SAVE  <- 10   # replicates per combo whose trees/VCVs are saved
mat_dir <- "sim_output/matrices"
dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

results <- vector("list", Reps)
saved   <- list()
for (r in seq_len(Reps)) {
  set.seed(seed + r)   # per-replicate seed: any replicate can be regenerated alone
  keep <- r <= N_SAVE
  res  <- one_rep(Nspec = p$Nspec, Lambda = p$Lambda, X_arch = p$X_arch,
                  Y_arch = p$Y_arch, Beta = p$Beta, return_matrices = keep)
  if (keep) { saved[[r]] <- res$vcvs; res <- res$results }
  results[[r]] <- res
}
# write the .rds before the .csv, so a finished .csv means both exist
saveRDS(saved, file.path(mat_dir, sub("\\.csv$", ".rds", basename(outfile))), compress = "xz")

out <- do.call(rbind, results)
out$rep      <- seq_len(Reps)
out$seed     <- seed
out$rep_seed <- seed + seq_len(Reps)   # set.seed(rep_seed) regenerates that replicate
out$Nspec    <- p$Nspec
out$Lambda   <- p$Lambda
out$X_arch   <- p$X_arch
out$Y_arch   <- p$Y_arch
out$Beta     <- p$Beta

# write to a temp name, then rename, so a job killed mid-write
# never leaves a partial CSV that looks finished
tmp <- paste0(outfile, ".tmp")
write.csv(out, tmp, row.names = FALSE)
file.rename(tmp, outfile)