#!/bin/bash
# Submits a contiguous range of combo IDs [START_ID, END_ID] as multiple
# SLURM array jobs, chunked to stay under this cluster's MaxArraySize (1000).
# The offset for each chunk is written directly into a generated sbatch
# script as a literal argument -- NOT passed via --export -- since
# --export=COMBO_OFFSET=... was not propagating into the job environment
# on this cluster.
#
# Usage: ./submit_block.sh START_ID END_ID


START_ID=$1
END_ID=$2
CHUNK=999

if [ -z "$START_ID" ] || [ -z "$END_ID" ]; then
  echo "Usage: $0 START_ID END_ID"
  exit 1
fi

mkdir -p chunks_weak logs_weak sim_output_weak

current=$START_ID
while [ "$current" -le "$END_ID" ]; do
  chunk_end=$((current + CHUNK - 1))
  if [ "$chunk_end" -gt "$END_ID" ]; then
    chunk_end=$END_ID
  fi
  local_size=$((chunk_end - current + 1))
  offset=$((current - 1))   # local index 1 + offset = current (real combo id)

  script="chunks/weak_${current}_${chunk_end}.sbatch"
  cat > "$script" << INNER
#!/bin/bash
#SBATCH --job-name=weak_${current}_${chunk_end}
#SBATCH --array=1-${local_size}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10:00:00
#SBATCH --output=logs_weak/weak_%A_%a.out
#SBATCH --error=logs_weak/weak_%A_%a.err
#SBATCH --partition=smith

module load R/4.4
Rscript run_task_weak.R \$SLURM_ARRAY_TASK_ID ${offset}
INNER

  echo "Submitting combo IDs ${current}-${chunk_end} (local array 1-${local_size}, offset ${offset}) via ${script}"
  sbatch "$script"

  current=$((chunk_end + 1))
done