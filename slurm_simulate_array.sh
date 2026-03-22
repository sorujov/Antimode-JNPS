#!/bin/bash -l
#SBATCH --job-name=antimode_array
#SBATCH --output=logs/simulate_array_%A_%a.out
#SBATCH --error=logs/simulate_array_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=6-00:00:00
#SBATCH --partition=longrun
#SBATCH --array=0-9               # 10 chunks × 50 reps = 500 total
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sorujov@ada.edu.az

# ── Activate conda environment ────────────────────────────────────────
source /share/common/anaconda/etc/profile.d/conda.sh
conda activate r_analysis_env

# Prevent numpy/scipy/OpenBLAS from spawning extra threads per worker.
# We parallelise at the process level (Pool), so each worker must use 1 thread.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ── Run from project root ─────────────────────────────────────────────
cd /share/castor/home/orujov/Antimode-JNPS
mkdir -p logs results

echo "Array job ${SLURM_ARRAY_JOB_ID}, chunk ${SLURM_ARRAY_TASK_ID}/9"
echo "Node: $(hostname) | CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

python src/simulate.py --full \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --n-chunks 10

echo "Chunk done: $(date)"
