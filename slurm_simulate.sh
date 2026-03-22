#!/bin/bash -l
#SBATCH --job-name=antimode_sim
#SBATCH --output=logs/simulate_%j.out
#SBATCH --error=logs/simulate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=6-00:00:00
#SBATCH --partition=longrun
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

echo "Starting full simulation on $(hostname) with $SLURM_CPUS_PER_TASK CPUs"
echo "Start: $(date)"

python src/simulate.py --full

echo "Simulation done: $(date)"
echo "Generating LaTeX tables..."

python src/csv_to_tex.py

echo "All done: $(date)"
