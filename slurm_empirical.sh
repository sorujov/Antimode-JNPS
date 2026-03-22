#!/bin/bash -l
#SBATCH --job-name=empirical
#SBATCH --output=logs/empirical_%j.out
#SBATCH --error=logs/empirical_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --partition=shortrun
#SBATCH --nodelist=sn1

source /share/common/anaconda/etc/profile.d/conda.sh
conda activate r_analysis_env

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /share/castor/home/orujov/Antimode-JNPS
mkdir -p logs results

echo "Empirical illustration on $(hostname)"
echo "Start: $(date)"

python src/empirical_faithful.py

echo "Done: $(date)"
