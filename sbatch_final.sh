#!/bin/bash
#SBATCH --partition=course
#SBATCH --qos=course
#SBATCH --job-name=final_experiment
#SBATCH --output=final_experiment%J.out
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:rtx_3090:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Load the cluster's Anaconda module
module load anaconda

# Activate your specific environment
# Ensure the environment is deactivated on the manager node before sbatch [cite: 52]
source activate spanner_env

pip install -r requirements.txt

# Run the grid search with unbuffered output
# -u is used to watch live progress and prevent buffered printing [cite: 801, 805]
python -u final_experiment.py
