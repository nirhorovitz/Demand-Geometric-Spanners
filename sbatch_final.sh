#!/bin/bash
#SBATCH --partition=course
#SBATCH --qos=course
#SBATCH --job-name=final_experiment
#SBATCH --output=final_experiment_%J.out
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:rtx_3090:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# -u is used to watch live progress and prevent buffered printing 
python -u final_experiment.py --n 4000
python -u final_experiment.py --n 4000 --force-dgf
python -u final_experiment.py --n 8000 
python -u final_experiment.py --n 8000 --force-dgf
