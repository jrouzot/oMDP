#!/bin/sh 

#SBATCH --mem-per-cpu=10gb 
#SBATCH --array=0-19
#SBATCH --time=01:10:00 
#SBATCH --output results/Ours/baseline/%a.out 
#SBATCH --error errors/Ours/baseline/%a.err 

./build/bin/global-constraint-transfert-scheduling --global 0 --single_window 0 --dense_ranking 0 --symmetry 0 --time_limit 3600000 instances/4-buffers/${SLURM_ARRAY_TASK_ID}.txt