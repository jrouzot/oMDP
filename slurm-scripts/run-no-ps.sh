#!/bin/sh 

#SBATCH --mem-per-cpu=10gb 
#SBATCH --array=0-19
#SBATCH --time=01:10:00 
#SBATCH --output results/Ours/no-ps/%a.out 
#SBATCH --error errors/Ours/no-ps/%a.err 

./build/bin/global-constraint-transfert-scheduling --symmetry 0 --time_limit 3600000 instances/4-buffers/${SLURM_ARRAY_TASK_ID}.txt