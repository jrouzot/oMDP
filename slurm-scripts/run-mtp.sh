#!/bin/sh 

#SBATCH --mem-per-cpu=10gb 
#SBATCH --array=1-4
#SBATCH --time=01:10:00 
#SBATCH --output results/Ours/mtp/%a.out 
#SBATCH --error errors/Ours/mtp/%a.err 

./build/bin/global-constraint-transfert-scheduling --time_limit 3600000 instances/mtp/MTP01${SLURM_ARRAY_TASK_ID}.txt