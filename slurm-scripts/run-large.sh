#!/bin/sh 

#SBATCH --mem-per-cpu=10gb 
#SBATCH --array=0-239
#SBATCH --time=01:10:00 
#SBATCH --output results/Ours/large/%a.out 
#SBATCH --error errors/Ours/large/%a.err 

./build/bin/global-constraint-transfert-scheduling --time_limit 3600000 instances/generated-large/${SLURM_ARRAY_TASK_ID}.txt