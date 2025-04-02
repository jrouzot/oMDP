#!/bin/sh 

#SBATCH --mem-per-cpu=10gb 
#SBATCH --array=0-239
#SBATCH --time=01:10:00 
#SBATCH --output results/RepairDescent/small/%a.out 
#SBATCH --error errors/RepairDescent/small/%a.err 

./build/bin/global-constraint-transfert-scheduling --repairdescent 1 --time_limit 3600000 instances/generated-small/${SLURM_ARRAY_TASK_ID}.txt