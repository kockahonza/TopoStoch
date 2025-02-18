#!/bin/bash
#SBATCH -p BIOP
#SBATCH --time 96:00:00
#requesting one task of ? cpus (on one node implied)
#SBATCH --cpus-per-task=128

echo echo "echo testing ->" $SLURM_CPUS_PER_TASK $SLURM_CPUS_ON_NODE

time julia -t $SLURM_CPUS_PER_TASK -L ../../scripts/m3_macroscan.jl -e 'macroscan(5)'
