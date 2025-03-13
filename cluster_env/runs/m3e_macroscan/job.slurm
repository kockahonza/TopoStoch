#!/bin/bash
#SBATCH -p BIOP
#SBATCH --time 48:00:00
#requesting one task of ? cpus (on one node implied)
#SBATCH --cpus-per-task=128

echo echo "echo testing ->" $SLURM_CPUS_PER_TASK $SLURM_CPUS_ON_NODE

time julia -t 20 -L job.jl -e 'macroscans(; testrun=true)'
# time julia -t $SLURM_CPUS_PER_TASK -L job.jl -e 'macroscans()'
