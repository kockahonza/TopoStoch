#!/bin/bash
#SBATCH -p BIOP
#SBATCH --time 36:00:00
#requesting one task of 20 cpus (on one node implied)
#SBATCH --cpus-per-task=20

echo echo "echo testing ->" $SLURM_CPUS_PER_TASK $SLURM_CPUS_ON_NODE

julia -t 20 -L ../../scripts/m3_bigscan2.jl -e 'bigrun2(4; testrun=true)'
