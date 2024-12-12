#!/bin/bash
#SBATCH -p BIOP
#SBATCH --time 36:00:00
#requesting one task of 20 cpus (on one node implied)
#SBATCH --cpus-per-task=20

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    SLURM_CPUS_PER_TASK=10
fi

echo echo "echo testing ->" $SLURM_CPUS_PER_TASK $SLURM_CPUS_ON_NODE

julia -t $SLURM_CPUS_PER_TASK -L ../../scripts/gmca_gil1.jl -e 'run_extreme(; steps=100000)'
