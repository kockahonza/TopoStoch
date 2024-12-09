#!/bin/bash
#SBATCH -p BIOP
#SBATCH --time 24:00:00
#requesting one node
#SBATCH -N1
#requesting 12 cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

echo echo "->" $SLURM_CPUS_PER_TASK $SLURM_CPUS_ON_NODE

julia -t $SLURM_CPUS_PER_TASK scripts/gmca_gil1.jl
