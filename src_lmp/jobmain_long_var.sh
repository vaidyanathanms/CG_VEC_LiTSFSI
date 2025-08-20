#!/bin/bash

#SBATCH --job-name py_jobname
#SBATCH --nodes=py_nnodes
#SBATCH --ntasks-per-node=py_ncores
#SBATCH --time=py_tottime:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition py_partition

# Load LAMMPS module
module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel
module load lammps//062322-intel-mpich

# Change directory and signal job start
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD

# Create output directories
mkdir -p outdir
mkdir -p trajfiles

#srun --mpi=pmi2 lmp -in in.nvt -e screen
#wait
srun --mpi=pmi2 lmp -in in.nvt_main -e screen
