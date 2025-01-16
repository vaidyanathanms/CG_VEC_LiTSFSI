#!/bin/bash

#SBATCH --job-name py_jobname
#SBATCH --nodes=py_nnodes
#SBATCH --time=py_tottime:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

# Load LAMMPS module
module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel
module load lammps

# Change directory and signal job start
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD

# Create output directories
mkdir -p outdir
mkdir -p trajfiles

srun -n py_nptot ./lmp_mpi -in in.nvt -e screen
wait
srun -n py_nptot ./lmp_mpi -in in.nvt_main -e screen
wait
srun -n py_nptot ./lmp_mpi -in in.npt -e screen
wait
srun -n py_nptot ./lmp_mpi -in in.rdf -e screen
