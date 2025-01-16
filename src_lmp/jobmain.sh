#!/bin/bash


#SBATCH --job-name si40_0.05_1
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=36
#SBATCH --time=47:30:00
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
module load lammps

# Change directory and signal job start
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD

# Create output directories
mkdir -p outdir
mkdir -p trajfiles

# Run jobs
srun -n 216 lmp -in in.init -e screen
wait
srun -n 216 lmp -in in.nve -e screen
wait
srun -n 216 lmp -in in.nvt -e screen
wait
srun -n 216 lmp -in in.nvt_main -e screen
