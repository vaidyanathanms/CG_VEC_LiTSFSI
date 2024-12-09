#!/bin/bash

#SBATCH --account=iontransport
#SBATCH --time=47:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -p shared
#SBATCH -J si40_0.05_1
#SBATCH -o out.%J
#SBATCH -e err.%J

module purge

module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel
module load lammps


cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

srun -n 1 lmp -in in.init -e screen
