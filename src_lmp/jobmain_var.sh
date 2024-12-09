#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH --time=py_tottime:30:00
#SBATCH --nodes=py_nnodes
#SBATCH --ntasks-per-node=py_ncores
#SBATCH --mem=2G
#SBATCH -J py_jobname
#SBATCH -o out.%J
#SBATCH -e err.%J

module rm PE-intel
module load PE-gnu

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

mpirun -np py_nptot ./lmp_mpi -in in.init -e screen
wait
mpirun -np py_nptot ./lmp_mpi -in in.nve -e screen
wait
mpirun -np py_ntot ./lmp_mpi -in in.nvt -e screen
wait
mpirun -np py_ntot ./lmp_mpi -in in.nvt_main -e screen
wait
mpirun -np py_nptot ./lmp_mpi -in in.npt -e screen
