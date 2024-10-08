#!/bin/bash

#SBATCH -A birthright
#SBATCH -p high_mem_cd
#SBATCH -t 08:30:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=2G
#SBATCH -J test_CG_VEC
#SBATCH -o out.%J
#SBATCH -e err.%J


module load PE-gnu

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

mpirun -np 32 ./lmp_mpi -in in.init -e screen
wait
mpirun -np 32 ./lmp_mpi -in in.nve -e screen
wait
mpirun -np 32 ./lmp_mpi -in in.nvt -e screen
wait
mpirun -np 32 ./lmp_mpi -in in.npt -e screen
wait
mpirun -np 32 ./lmp_mpi -in in.rg -e screen
wait
mpirun -np 32 ./lmp_mpi -in in.rdf -e screen
