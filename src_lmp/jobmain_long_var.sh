#!/bin/bash

#SBATCH -A birthright
#SBATCH -p high_mem_cd
#SBATCH --time=py_tottime:30:00
#SBATCH --nodes=py_nnodes
#SBATCH --ntasks-per-node=py_ncores
#SBATCH --mem=24G
#SBATCH -J py_jobname
#SBATCH -o out.%J
#SBATCH -e err.%J

module rm PE-intel
module load PE-gnu

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

mpirun -np py_nptot ./lmp_mpi -in in.nvt -e screen
wait
mpirun -np py_nptot ./lmp_mpi -in in.npt -e screen
wait
mpirun -np py_nptot ./lmp_mpi -in in.rdf -e screen
