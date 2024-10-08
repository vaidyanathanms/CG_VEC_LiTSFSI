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

echo "begin job.."
echo $PWD

export OMP_NUM_THREADS=py_ncores

mkdir -p anaoutdir

./ana.o py_anainp
