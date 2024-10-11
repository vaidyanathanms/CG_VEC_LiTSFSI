#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH --time=py_tottime:30:00
#SBATCH --nodes=py_nnodes
#SBATCH --ntasks-per-node=py_ncores
#SBATCH --mem=24G
#SBATCH -J py_jobname
#SBATCH -o allresults/out.%J
#SBATCH -e allresults/err.%J

echo "begin job.."
echo $PWD

export OMP_NUM_THREADS=py_ncores

./ana.o py_anainp
wait

# Clean-up

mv rgavg_pyconfig allresults
mv rgall_pyconfig allresults
mv rdf_pyconfig allresults
mv log.pyconfig allresults
mv anainp_pycase.txt allresults
mv jobana_pycase.sh allresults

echo "Analysis completed ..."
