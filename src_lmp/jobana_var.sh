#!/bin/bash

#SBATCH --job-name py_jobname
#SBATCH --nodes=py_nnodes
#SBATCH --cpus-per-task=py_ncores
#SBATCH --time=py_tottime:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

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
