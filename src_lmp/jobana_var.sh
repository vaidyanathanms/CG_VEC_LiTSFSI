#!/bin/bash

#SBATCH --job-name py_jobname
#SBATCH --nodes=py_nnodes
#SBATCH --ntasks-per-node=py_ncores
#SBATCH --time=py_tottime:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

# Signal job start job
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD


# Run jobs
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
mv std.err_%j allresults
mv std.err_%j allresults

wait

rm anainp_var.txt

echo "Analysis completed ..."
