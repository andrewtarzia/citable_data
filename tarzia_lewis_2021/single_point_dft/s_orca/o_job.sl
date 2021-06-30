#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=128
#SBATCH --time=4:0:00
#SBATCH --error="%x.e%j"
#SBATCH --output="%x.o%j"

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=e05-discov-kim
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the batch environment
module load epcc-job-env

# Usage of this script:
#sbatch -J jobname job-orca-SLURM.sh , where jobname is the name of your ORCA inputfile (jobname.inp).

# Jobname below is set automatically when submitting like this: sbatch -J jobname job-orca.sh
#Can alternatively be set manually below. job variable should be the name of the inputfile without extension (.inp)
job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})

export OMP_NUM_THREADS=1

#Setting OPENMPI paths here:
export PATH=/work/e05/e05/atarzia/openmpi-3.1.4/bin:$PATH
export LD_LIBRARY_PATH=/work/e05/e05/atarzia/openmpi-3.1.4/lib:$LD_LIBRARY_PATH

# Here giving the path to the ORCA binaries and giving communication protocol
#You can also load module here.
export orcadir=/work/e05/e05/atarzia/orca_4_2_1_linux_x86-64_openmpi314
export RSH_COMMAND="/usr/bin/ssh -x"
export PATH=$orcadir:$PATH
export LD_LIBRARY_PATH=$orcadir:$LD_LIBRARY_PATH


# Creating nodefile in scratch
echo $SLURM_NODELIST > $job.nodes


# Copy job and node info to beginning of outputfile
echo "Job execution start: $(date)" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Shared library path: $LD_LIBRARY_PATH" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  $SLURM_SUBMIT_DIR/$job.out
echo $SLURM_NODELIST >> $SLURM_SUBMIT_DIR/$job.out

#Start ORCA job. ORCA is started using full pathname (necessary for parallel execution). Output file is written directly to submit directory on frontnode.
$orcadir/orca $job.in >>  $SLURM_SUBMIT_DIR/$job.out

