#!/bin/bash
#SBATCH -J dld_AACSE        #job name to remember
#SBATCH -n 5	#number of CPU cores you request for the job
#SBATCH -N 1 	#number of computing nodes, jobs could run across nodes, 
#but not recommended if not have to.
#SBATCH -A xtyang  #queue to submit the job, our lab queue.
#SBATCH --mem-per-cpu 5000 	#requested memory per CPU
#SBATCH -t 7-00:00			#requested time day-hour:minute
#SBATCH -o %x.out  #path and name to save the output file.
#SBATCH -e %x.err 	#path to save the error file.

module purge			#clean up the modules
module load rcac		#reload rcac modules.
module use /depot/xtyang/etc/modules
module load conda-env/seisgo-py3.7.6

mpirun -n $SLURM_NTASKS python download_MPI.py  
 					#run line, change file name
