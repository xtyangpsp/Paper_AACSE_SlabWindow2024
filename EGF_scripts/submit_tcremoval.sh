#!/bin/bash
#SBATCH -J tc_test                # Job name
#SBATCH -n 40                     # Number of cores
#SBATCH -A xtyang		 # Partition: shared OR xtyang
#SBATCH --mem-per-cpu 5000        # Memory request
#SBATCH -t 1-23:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o /depot/xtyang/data/projects/vsassard/%x.out        # Standard output
#SBATCH -e /depot/xtyang/data/projects/vsassard/%x.err        # Standard error

module purge
module load rcac
module use /depot/xtyang/etc/modules
module load conda-env/seisgo-py3.7.6

#source activate seisgo # ~/noisepy/bin/activate

mpirun -n $SLURM_NTASKS python tcefficiency_MPI.py

#to save storage remove uncorrected data that won't be useful anymore, switch to True only if confident in parameters set in tcefficiency_MPI.py
#cleanrawdata = False
#if cleanrawdata == True:
#  rm AACSE_region/Raw_uncorrected_test/*
#  print("Uncorrected data removed")
#del AACSE_region/Raw_uncorrected_test/*.h5
