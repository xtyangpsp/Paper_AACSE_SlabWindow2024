#!/bin/bash
#SBATCH -J dipping_slab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A xtyang
#SBATCH --mem-per-cpu 20G
#SBATCH -t 6:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err
#module load intel
#module load netcdf-fortran/4.5.3
for base in slab_100km_45degrees
do
	./run.all.sh ${base}
done
