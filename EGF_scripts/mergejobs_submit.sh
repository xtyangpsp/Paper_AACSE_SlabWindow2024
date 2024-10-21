#!/bin/bash

submitdir='merge_submit'

#sleep 7200

for j in `ls $submitdir/*.submit`
do
	echo $j	
	sbatch $j
	sleep 0.05
done
