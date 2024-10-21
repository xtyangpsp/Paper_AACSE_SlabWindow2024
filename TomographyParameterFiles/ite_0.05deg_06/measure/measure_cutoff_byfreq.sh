#!/bin/bash
if [ $# -lt 1 ]
        then
                echo 'USAGE: ./measure_cutoff_byfreq.sh inputresultfile'
                exit 1
fi

dfile=$1 #'ite_0.05deg_02_measure_result_craton_dt6cc0.75snr7.dat'
#qc: dtmax cc_min snr_min
f1qc='10 0.7 5'
f2qc='10 0.7 5'
f3qc='10 0.7 5'
f4qc='5 0.7 5'
f5qc='5 0.7 5'
f6qc='5 0.7 5'
f7qc='5 0.7 5'
f8qc='5 0.7 5'

ftag='f1'
echo $ftag $f1qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh  ${fnm_meas} $f1qc > ${dfile}_${ftag}

ftag='f2'
echo $ftag $f2qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f2qc > ${dfile}_${ftag}

ftag='f3'
echo $ftag $f3qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f3qc > ${dfile}_${ftag}

ftag='f4'
echo $ftag $f4qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f4qc > ${dfile}_${ftag}

ftag='f5'
echo $ftag $f5qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f5qc > ${dfile}_${ftag}

ftag='f6'
echo $ftag $f6qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f6qc > ${dfile}_${ftag}

ftag='f7'
echo $ftag $f7qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f7qc > ${dfile}_${ftag}

ftag='f8'
echo $ftag $f8qc
fnm_meas=".measure_temp_"${ftag}
awk '{ if ( $7 == "'${ftag}'" ) print }' ${dfile} > ${fnm_meas}
./measure_cutoff_flat.sh ${fnm_meas} $f8qc > ${dfile}_${ftag}


