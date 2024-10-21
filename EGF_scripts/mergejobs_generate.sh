#!/bin/bash
sourcelist='allsources.list'

submitdir='merge_submit'
mkdir $submitdir

submit_template='submit_merge_pairs_bysources_template.sh' 

for s in `cat $sourcelist`
do
	sed -e 's/TEMPLATESTATION/'$s'/g' ${submit_template} > $submitdir/${s}.submit 

done



