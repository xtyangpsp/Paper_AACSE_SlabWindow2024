#!/bin/bash
if [ $# -lt 1 ]
	then
		echo './run.ckb.syn.sh infile'
		exit 1
fi
infile=$1

bin_ckb_syn='../../Codes/InvCkb/inv.ckb/inv_ckb_synthetic'

# input list file name

list_Gd_rcd="inv_Gd_list_observed"

fnm_ckb=${infile}
list_Gd='inv_Gd_list_'${infile}'.ori'

# num_cmp should be 2 when both Vp and Vs are inverted
num_cmp=2
# see ./run.kernel.assm.sh in sim.kernel
ker_thres="5e-3"

# maximum absolute synthetic delay time
data_thres=15

# ----------------------------------------------------------------------
${bin_ckb_syn} << EOF
$num_cmp
$fnm_ckb
$list_Gd_rcd
$list_Gd
$ker_thres $data_thres
EOF
