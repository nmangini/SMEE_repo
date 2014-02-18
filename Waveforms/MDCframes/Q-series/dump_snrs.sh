#!/bin/bash

if [ $# -ne 2 ] 
then
    echo "usage: sh dump_snrs.sh FRAMENAME XMLFILE"
    exit
fi

gwfname=${1}
xmlname=${2}
txtname=`echo ${2} | sed 's/xml/txt/g'`

# Dump out times, distances
lwtprint ${xmlname} \
    -t sim_inspiral -c geocent_end_time,h_end_time,l_end_time,\
    v_end_time,distance,eff_dist_h,eff_dist_l,eff_dist_v \
    > ${txtname}

# Compute SNRs
export MATLABPATH="${SMEEBBH_PREFIX}/mdc_script":$MATLABPATH
echo matlab -nojvm -nodisplay -r "dump_mdc_snrs('${gwfname}','H1:STRAIN','${txtname}',10,8192,'nominal')"
matlab -nojvm -nodisplay -r "dump_mdc_snrs('${gwfname}','H1:STRAIN','${txtname}',10,8192,'nominal')"

# Cleanup
rm ${txtname}
reset


