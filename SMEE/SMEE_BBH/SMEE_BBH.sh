#!/bin/bash

if [ $# -ne 4 ] 
then
    echo "insufficient args"
    exit
fi

signal=${1}
jobnum=${2}
waveformN=${3}
model=${4}

# Make directory if it doesn't already exist
mkdir -p /data/nmangini/SMEE_BBH/SMEE/Results/${signal}${waveformN}-${model}-seed${jobnum}-snr10-pcs8

# Your matlab path (matlab willl now know about all the codes in this directory)
export MATLABPATH=${MATLABPATH}:"/data/nmangini/SMEE_BBH/SMEE"

matlab -singleCompThread -nosplash -nojvm -nodisplay -r "SMEE_BBH('${signal}${waveformN}-${model}-seed${jobnum}-snr10-pcs8','aligo','${model}','${signal}',${waveformN},10,${jobnum},1,8,0,'SNR',10)"
