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

targetSNR=50
Fmin=10
Fmax=100
NumPCA=8

Ndet=1 # DON"T CHANGE THIS
doplots=0
dodistance=1

# Make directory if it doesn't already exist
mkdir -p ${PWD}/${signal}${waveformN}-${model}-seed${jobnum}

# Your matlab path (matlab willl now know about all the codes in this directory)
export MATLABPATH=${MATLABPATH}:"${SMEEBBH_PREFIX}/SMEE"


echo matlab -nosplash -nojvm -nodisplay -r "SMEE_BBH('${signal}${waveformN}-${model}-seed${jobnum}','aligo','${model}','${signal}',${waveformN},${Fmin},${Fmax},${jobnum},${Ndet},${NumPCA},${doplots},'SNR',${targetSNR},${dodistance})"
matlab -nosplash -nojvm -nodisplay -r "SMEE_BBH('${signal}${waveformN}-${model}-seed${jobnum}','aligo','${model}','${signal}',${waveformN},${Fmin},${Fmax},${jobnum},${Ndet},${NumPCA},${doplots},'SNR',${targetSNR},${dodistance})"
