#!/bin/bash

if [ $# -ne 3 ] 
then
    echo "insufficient args"
    exit
fi

# The name of the waveform catalogue
model=${1}

# The name of the signal catalogue
signal=${2}

# The waveform number
waveformN=${3}

# A random number for unique naming
runid=${RANDOM}

# The location for the output files
outdir=${waveform}

# The name of the condor submission file for this run
subfile=cat-${model}_inj-${signal}_${runid}.sub
shellfile=cat-${model}_inj-${signal}_${runid}.sh

echo '
#!/bin/bash 
####################
# SMEE BBH
####################
# Run this interactively

' >> ${shellfile}


# The header for the sub file
echo '
####################
# SMEE BBH
####################

executable = SMEE_BBH.sh
universe   = vanilla 

output     = condor_logs/$(Process).out
error      = condor_logs/$(Process).err
log        = condor_logs/$(Process).log

getenv=True
' >> ${subfile}

if [ ! -d condor_logs ]
then
    echo "making condor log directory"
    mkdir condor_logs
fi

#for job_number in `seq 11 20`
for job_number in 10
do
    echo """
    arguments = ${signal} ${job_number} ${waveformN} ${model}
    queue
    """ >> ${subfile}

    echo """
    sh SMEE_BBH.sh ${signal} ${job_number} ${waveformN} ${model}
    """ > ${shellfile}

done

        
