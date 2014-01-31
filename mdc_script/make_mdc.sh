#The purpose of this script is to generate MDC frames for all waveforms in a
#given directory

INPUTDIR=$1

# Get directory names
DIRS = `ls -l INPUTDIR | egrep '^d' | awk '{print $9}'`

#for DIR in $DIRS
#do
#    echo ${DIR}
#done
for DIR in $DIRS
do
    lalapps_inspinj \
        --seed 42 \
        --real8-ninja2 \
        --f-lower 10 \
        --waveform GATech \
        --gps-start-time 900000005 \
        --gps-end-time 900004101 \
        --time-step 16 \
        --time-interval 1 \
        --l-distr fixed \
        --d-distr uniform \
        --i-distr fixed \
        --fixed-inc 0.0 \
        --longitude 0.0 \
        --latitude 0.0 \
        --polarization 0.0 \
        --min-distance 1000000 \
        --max-distance 1000000 \
        --m-distr nrwaves \
        --min-mass1 0.5 \
        --max-mass1 0.5 \
        --min-mass2 0.5 \
        --max-mass2 0.5 \
        --min-mtotal 250 \
        --max-mtotal 250 \
        --nr-file sim-${DIR}.xml \
        --disable-spin \
        --output GHLTV-GATech-${DIR}-900000005-4096.xml
done
#lalapps_mdc_ninja \
#    --verbose \
#    --debug-level 1\
#    --injection-type NR \
#    --injection-file $1 \
#    --gps-start-time 900000005 \
#    --gps-end-time 900004101 \
#    --ifo H1 \
#    --sample-rate 16384 \
#    --write-mdc-log \
#    --frame-type NINETY_HIGHER_Gpc \
#    --set-name MAYA_q4 \
#    --mdc-log NR_GATech.log \
#    --write-frame \
#    --freq-low-cutoff 10 \
#    --snr-low 0 \
#    --snr-high 1e6 \
#    --out-xml-file out.xml \
#    --fr-out-dir ./MDC/frames \
    --double-precision 
