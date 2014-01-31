WF=$1
lalapps_inspinj \
 --seed 76 \
 --f-lower 10 \
 --waveform GATech \
 --gps-start-time 900000005 \
 --gps-end-time 900000110 \
 --time-step 10 \
 --time-interval 0 \
 --l-distr fixed \
 --d-distr uniform \
 --i-distr fixed \
 --fixed-inc 0.0 \
 --longitude 0.0 \
 --latitude 0.0 \
 --polarization 0.0 \
 --min-distance 1000000 \
 --max-distance 1000000 \
 --nr-file /home/nmangini/GATech_Studies/MDCframe_generation/Q-series/sim-$WF.xml \
 --m-distr nrwaves \
 --min-mtotal 250 \
 --max-mtotal 250 \
 --min-mass1 0.5 \
 --max-mass1 0.5 \
 --min-mass2 0.5 \
 --max-mass2 0.5 \
 --disable-spin \
 --output GHLTV-GATech-$WF-900000005-105.xml

