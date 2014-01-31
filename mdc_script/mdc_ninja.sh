WF=$1
Q=$2
/home/nmangini/ninja-mdcs/lalapps_mdc_ninja \
 --verbose \
 --injection-type NR \
 --injection-file GHLTV-GATech-$WF-900000005-105.xml \
 --gps-start-time 900000005 \
 --gps-end-time 900000110 \
 --ifo H1 \
 --sample-rate 16384 \
 --write-mdc-log \
 --frame-type Q-series_q$Q \
 --set-name Q-series \
 --mdc-log $WF.log \
 --write-frame \
 --freq-low-cutoff 10 \
 --snr-low 0 \
 --snr-high 1e6 \
 --out-xml-file out-$WF.xml \
 --fr-out-dir ./MDC/frames/ \
 --double-precision 
