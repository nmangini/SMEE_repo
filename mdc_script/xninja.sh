# Use this script tp create xml files containing mass and spin parameters of the
# given waveform
# Arguments:
#           --WF - The waveform to be used
#           --Q  - Mass ratio of the waveform

WF=$1
Q=$2

/home/nmangini/ninja-mdcs/lalapps_ninja \
 --datadir ./ \
 --pattern ./$WF.gwf \
 --min-mass-ratio $Q \
 --max-mass-ratio $Q \
 --freq-lo 10 \
 --outfile sim-$WF.xml
