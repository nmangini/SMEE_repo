# This script creates the frame file and xml of the given waveforms
# Arguments:
#           --INPUTDIR - Directory containing the waveforms
#           --INPUTCAT - Name of waveform catalogue
INPUTDIR=$1
INPUTCAT=$2

DIRS=`ls -l $INPUTDIR | egrep '^d' | awk '{print $9}'`
FULLPATH=/home/nmangini/GATech_Studies/Waveforms/$INPUTCAT/

for DIR in $DIRS
do
    # fr_ninja creates a .gwf file pointing to the waveform files in the input
    # directory
    /home/nmangini/ninja-mdcs/lalapps_fr_ninja \
        --verbose \
        --format NINJA1 \
        --double-precision \
        --nr-meta-file $FULLPATH${DIR}/Q_${DIR}.bbh \
        --nr-data-dir $FULLPATH${DIR} \
        --output ${DIR}.gwf
done
