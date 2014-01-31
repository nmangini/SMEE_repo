#! /bin/sh
# Make env setup script
# James Clark, james.clark@ligo.org

# Get the location of the git repository by finding the full path of this file
# and stripping off the name
SMEEBBH_PREFIX=`python -c "import os, sys; print os.path.realpath('${0}')" | sed 's|/setup.sh||g'`
FRGETVECT_PATH="/home/jclark/opt/xpipeline/share/xpipeline/matlab/"

echo "export SMEEBBH_PREFIX=${SMEEBBH_PREFIX}" > smee_env.sh
echo "export FRGETVECT_PATH=${FRGETVECT_PATH}" >> smee_env.sh


