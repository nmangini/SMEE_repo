This README serves as a simple example on how to run SMEE end to end using 
the Q-series BBH merger waveforms. The goal of SMEE is to have a fast,
unmodeled characterization of a gravitational wave from a binary black hole 
source which will tell us if the system is spinning, precessing, etc. First, 
the principle components need to be created of the desired waveform catalogue 
along with some post processing of the signals and injections (assumed to be
already created). Then SMEE can be run with the principle components using any
injections dired (MDC frames). After SMEE runs plots can be made of the
output files to see how well the code performed

NOTE: SMEE is hard coded to run with a specific directory setup:
 ~/SMEE_BBH/
 ~/SMEE_BBH/SMEE/
 ~/SMEE_BBH/SMEE/Results
 ~/SMEE_BBH/PCA/
 ~/SMEE_BBH/PCA/Q-series
 ~/SMEE_BBH/PCA/HR-series
 ~/SMEE_BBH/PCA/RO3-series
 ~/SMEE_BBH/Waveforms/Q-series
 ~/SMEE_BBH/Waveforms/HR-series
 ~/SMEE_BBH/Waveforms/RO3-series

CHANGE PATH NAMES IN THESE LOCATION:
    ~/SMEE_repo/PCA/Q-series/Q_allmodes_loadmatrix.m 
        line 21
    ~/SMEE_repo/Utilities/load_frame.m 
        lines 14, 24, 35
    ~/SMEE_repo/SMEE/SMEE_BBH.m
        lines 124, 126, 660

AN EXAMPLE RUN
    >Run ~/SMEE_repo/SMEE/SMEE_example.m
    >All the functions live in ~/SMEE_repo/Utilities
    >All plots are commented out in this repository as to run this code in a
    terminal
