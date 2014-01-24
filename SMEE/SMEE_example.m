% This script serves as a simple example on how to run SMEE end to end using 
% the Q-series BBH merger waveforms. First, the principle components need to
% be created of the desired waveform catalogue along with some post
% prossessing of the signals and injections (assumed to be already
% created). Then SMEE can be run with the principle components using any
% injections dired (MDC frames). After SMEE runs plots can be made of the
% output files to see how well the code performed

% NOTE: SMEE is hard coded to run with a specific directory setup:
    % ~/SMEE_BBH/
    % ~/SMEE_BBH/SMEE/
    % ~/SMEE_BBH/SMEE/Results
    % ~/SMEE_BBH/PCA/
    % ~/SMEE_BBH/PCA/Q-series
    % ~/SMEE_BBH/Waveforms/Q-series
   
repo_root='/data/jclark/SMEEBBH/SMEE_repo';

addpath(genpath(repo_root))
addpath('/home/jclark/opt/xpipeline/share/xpipeline/matlab')

% Creating the principle components:
    % See these scripts's comments for explanation of all steps
    disp('creating principle components...')
    cd([repo_root,'/PCA/Q-series/Q_allmodes'])
    load_frame('Q')
    eval('Q_allmodes_loadmatrix')
    eval('Q_allmodes_pca')
%    Q_allmodes_match % optional
    disp('resampling waveforms to data rate (usually 16384 Hz)')
    eval('resampleQ')
    
% Pre-processing (waveform alignment, truncation, ...):
    cd([repo_root,'SMEE'])
    disp('Loading PCs and Injections')
    eval('load_frame')
    copyfile('finalRvectorsPC_Q-series.mat','~/SMEE_repo/SMEE/')
    eval('align_MDCs')
    copyfile('final-MDC_Q-series.mat','~/SMEE_repo/SMEE/')

% Running SMEE:
    disp('executing SMEE analysis...')
    SMEE_BBH('EXAMPLE-Q1-Q-seed13-snr10-pcs8','aligo','Q','Q',1,10,13,1,8,0,'SNR',10);
    
    % NOTE: This script can take anywhere from 15 minutes to a couple hours
    % to finish running depending on the complexity of the signals and the
    % number of PCs used.

% Plotting Results
    % These scripts are still under development to be more general/robust
    % so some of these scripts may not work
    %cd /data/nmangini/SMEE_BBH/SMEE/Results/
    %collect_results.m % plots Bayes factors
    %reconstruction.m % reconstructs signal from beta factors and principle 
                     % components
    %plot_params.m % plots beta values of nested sampling loop
