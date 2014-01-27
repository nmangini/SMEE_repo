% End-to-end SMEE example

% Setting up paths
repo_root='~/SMEE_repo';
addpath(genpath(repo_root))
%addpath('/home/jclark/opt/xpipeline/share/xpipeline/matlab')

% Creating the principle components:
    % See these scripts's comments for explanation of all steps
    disp('################################')
    disp('Creating principle components...')
    disp('################################')
    cd([repo_root,'/PCA/Q-series/Q_allmodes'])
    load_frame('Q')
    eval('Q_allmodes_loadmatrix')
    eval('Q_allmodes_pca')
    %Q_allmodes_match % optional
    disp('Resampling waveforms to data rate (usually 16384 Hz)')
    eval('resampleQ')
    
% Pre-processing (waveform alignment, truncation, ...):
    cd([repo_root,'SMEE'])
    disp('Loading PCs and Injections')
    eval('load_frame')
    copyfile('finalRvectorsPC_Q-series.mat','~/SMEE_repo/SMEE/')
    eval('align_MDCs')
    copyfile('final-MDC_Q-series.mat','~/SMEE_repo/SMEE/')

% Running SMEE:
    disp('##########################')
    disp('Executing SMEE analysis...')
    disp('##########################')
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
