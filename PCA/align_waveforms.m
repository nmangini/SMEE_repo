% Load waveforms
load('data_Q-series')
load('recon_Q-series')
load('data_HR-series')
load('recon_HR-series')
load('data_RO3-series')
load('recon_RO3-series')
load('time_RO3-series')

% Find relevent PC times
for i = 1:13
    [Qmax(i),Qmaxindex(i)] = max(abs(QreconstructedfullWF(:,i)));
end
for i = 1:15
    [HRmax(i),HRmaxindex(i)] = max(abs(HRreconstructedfullWF(:,i)));
end
for i = 1:20
    [RO3max(i),RO3maxindex(i)] = max(abs(RO3reconstructedfullWF(:,i)));
end
peakstart_diff = RO3maxindex(1)-1;

% Cut begining of waveform data as necessary
for k = 1:13
    Q_cut(:,k) = QfullWF(Qmaxindex(1)-peakstart_diff:length(QreconstructedfullWF),k);
end
for k = 1:15
   HR_cut(:,k) = HRfullWF(HRmaxindex(1)-peakstart_diff:length(HRreconstructedfullWF),k);
end
%for k = 1:20
    %RO3_cut(:,k) = RO3fullWF(RO3maxindex(1)-peakstart_diff:length(RO3fullWF),k);
    RO3_cut = RO3reconstructedfullWF;
%end

% Cut end of PC data to match waveform duration
for l = 1:13
   Q_recon(:,l) = Q_cut(1:length(RO3T),l);
end
for l = 1:15
   HR_recon(:,l) = HR_cut(1:length(RO3T),l);
end
for l = 1:20
   RO3_recon(:,l) = RO3_cut(1:length(RO3T),l);
end


% Save
Q_savefile = 'recon_Q';
save(Q_savefile,'Q_recon')
HR_savefile = 'recon_HR';
save(HR_savefile,'HR_recon')
RO3_savefile = 'recon_RO3';
save(RO3_savefile,'RO3_recon')