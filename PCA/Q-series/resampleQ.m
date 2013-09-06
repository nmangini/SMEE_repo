% This script resamples the vectorsPC waveforms to 16384 to match the
% sampling rate of the MDC injection. Then the MDC injection is modified to
% match the length of the PC vectors

tic
% Load data
load('vectors_Q-series.mat')
load('time_Q-series.mat')
load('injection_Q-series.mat')
load('Tinjection_Q-series.mat')

% Variables and matrix set up
fs = 16384;                 % New sampling frequency
a = zeros(length(QT),13);
for i = 1:13
    a(:,i) = QT;            % Time matrix to be added to the PC matrix
end
WF = vectorsPC;
WF(:,:,2) = a;              % 3-D pc matrix:(PCs) x (waveform number) x (time) 

% Resample PC vectors
for i = 1:13
    for cnt=1:length(WF)
        %new time vector with equal spacing equal to 1/16384 seconds
        k=1:length(WF(:,i,2));
        if length(unique(WF(:,i,2))) < length(k)
            t=min(WF(:,i,2)):1/fs:1;
            k=find(WF(:,i,2)<=1);
        else
            t=min(WF(:,i,2)):1/fs:WF(end,i,2);
        end
        RvectorsPC(:,i)=interp1(WF(k,i,2),WF(k,i,1),t); %RvectorsPC now has a fixed sample rate
    end
end
t = t';


% Find times to cut MDC injection (using first PC)
[peak, peakindex] = max(abs(RvectorsPC(:,1)));
[MDCpeak, MDCpeakindex] = max(abs(injections(:,1)));
startindex = peakindex - 1;
endindex = length(RvectorsPC)-peakindex;

% Adjust the MDC injection to match the length of the PC vectors
MDCinj = zeros(length(RvectorsPC),13);
for i = 1:size(injections,2)
    MDCinj(:,i) = injections(MDCpeakindex-startindex:MDCpeakindex+endindex,i);
end

% Plotting
figure(1)
plot(t,RvectorsPC(:,1))
grid on
yL = get(gca,'YLim');
line([t(peakindex) t(peakindex)],yL,'Color','r')
line([t(1) t(1)], yL,'Color','g')
line([t(length(RvectorsPC)) t(length(RvectorsPC))],yL,'Color','y')
title('1st PC, resmapled')
xlabel('t (s)')
ylabel('Strain')

figure(2)
plot(ts,injections(:,1))
grid on
yL = get(gca,'YLim');
line([ts(MDCpeakindex) ts(MDCpeakindex)], yL,'Color','r')
line([ts(MDCpeakindex)-ts(startindex) ts(MDCpeakindex)-ts(startindex)], yL,'Color','g')
line([ts(MDCpeakindex)+ts(endindex) ts(MDCpeakindex)+ts(endindex)],yL,'Color','y')
title('Originial MDC Injection')
xlabel('t (s)')
ylabel('Strain')

figure(3)
plot(t,MDCinj(:,1))
title('Adjusted MDC Injection')
grid on
yL = get(gca,'YLim');
line([t(peakindex) t(peakindex)],yL,'Color','r')
line([t(1) t(1)], yL,'Color','g')
line([t(length(RvectorsPC)) t(length(RvectorsPC))],yL,'Color','y')
xlabel('t (s)')
ylabel('Strain')

% Save files
PC_savefile = 'RvectorsPC_Q-series';
save(PC_savefile,'RvectorsPC')
MDC_savefile = 'MDC_Q-series';
save(MDC_savefile,'MDCinj')
time_savefile = 'resampledtime_Q-series';
save(time_savefile,'t')
toc