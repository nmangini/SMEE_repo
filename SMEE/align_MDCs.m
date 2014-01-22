% This script makes the MDC injections have the same time duration as the
% shortest MDC and aligns the peak in each waveform. The PC durations are
% then matched to the MDC durations

% Load MDCs
load('MDC_Q-series');
QMDC = MDCinj;
load('resampledtime_Q-series');
Qt = t-t(1);
%load('MDC_HR-series');
%HRMDC = MDCinj;
%load('resampledtime_HR-series');
%HRt = t-t(1);
%load('MDC_RO3-series');
%RO3MDC = MDCinj;
%load('resampledtime_RO3-series');
%RO3t = t-t(1);

% Load PCs
load('RvectorsPC_Q-series')
QPCs = RvectorsPC;
%load('RvectorsPC_HR-series')
%HRPCs = RvectorsPC;
%load('RvectorsPC_RO3-series')
%RO3PCs = RvectorsPC;

% Find relevent MDC times
[Qmax,Qmaxindex] = max(abs(QMDC(:,1)));
Qpeaktime = Qt(Qmaxindex);
%[HRmax,HRmaxindex] = max(abs(HRMDC(:,1)));
%HRpeaktime = HRt(HRmaxindex);
%[RO3max,RO3maxindex] = max(abs(RO3MDC(:,1)));
%RO3peaktime = RO3t(RO3maxindex);
%peakstart_diff = RO3maxindex-1;
peakstart_diff = 12076; %!!!!!ONLY FOR EXAMPLE SCRIPT!!!!!

% Cut begining of MDC data as necessary
Qt2 = Qt(Qmaxindex-peakstart_diff:length(Qt));
for i = 1:13
    QMDC_cut(:,i) = QMDC(Qmaxindex-peakstart_diff:length(QMDC),i);
end
%HRt2 = HRt(HRmaxindex-peakstart_diff:length(HRt));
%for i = 1:15
%    HRMDC_cut(:,i) = HRMDC(HRmaxindex-peakstart_diff:length(HRMDC),i);
%end

% Cut end of MDC data as necessary
Qt3 = Qt2(1:19661);
for j = 1:13
    QMDC_final(:,j) = QMDC_cut(1:19661,j);
end
%HRt3 = HRt2(1:19661);
%for j = 1:15
%    HRMDC_final(:,j) = HRMDC_cut(1:19661,j);
%end
%RO3t2 = RO3t(1:19661);
%for j = 1:20
%    RO3MDC_final(:,j) = RO3MDC(1:19661,j);
%end

% Find relevent PC times
[QPCmax,QPCmaxindex] = max(abs(QPCs(:,1)));
%[HRPCmax,HRPCmaxindex] = max(abs(HRPCs(:,1)));
%[RO3PCmax,RO3PCmaxindex] = max(abs(RO3PCs(:,1)));
%peakstart_diffPC = RO3PCmaxindex-1;
peakstart_diffPC = 12076; %!!!!!ONLY FOR EXAMPLE SCRIPT!!!!!
% Cut begining of PC data as necessary
QPCs_cut = zeros(length(QMDC_cut),13);
for k = 1:13
    QPCs_cut(:,k) = QPCs(QPCmaxindex-peakstart_diffPC:length(QPCs),k);
end
%for k = 1:15
%    HRPCs_cut(:,k) = HRPCs(HRPCmaxindex-peakstart_diffPC:length(HRPCs),k);
%end

% Cut end of PC data to match MDC Duration
for l = 1:13
    %QPCs_final(:,l) = QPCs_cut(1:length(RO3t2),l);
    QPCs_final(:,l) = QPCs_cut(1:19661,l); %!!!!!ONLY FOR EXAMPLE SCRIPT!!!!!
end
%for l = 1:15
%    HRPCs_final(:,l) = HRPCs_cut(1:length(RO3t2),l);
%end
%for l = 1:20
%    RO3PCs_final(:,l) = RO3PCs(1:length(RO3t2),l);
%end

% Plots
%figure(1)
%subplot(3,1,1)
%    plot(Qt,QMDC(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('Q MDC')
%    grid on
%    yL = get(gca,'YLim');
%    line([Qpeaktime Qpeaktime],yL,'Color','r')
%    line([Qt(Qmaxindex-peakstart_diff) Qt(Qmaxindex-peakstart_diff)],yL,'Color','g')
%    xL = get(gca,'xlim');
%subplot(3,1,2)
%    plot(HRt,HRMDC(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('HR MDC')
%    grid on
%    yL = get(gca,'YLim');
%    line([HRpeaktime HRpeaktime],yL,'Color','r')
%    line([HRt(HRmaxindex-peakstart_diff) HRt(HRmaxindex-peakstart_diff)],yL,'Color','g')
%    xlim(xL)
%subplot(3,1,3)
%    plot(RO3t,RO3MDC(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('RO3 MDC')
%    grid on
%    yL = get(gca,'YLim');
%    line([RO3peaktime RO3peaktime],yL,'Color','r')
%    xlim(xL)
%    
%figure(2)
%subplot(3,1,1)
%    plot(Qt3-Qt3(1),QMDC_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('Q MDC')
%    grid on
%subplot(3,1,2)
%    plot(HRt3-HRt3(1),HRMDC_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('HR MDC')
%    grid on
%subplot(3,1,3)
%    plot(RO3t2,RO3MDC_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('RO3 MDC')
%    grid on
    
%figure(3)
%subplot(3,3,1)
%    plot(Qt3-Qt3(1),QPCs_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('Q PC 1')
%    grid on
%subplot(3,3,4)
%    plot(Qt3-Qt3(1),QPCs_final(:,2))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('Q PC 2')
%    grid on
%subplot(3,3,7)
%    plot(Qt3-Qt3(1),QPCs_final(:,3))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('Q PC 3')
%    grid on
%subplot(3,3,2)
%    plot(HRt3-HRt3(1),HRPCs_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('HR PC 1')
%    grid on
%subplot(3,3,5)
%    plot(HRt3-HRt3(1),HRPCs_final(:,2))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('HR PC 2')
%    grid on
%subplot(3,3,8)
%    plot(HRt3-HRt3(1),HRPCs_final(:,3))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('HR PC 3')
%    grid on
%subplot(3,3,3)
%    plot(RO3t2-RO3t2(1),RO3PCs_final(:,1))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('RO3 PC 1')
%    grid on
%subplot(3,3,6)
%    plot(RO3t2-RO3t2(1),RO3PCs_final(:,2))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('RO3 PC 2')
%    grid on
%subplot(3,3,9)
%    plot(RO3t2-RO3t2(1),RO3PCs_final(:,3))
%    xlabel('t (s)')
%    ylabel('Strain')
%    title('RO3 PC 3')
%    grid on
    
% Save
QMDC_savefile = 'final-MDC_Q-series';
MDC_final = QMDC_final;
save(QMDC_savefile,'MDC_final')
%HRMDC_savefile = 'final-MDC_HR-series';
%MDC_final = HRMDC_final;
%save(HRMDC_savefile,'MDC_final')
%RO3MDC_savefile = 'final-MDC_RO3-series';
%MDC_final = RO3MDC_final;
%save(RO3MDC_savefile,'MDC_final')
QPC_savefile = 'finalRvectorsPC_Q-series';
PCs_final = QPCs_final;
save(QPC_savefile,'PCs_final')
%HRPC_savefile = 'finalRvectorsPC_HR-series';
%PCs_final = HRPCs_final;
%save(HRPC_savefile,'PCs_final')
%RO3PC_savefile = 'finalRvectorsPC_RO3-series';
%PCs_final = RO3PCs_final;
%save(RO3PC_savefile,'PCs_final')
