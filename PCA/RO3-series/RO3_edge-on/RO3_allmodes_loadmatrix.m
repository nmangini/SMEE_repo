% This script creates a matrix of waveforms for the RO3-series catalogue,
% aligns the waveforms, zeros the data after the ringdown signal, shortens
% the signal duration, and combines the spherical harmonic modes

tic
% Set up variables and solution arrays
maxlength = 16493;
qvalue = [1.5; 1.5; 1.5; 1.5; 2; 2; 2; 2; 2; 2; 2; 2.5; 2.5; 2.5; 2.5; 2.5; 2.5; 4; 4; 4];
avalue = [0.4; 0.6; 0.6; 0.6; 0.4; 0.6; 0.6; 0.6; 0.6; 0.6; 0.6; 0.4; 0.4; 0.4; 0.6; 0.6; 0.6; 0.6; 0.6; 0.6];
othvalue = [60; 45; 60; 90; 60; 45; 60; 90; 135; 180; 270; 45; 60; 90; 45; 60; 90; 45; 60; 90];
qvaluestr = cellstr(num2str(qvalue,'%4.2f')); % String of q values
avaluestr = cellstr(num2str(avalue,'%01.2f')); % String of a values
othvaluestr = cellstr(num2str(othvalue,'%03.0f')); % String of o values
RO3rawWF = zeros(maxlength,length(qvalue),3,9); % rawWF = (length of waveform)X(q value)X(l value)X(m value)
theta = 90*pi/180;
phi = 0*pi/180;

% Read in each waveform
fprintf('Reading in waveforms...')
for wfcnt1 = 1:length(qvalue)
    for lcnt1 = 2:4;
        for mcnt1 = -lcnt1:lcnt1
            fname = strcat('/Users/Nick/LIGO/SMEE_repo/Waveforms/RO3-series/RO3_D10_q',qvaluestr{wfcnt1},...
                '_a',avaluestr{wfcnt1},'_oth.',othvaluestr{wfcnt1},'_M120/rStrain_FFI_l', num2str(lcnt1), '_m',...
                num2str(mcnt1) ,'_r75.00.ampphase.asc');
            importfile(fname);
            RO3rawWF(1:length(data(:,2)),wfcnt1,lcnt1-1,mcnt1+5) = data(:,2);
        end
    end
    if length(data(:,1)) > maxlength-1
        RO3rawT = data(:,1);
    end
end
fprintf('done\n')

% Set up matrices, solution arrays, and variables for waveform processing
x = size(RO3rawWF);
L = x(1);
s = x(2);
peak = zeros(s,1);
peaktime = zeros(1,s);

% Find the peak time in each (2,2) waveform
fprintf('Finding the peak time of the (2,2) mode...')
starttime = 0;
i = 0;
while starttime == 0
    i = i + 1;
    if floor(RO3rawT(i)) == 200 % finds when t = 200 ms
        starttime = i;
    end
end

for wfcnt2 = 1:s
    for tcnt1 = 1:L
        if RO3rawWF(tcnt1,wfcnt2,1,3)^2 > peak(wfcnt2)
            peak(wfcnt2) = RO3rawWF(tcnt1,wfcnt2,1,3).^2;
            peaktime(wfcnt2) = tcnt1;
        end
    end
end
starttime = starttime + (max(peaktime) - min(peaktime));
fprintf('done\n')

% Align waveforms along thier peak
fprintf('Aligning waveforms...')
shiftedWF = zeros(L,s,3,9);
for wfcnt3 = 1:s
    for lcnt2 = 2:4
        for mcnt2 = -lcnt2:lcnt2
            for tcnt2 = (max(peaktime)-peaktime(wfcnt3)+1):L
                shiftedWF(tcnt2,wfcnt3,lcnt2-1,mcnt2+5) = RO3rawWF(tcnt2-(max(peaktime) - peaktime(wfcnt3)),wfcnt3,lcnt2-1,mcnt2+5);
                if RO3rawT(tcnt2) > (max(peaktime) + 200)
                    shiftedWF(tcnt2,wfcnt3,lcnt2-1,mcnt2+5) = 0; % sets data after ringdown to 0
                end
            end
        end
    end
end
fprintf('done\n')

% Removes data prior to 200 ms
fprintf('Removing data prior to 200 ms...')
WF = zeros(length(shiftedWF)-starttime,s,x(3),x(4));
for i = starttime:length(shiftedWF)
    for wfcnt4 = 1:s
        for lcnt3 = 2:4
            for mcnt3 = -lcnt3:lcnt3
                WF(i-starttime+1,wfcnt4,lcnt3-1,mcnt3+5) = shiftedWF(i,wfcnt4,lcnt3-1,mcnt3+5);
            end
        end
    end
    RO3T(i-starttime+1) = RO3rawT(i);
end

% Convert from NR units to physical time units
RO3T = RO3T';
RO3T = RO3T*4.92549095e-6*250; %250 M_sun
fprintf('done\n')

% Combine all spherical harmonic modes
fprintf('Combining spherical harmonic modes...')
RO3fullWF = zeros(length(WF),s);
for wfcnt5 = 1:s
    for lcnt4 = 2:4
        for mcnt4 = -lcnt4:lcnt4
            RO3fullWF(:,wfcnt5) = RO3fullWF(:,wfcnt5) + calcSWSH(-2,lcnt4,mcnt4,theta,phi)*WF(:,wfcnt5,lcnt4-1,mcnt4+5);
        end
    end
end
fprintf('done\n')

% Plots
figure(1);
plot(RO3rawT,RO3rawWF(:,1,1,3),'b')
hold on
plot(RO3rawT,RO3rawWF(:,15,1,3),'r')
hold off
xlabel('t/M')
ylabel('Strain')
title('Raw (2,2) RO3-series Waveforms')
legend('q = 1.00 a = 0.0','q = 4.00 a = 0.0')
grid on

figure(2);
plot(RO3rawT,shiftedWF(:,4,1,3),'b')
hold on
plot(RO3rawT,shiftedWF(:,20,1,3),'r')
hold off
xlabel('t/M')
ylabel('Strain')
title('Timeshifted (2,2) RO3-series Waveforms')
legend('q = 1.00 a = 0.0 oth = 0.090','q = 4.00 a = 0.0 oth = 0.090')
grid on

figure(3);
plot(RO3T,WF(:,4,1,3),'b')
hold on
plot(RO3T,WF(:,20,1,3),'r')
hold off
xlabel('t (s)')
ylabel('Strain')
title('Processed (2,2) RO3-series Wavefroms')
legend('q = 1.00 a = 0.0 oth = 0.90','q = 4.00 a = 0.0 oth = 0.90')
grid on

figure(4);
plot(RO3T,RO3fullWF(:,4),'b')
hold on
plot(RO3T,RO3fullWF(:,8),'r')
plot(RO3T,RO3fullWF(:,20),'g')
hold off
xlabel('t (s)')
ylabel('Strain')
title('RO3-series Waveforms, all modes')
grid on
legend('q = 1.50 a = 0.60 oth = 0.090','q= 2.00 a = 0.60 oth = 0.090', 'q = 4.00 a = 0.60 oth = 0.090')

% Save waveform matrix
waveform_savefile = 'data_RO3-series';
save(waveform_savefile,'RO3fullWF');
time_savefile = 'time_RO3-series';
save(time_savefile,'RO3T');

% Clean up workspace
clearvars -except RO3rawWF RO3rawT RO3T RO3fullWF qvalue avalue
toc

% Next, use RO3_allmodes_pca.m to calculate the principle components and beta
% values