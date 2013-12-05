% This script creates a matrix of waveforms for the HR-series catalogue,
% aligns the waveforms, zeros the data after the ringdown signal, shortens
% the signal duration, and combines the spherical harmonic modes

tic
% Set up variables and solution arrays
maxlength = 31238;
qvalue = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1.5; 1.5; 1.5; 2; 2.5; 4];
avalue = [0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0; 0.2; 0.4; 0; 0; 0];
mvalue = [1, 2, -1, -2; 2, 3, -2, -3; 3, 4, -3, -4];
q1mvalue = [0, 2, 0, -2; 2, 0, -2, -0; 0, 4, 0, -4]; % used in combining spherical harmonic modes
qvaluestr = cellstr(num2str(qvalue,'%4.2f')); % String of q values
avaluestr = cellstr(num2str(avalue,'%4.1f')); % String of a values
HRrawWF = zeros(maxlength,length(qvalue),3,4); % rawWF = (length of waveform)X(q value)X(l value)X(m value)
theta = 0*pi/180;
phi = 0*pi/180;

% Read in each waveform
fprintf('Reading in waveforms...')
for wfcnt1 = 1:length(qvalue)
    for lcnt1 = 2:4;
        for mcnt1 = 1:3
            fname = strcat('/Users/Nick/LIGO/Glasgow_2013/Waveforms/HR-series/D11_q',qvaluestr{wfcnt1},...
                '_a',avaluestr{wfcnt1},'_m200/rStrain_FFI_l', num2str(lcnt1), '_m',...
                num2str(mvalue(lcnt1-1,mcnt1)) ,'_r75.00.ampphase.asc');
            importfile(fname);
            HRrawWF(1:length(data(:,2)),wfcnt1,lcnt1-1,mcnt1) = data(:,2);
        end
    end
end
HRrawT = data(:,1);
fprintf('done\n')

% Set up matrices, solution arrays, and variables for waveform processing
x = size(HRrawWF);
L = x(1);
s = x(2);
l_size = x(3);
m_size = x(4);
peak = zeros(s,1);
peaktime = zeros(1,s);

% Find the peak time in each (2,2) waveform
fprintf('Finding the peak time of the (2,2) mode...')
starttime = 0;
i = 0;
while starttime == 0
    i = i + 1;
    if floor(HRrawT(i)) == 200 % finds when t = 200 ms
        starttime = i;
    end
end

for wfcnt2 = 1:s
    for tcnt1 = 1:L
        if HRrawWF(tcnt1,wfcnt2,1,2)^2 > peak(wfcnt2)
            peak(wfcnt2) = HRrawWF(tcnt1,wfcnt2,1,2).^2;
            peaktime(wfcnt2) = tcnt1;
        end
    end
end
starttime = starttime + (max(peaktime) - min(peaktime));
fprintf('done\n')

% Align waveforms along thier peak
fprintf('Aligning waveforms...')
shiftedWF = zeros(L,s,3,4);
for wfcnt3 = 1:s
    for lcnt2 = 1:3
        for mcnt2 = 1:4
            for tcnt2 = (max(peaktime)-peaktime(wfcnt3)+1):L
                shiftedWF(tcnt2,wfcnt3,lcnt2,mcnt2) = HRrawWF(tcnt2-(max(peaktime)-peaktime(wfcnt3)),wfcnt3,lcnt2,mcnt2);
                if HRrawT(tcnt2) > (max(peaktime) + 200)
                    shiftedWF(tcnt2,wfcnt3,lcnt2,mcnt2) = 0; % sets data after ringdown to 0
                end
            end
        end
    end
end
fprintf('done\n')
plot(HRrawT,shiftedWF(:,1,1,2))

% Removes data prior to 200 ms
fprintf('Removing data prior to 200 ms...')
WF = zeros(length(shiftedWF)-starttime,s,x(3),x(4));
for i = starttime:length(shiftedWF)
    for wfcnt4 = 1:s
        for lcnt3 = 1:l_size
            for mcnt3 = 1:m_size
                WF(i-starttime+1,wfcnt4,lcnt3,mcnt3) = shiftedWF(i,wfcnt4,lcnt3,mcnt3);
            end
        end
    end
    HRT(i-starttime+1) = HRrawT(i);
end

% Convert from NR units to physical time units
HRT = HRT';
HRT = HRT*4.92549095e-6*250; %250 M_sun
fprintf('done\n')

% Combine all spherical harmonic modes
fprintf('Combining spherical harmonic modes...')
HRfullWF = zeros(length(WF),s);
for wfcnt5 = 1:9
    for lcnt4 = 1:3;
        for mcnt4 = 1:4;
            if q1mvalue(lcnt4,mcnt4) ~= 0    
                HRfullWF(:,wfcnt5) = HRfullWF(:,wfcnt5) + calcSWSH(-2,(lcnt4+1),mvalue(lcnt4,mcnt4),theta,phi)*WF(:,wfcnt5,lcnt4,mcnt4);
            end
        end
    end
end
for wfcnt6 = 10:s;
    for lcnt5 = 1:3;
        for mcnt5 = 1:4;
            HRfullWF(:,wfcnt6) = HRfullWF(:,wfcnt6) + calcSWSH(-2,(lcnt5+1),mvalue(lcnt5,mcnt5),theta,phi)*WF(:,wfcnt6,lcnt5,mcnt5);
        end
    end
end
fprintf('done\n')

% Plots
figure(1);
plot(HRrawT,HRrawWF(:,1,1,2),'b')
hold on
plot(HRrawT,HRrawWF(:,15,1,2),'r')
hold off
xlabel('t/M')
ylabel('Strain')
title('Raw (2,2) HR-series Waveforms')
legend('q = 1.00 a = 0.0','q = 4.00 a = 0.0')
grid on

figure(2);
plot(HRrawT,shiftedWF(:,1,1,2),'b')
hold on
plot(HRrawT,shiftedWF(:,15,1,2),'r')
hold off
xlabel('t/M')
ylabel('Strain')
title('Timeshifted (2,2) HR-series Waveforms')
legend('q = 1.00 a = 0.0','q = 4.00 a = 0.0')
grid on

figure(3);
plot(HRT,WF(:,1,1,2),'b')
hold on
plot(HRT,WF(:,15,1,2),'r')
hold off
xlabel('t (s)')
ylabel('Strain')
title('Processed (2,2) HR-series Wavefroms')
legend('q = 1.00 a = 0.0','q = 4.00 a = 0.0')
grid on

figure(4);
plot(HRT,HRfullWF(:,1),'b')
hold on
plot(HRT,HRfullWF(:,13),'r')
plot(HRT,HRfullWF(:,15),'g')
legend('q = 1.00 a = 0.0','q = 2.50 a = 0.0','q = 4.00 a = 0.0')
hold off
xlabel('t (s)')
ylabel('Strain')
title('HR-series Waveforms, all modes')
grid on

% Save waveform matrix
waveform_savefile = 'data_HR-series';
save(waveform_savefile,'HRfullWF');
time_savefile = 'time_HR-series';
save(time_savefile,'HRT');

% Clean up workspace
clearvars -except HRrawWF HRrawT HRT HRfullWF qvalue avalue
toc

% Next, use HR_allmodes_pca.m to calculate the principle components and beta
% values