% This script creates a matrix of waveforms for the Q-series catalogue,
% aligns the waveforms, zeros the data after the ringdown signal, shortens
% the signal duration, and combines the spherical harmonic modes

tic
% Set up variables and solution arrays
qvalue = [1; 1.15; 1.3; 1.45; 1.5; 1.6; 1.75; 1.9; 2; 2.05; 2.20; 2.35; 2.50];
mvalue = [1, 2, -1, -2; 2, 3, -2, -3; 3, 4, -3, -4];
QrawWF = zeros(10959,length(qvalue),3,4); % rawWF = (length of waveform)X(q value)X(l value)X(m value)
theta = 0*pi/180;
phi = 0*pi/180;
MASS = 250; % 250 M_sun


% Read in each waveform
fprintf('Reading in waveforms...')
qvaluestr = cellstr(num2str(qvalue,'%4.2f')); % String of q values
for qcnt1 = 1:length(qvalue)
    for lcnt1 = 2:4;
        for mcnt1 = 1:3
            fname = strcat('~/SMEE_repo/Waveforms/Q-series/D10_a0.0_q',...
                qvaluestr{qcnt1}, '_m103_Qs/rStrain_FFI_l', num2str(lcnt1), '_m',...
                num2str(mvalue(lcnt1-1,mcnt1)) ,'_r75.00.asc');
            importfile(fname);
            QrawWF(:,qcnt1,lcnt1-1,mcnt1) = data(:,2);
        end
    end
end
QrawT = data(:,1);
fprintf('done\n')

% Set up matrices, solution arrays, and variables for waveform processing
x = size(QrawWF);
L = x(1);
s = x(2);
l_size = x(3);
m_size = x(4);
peak = zeros(1,s);
peaktime = zeros(1,s);

% Find the peak time in each (2,2) waveform
fprintf('Finding the peak time of the (2,2) mode...')
starttime = 0;
i = 0;
while starttime == 0
    i = i + 1;
    if floor(QrawT(i)) == 100 % finds when t = 200 ms
        starttime = i;
    end
end

for qcnt2 = 1:s
    for tcnt1 = 1:L
        if QrawWF(tcnt1,qcnt2,1,2)^2 > peak(qcnt2)
            peak(qcnt2) = QrawWF(tcnt1,qcnt2,1,2).^2;
            peaktime(qcnt2) = tcnt1;
        end
    end
end
starttime = starttime + (peaktime(13) - peaktime(1));
fprintf('done\n')

% Align waveforms along thier peak
fprintf('Aligning waveforms...')
shiftedWF = zeros(L,s,3,4);
for qcnt3 = 1:s
    for lcnt2 = 1:l_size
        for mcnt2 = 1:m_size
            for tcnt2 = (peaktime(13)-peaktime(qcnt3)+1):L
                shiftedWF(tcnt2,qcnt3,lcnt2,mcnt2) = QrawWF(tcnt2-(peaktime(13) - peaktime(qcnt3)),qcnt3,lcnt2,mcnt2);
                if QrawT(tcnt2) > (QrawT(peaktime(qcnt2)) + 100)
                    shiftedWF(tcnt2,qcnt3,lcnt2,mcnt2) = 0; % sets data after ringdown to 0
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
    for qcnt4 = 1:s
        for lcnt3 = 1:l_size
            for mcnt3 = 1:m_size
                WF(i-starttime+1,qcnt4,lcnt3,mcnt3) = shiftedWF(i,qcnt4,lcnt3,mcnt3);
            end
        end
    end
    QT(:,i-starttime+1) = QrawT(i); 
end

% Convert from NR units to physical time units
QT = QT';
QT = QT*4.92549095e-6*MASS;
fprintf('done\n')

% Combine all spherical harmonic modes
fprintf('Combining spherical harmonic modes...')
q1mvalue = [0, 2, 0, -2; 2, 0, -2, -0; 0, 4, 0, -4]; % used in combining spherical harmonic modes
QfullWF = zeros(length(WF),s);
for lcnt4 = 1:l_size;
    for mcnt4 = 1:m_size;
        if q1mvalue(lcnt4,mcnt4) ~= 0    
            QfullWF(:,1) = QfullWF(:,1) + calcSWSH(-2,(lcnt4+1),mvalue(lcnt4,mcnt4),theta,phi)*WF(:,1,lcnt4,mcnt4);
        end
   end
end
for qcnt5 = 2:s;
    for lcnt5 = 1:l_size;
        for mcnt5 = 1:m_size;
            QfullWF(:,qcnt5) = QfullWF(:,qcnt5) + calcSWSH(-2,(lcnt5+1),mvalue(lcnt5,mcnt5),theta,phi)*WF(:,qcnt5,lcnt5,mcnt5);
        end
    end
end
fprintf('done\n')

% Plots
%figure(1);
%plot(QrawT,QrawWF(:,1,1,2),'b')
%hold on
%plot(QrawT,QrawWF(:,13,1,2),'r')
%hold off
%xlabel('t/M')
%ylabel('Strain')
%title('Raw (2,2) Q-series Waveforms')
%legend('q = 1.00','q = 2.50')
%grid on

%figure(2);
%plot(QrawT,shiftedWF(:,1,1,2),'b')
%hold on
%plot(QrawT,shiftedWF(:,13,1,2),'r')
%hold off
%xlabel('t/M')
%ylabel('Strain')
%title('Timeshifted (2,2) Q-series Waveforms')
%legend('q = 1.00','q = 2.50')
%grid on

%figure(3);
%plot(QT,WF(:,1,1,2),'b')
%hold on
%plot(QT,WF(:,13,1,2),'r')
%hold off
%xlabel('t (s)')
%ylabel('Strain')
%title('Processed (2,2) Q-series Wavefroms')
%legend('q = 1.00','q = 2.50')
%grid on

%figure(4);
%plot(QT,QfullWF(:,1),'b')
%hold on
%plot(QT,QfullWF(:,5),'r')
%plot(QT,QfullWF(:,13),'g')
%legend('q = 1.00','q = 1.50','q = 2.50')
%hold off
%xlabel('t (s)')
%ylabel('Strain')
%title('Q-series Waveforms, all modes')
%grid on

% Save waveform matrix
waveform_savefile = 'data_Q-series';
save(waveform_savefile,'QfullWF');
time_savefile = 'time_Q-series';
save(time_savefile,'QT');

% Clean up work space
%clearvars -except QrawT QrawWF QT QfullWF qvalue
toc

% Next, use Q_allmodes_pca.m to calculate the principle components and beta
% values
