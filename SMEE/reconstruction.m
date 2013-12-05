% This script will plot the reconstructed waveforms against the original
% injections

% Load PCs, injections, and reconstruction data
load('/data/nmangini/SMEE_BBH/SMEE/final-MDC_Q-series.mat')
QMDC = MDC_final;
load('/data/nmangini/SMEE_BBH/SMEE/final-MDC_HR-series.mat')
HRMDC = MDC_final;
load('/data/nmangini/SMEE_BBH/SMEE/final-MDC_RO3-series.mat')
RO3MDC = MDC_final;
load('/data/nmangini/SMEE_BBH/SMEE/finalRvectorsPC_Q-series.mat')
QPC = PCs_final*1e-15;
load('/data/nmangini/SMEE_BBH/SMEE/finalRvectorsPC_HR-series.mat')
HRPC = PCs_final*1e-15;
load('/data/nmangini/SMEE_BBH/SMEE/finalRvectorsPC_RO3-series.mat')
RO3PC = PCs_final*1e-15;
SMEE_file = dir('*.mat');
load(SMEE_file.name)

% Reconstruction
avgbeta = mean(postbetas);
time = linspace(0,length(QMDC)/16384,length(QMDC));
recon = 0;
for i = 1:length(avgbeta)
    if strcmp(model,'Q') == 1
        s = QPC(:,i)*avgbeta(i);
    elseif strcmp(model,'HR') == 1
        s = HRPC(:,i)*avgbeta(i);
    elseif strcmp(model,'RO3') == 1
        s = RO3PC(:,1)*avgbeta(i);
    end
    recon = recon + s;
end

% Plotting
plot(time,recon/max(recon(time>0.5)))
hold on
if strcmp(catalogue,'Q') == 1
    plot(time,QMDC(:,wv)/max(QMDC(:,wv)),'r')
elseif strcmp(catalogue,'HR') == 1
    plot(time,HRMDC(:,wv)/max(HRMDC(:,wv)),'r')
elseif strcmp(catalogue,'RO3') == 1
    plot(time,RO3MDC(:,wv)/max(RO3MDC(:,wv)),'r')
end
grid on
xlabel('t(s)','fontweight','bold')
ylabel('Strain','fontweight','bold')
title(strcat('Reconstructed ',catalogue,'Injection with ',model,'Waveform: Bayes = ',num2str(Bayes)),'fontweight','bold')
legend('Reconstruction','Injection')
ylim([-1.5,1.5])