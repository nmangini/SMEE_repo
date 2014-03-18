function collect_results(Signal, Waveform, Identifier)
% COLLECT_RESULTS This script will move through directories with a given
%   signal and model and gather the Bayes factor from each smee_results.txt
%   file. It will the then plot the Bayes factor for each SMEE run. The
%   directory syntax is <Signal><Waveform>-<Model>-<Identifier>
%       Signal -- The injected waveform (Q, HR, or RO3)
%       Waveform -- The number of the waveform injected:
%                   1:13 for Q-series
%                   1:15 for HR-series
%                   1:20 for RO3-series
%       Model -- The model waveform used to reconstruct the signal (Q, HR, or R03)
%       Identifier -- The tag associated with SMEE runs:
%                   seed - Random seed runs

% Input Managment
wv = num2str(Waveform);
if strcmp(Identifier,'seed') == 1
    rseed = 1:20;
    Bayes = zeros(3,length(rseed));
elseif strcmp(Identifier,'seed') == 0
    error('Incorrect Identifier')
end

% Variables
ResultsDir = pwd;
%ResultsDir = '/data/nmangini/SMEE_BBH/SMEE/Results/';
model1 = 'Q';
model2 = 'HR';
model3 = 'RO3';
ModelArray = {model1,model2,model3};

% Collect Results
for mod = 1:3
    for dir = 1:length(rseed)
        %load(strcat(ResultsDir,Signal,wv,'-',ModelArray{mod},'-',Identifier,...
        %num2str(rseed(dir)),'/smee_results_6PCs.txt'));
        %load(strcat(ResultsDir,'/',Signal,wv,'-',ModelArray{mod},'-',Identifier,...
        %num2str(rseed(dir)),'-snr10-pcs8','/smee_results_8PCs.txt'));
        load(strcat(ResultsDir,'/',Signal,wv,'-',ModelArray{mod},'-',Identifier,...
        num2str(rseed(dir)),'/smee_results.txt'));
        Bayes(mod,dir) = smee_results(3);
    end
end

% Plotting
p1 = plot(rseed,Bayes(1,:),'b');
hold on
p2 = plot(rseed,Bayes(2,:),'r');
%plot(rseed,Bayes(2,:),'r');
p3 = plot(rseed,Bayes(3,:),'g');
grid on
title(strcat(Signal,'-series Injection'))
if strcmp(Identifier,'seed') == 1
    xlabel('Random Seed')
end
ylabel('Bayes Factor')
legend([p1,p2,p3],{'Q-series','HR-series','RO3-series'},'Location','NorthWest')
