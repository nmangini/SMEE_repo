% This script calculates the match versus the number of pcs

% Normalize the full and reconstructed waveforms and calculate the match
for i = 1:length(qvalue);
    norm_HRfullWF(:,i) = HRfullWF(:,i)/(norm(HRfullWF(:,i)));
    for j = 1:length(qvalue);
        norm_HRreconstructedfullWF(:,i,j) = HRreconstructedfullWF(:,i,j)/(norm(HRreconstructedfullWF(:,i,j)));
        HRmatchfull(i,j) = dot(norm_HRfullWF(:,i),norm_HRreconstructedfullWF(:,i,j));      
    end
end

for k = 1:13
    minHRmatchfull(k) = min(HRmatchfull(:,k));
    avgHRmatchfull(k) = mean(HRmatchfull(:,k));
end

% Plots
plot(1:13,minHRmatchfull,'bo--')
hold on
plot(1:13,avgHRmatchfull,'ro--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Minimum Match vs. Number of PCs, HR-series')
legend('Minimum Match','Average Match')
axis([1 13 0 1])
grid on
hold off

% Clean up work space
clearvars -except HRrawT HRrawWF HRT HRfullWF qvalue HRvectorsPC HRreconstructedfullWF ...
    HRbeta minbeta maxbeta minHRmatchfull norm_HRreconstructedfullWF