% This script calculates the match versus the number of pcs

% Normalize the full and reconstructed waveforms and calculate the match
for i = 1:length(qvalue);
    norm_QfullWF(:,i) = QfullWF(:,i)/(norm(QfullWF(:,i)));
    for j = 1:length(qvalue);
        norm_QreconstructedfullWF(:,i,j) = QreconstructedfullWF(:,i,j)/(norm(QreconstructedfullWF(:,i,j)));
        Qmatchfull(i,j) = dot(norm_QfullWF(:,i),norm_QreconstructedfullWF(:,i,j));      
    end
end

for k = 1:13
    minQmatchfull(k) = min(Qmatchfull(:,k));
end

% Plots
plot(1:13,minQmatchfull,'bo--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Minimum Match vs. Number of PCs, Q-series')
axis([1 13 0 1])
grid on

% Clean up work space
clearvars -except QrawT QrawWF QT QfullWF qvalue QvectorsPC QreconstructedfullWF...
    Qbeta minbeta maxbeta minQmatchfull norm_QreconstructedfullWF