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
    avgQmatchfull(k) = mean(Qmatchfull(:,k));
end

% Plots
plot(1:13,minQmatchfull,'bo--')
hold on
plot(1:13,avgQmatchfull,'ro--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Minimum Match vs. Number of PCs, Q-series')
legend('Minimum Match','Average Match')
axis([1 13 0 1])
grid on
hold off

% Clean up work space
clearvars -except QrawT QrawWF QT QfullWF qvalue QvectorsPC QreconstructedfullWF...
    Qbeta minbeta maxbeta minQmatchfull norm_QreconstructedfullWF
