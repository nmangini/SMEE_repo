% This script calculates the match versus the number of pcs

% Normalize the full and reconstructed waveforms and calculate the match
for i = 1:length(qvalue);
    norm_RO3fullWF(:,i) = RO3fullWF(:,i)/(norm(RO3fullWF(:,i)));
    for j = 1:length(qvalue);
        norm_RO3reconstructedfullWF(:,i,j) = RO3reconstructedfullWF(:,i,j)/(norm(RO3reconstructedfullWF(:,i,j)));
        RO3matchfull(i,j) = dot(norm_RO3fullWF(:,i),norm_RO3reconstructedfullWF(:,i,j));      
    end
end

for k = 1:13
    minRO3matchfull(k) = min(RO3matchfull(:,k));
end

% Plots
plot(1:13,minRO3matchfull,'bo--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Minimum Match vs. Number of PCs, RO3-series')
axis([1 13 0 1])
grid on

% Clean up work space
clearvars -except RO3rawT RO3rawWF RO3T RO3fullWF qvalue RO3vectorsPC RO3reconstructedfullWF RO3beta minbeta maxbeta minRO3matchfull