% Load in waveforms and principle components
load('Q-series.mat')
load('recon_Q.mat')
load('HR-series.mat')
load('recon_HR.mat')
load('RO3-series.mat')
load('recon_RO3.mat')

% Normalize the full and reconstructed waveforms and calculate the match
% Q-series
for i = 1:13;
    norm_Q(:,i) = Q_final(:,i)/(norm(Q_final(:,i)));
    for j = 1:13;
        norm_Q_recon(:,j) = Q_recon(:,j)/(norm(Q_recon(:,j)));
        QQmatchfull(i,j) = dot(norm_Q(:,i),norm_Q_recon(:,j));      
    end
end
for k = 1:13
    minQQmatchfull(k) = min(QQmatchfull(:,k));
    avgQQmatchfull(k) = mean(QQmatchfull(:,k));
end
for i = 1:13;
    norm_Q(:,i) = Q_final(:,i)/(norm(Q_final(:,i)));
    for j = 1:13;
        norm_HR_recon(:,j) = HR_recon(:,j)/(norm(HR_recon(:,j)));
        QHRmatchfull(i,j) = dot(norm_Q(:,i),norm_HR_recon(:,j));      
    end
end
for k = 1:13
    minQHRmatchfull(k) = min(QHRmatchfull(:,k));
    avgQHRmatchfull(k) = mean(QHRmatchfull(:,k));
end
for i = 1:13;
    norm_Q(:,i) = Q_final(:,i)/(norm(Q_final(:,i)));
    for j = 1:13;
        norm_RO3_recon(:,j) = RO3_recon(:,j)/(norm(RO3_recon(:,j)));
        QRO3matchfull(i,j) = dot(norm_Q(:,i),norm_RO3_recon(:,j));      
    end
end
for k = 1:13
    minQRO3matchfull(k) = min(QRO3matchfull(:,k));
    avgQRO3matchfull(k) = mean(QRO3matchfull(:,k));
end

% Plotting
figure()
plot(1:13,minQQmatchfull,'bo--')
hold on
plot(1:13,avgQQmatchfull,'ro--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Match vs. Number of Q PCs, Q-series')
legend('Minimum Match','Average Match')
axis([1 13 0 1])
grid on
hold off

figure()
plot(1:13,minQHRmatchfull,'bo--')
hold on
plot(1:13,avgQHRmatchfull,'ro--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Match vs. Number of HR PCs, Q-series')
legend('Minimum Match','Average Match')
axis([1 13 0 1])
grid on
hold off

figure()
plot(1:13,minQRO3matchfull,'bo--')
hold on
plot(1:13,avgQRO3matchfull,'ro--')
xlabel('Number of PCs')
ylabel('Minimum Match')
title('Match vs. Number of RO3 PCs, Q-series')
legend('Minimum Match','Average Match')
axis([1 13 0 1])
grid on
hold off