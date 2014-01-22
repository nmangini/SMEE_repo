% This script calculates the principle components of the Q-series catalogue
% and the associated beta values.

% Variables
tic
qvalue = [1.00; 1.15; 1.30; 1.45; 1.50; 1.60; 1.75; 1.90; 2.00; 2.05; 2.20; 2.35; 2.50];
x = size(QfullWF);
L = x(1);
s = x(2);

% D = USV'
%   D is an mXn matrix containing the data
%   U and V are orthogonal matrices
%   S is a diagonal matrix
D = QfullWF;

% Step #1: Calculate the adjusted data set (A)
fprintf('Calculating principle components...')
% M = mean(D,2);
% A = zeros(L,s);
% for i = 1:L
%     A(i,:) = D(i,:) - M(i,:);
% end
A = D;

% Step #2: Calculate covariance matrix, C, from A
C = (A'*A)/(L-1);

% Step #3: Calculate eigenvectors (V) and eigenvalues (S) of C
[V,lambda2] = eig(C);
S = sqrt(lambda2);

% Step #3: Organize eigenvalues in descending order and thier corresponding
% eigenvectors
evalues = diag(S);
[evalues,index] = sort(evalues,'descend');
V = V(:,index);

% Step #4: Compute the eigenvectors of the real covariance matrix U and
% normalize
U = A*V;
for i = 1:s
    U(:,i) = U(:,i)/norm(U(:,i));
end
vectorsPC = U;
fprintf('done\n')

% Step #5: Calculate the beta values by projecting A onto U
fprintf('Calculating Beta values...')
Qbeta = zeros(L,s);
minbeta = zeros(s,1);
maxbeta = zeros(s,1);
for i = 1:s
    for j = 1:s
        Qbeta(i,:) = A(:,i)' * vectorsPC;
        QreconstructedfullWF(:,i,j) = (Qbeta(i,1:j)*(vectorsPC(:,1:j)'))';
        minbeta(i,:) = min((Qbeta(i,:)));
        maxbeta(i,:) = max((Qbeta(i,:)));
    end
end
fprintf('done\n')

% Plots
%figure(1);
%for i = 1:3
%    subplot(3,1,i)
%    plot(QT,vectorsPC(:,i))
%    grid on
%    xlabel('t (s)')
%    ylabel('Strain')
%    grid on
%end
%suptitle('First 3 (normalized) Principal Components');

%figure(2);
%for i = 1:3
%    subplot(3,1,i)
%    plot(qvalue,Qbeta(i,:),'x-')
%    xlabel('q')
%    ylabel(strcat({'\beta '},num2str(i)))
%    grid on
%end
%suptitle('Q-series \beta s')

% Save
evector_savefile = 'vectors_Q-series';
save(evector_savefile,'vectorsPC');
betas_savefile = 'betas_Q-series';
save(betas_savefile,'Qbeta', 'maxbeta', 'minbeta');

% Clean up work space
%clearvars -except QrawT QrawWF QT QfullWF qvalue vectorsPC QreconstructedfullWF Qbeta minbeta maxbeta M
toc

% Next use Q_allmodes_match.m to calculate the match
