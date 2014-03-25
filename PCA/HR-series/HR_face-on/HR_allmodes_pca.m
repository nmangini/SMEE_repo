% This script calculates the principle components of the HR-series catalogue
% and the associated beta values.

tic
% Variables
x = size(HRfullWF);
L = x(1);
s = x(2);

% D = USV'
%   D is an mXn matrix containing the data
%   U and V are orthogonal matrices
%   S is a diagonal matrix
fprintf('Calculating principle components...')
D = HRfullWF;

% Step #1: Calculate the adjusted data set (A)
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

% % Step #5: Calculate the beta values by projecting A onto U
fprintf('Calculating beta values...')
HRbeta = zeros(L,s);
minbeta = zeros(s,1);
maxbeta = zeros(s,1);
for i = 1:s
    for j = 1:s
        HRbeta(i,:) = A(:,i)' * vectorsPC;
        HRreconstructedfullWF(:,i,j) = (HRbeta(i,1:j)*(vectorsPC(:,1:j)'))';
        minbeta(i,:) = min((HRbeta(i,:)));
        maxbeta(i,:) = max((HRbeta(i,:)));
    end
end 
fprintf('done\n')

% Plots
figure(1);
for i = 1:3
    subplot(3,1,i)
    plot(HRT,vectorsPC(:,i))
    grid on
    xlabel('t (s)')
    ylabel('Strain')
end
suptitle('First 3 (normalized) Principal Components');

figure(2);
for i = 1:3
    subplot(3,1,i)
    plot(avalue(1:9),HRbeta((1:9),i),'x-')
    xlabel('a')
    ylabel(strcat({'\beta '},num2str(i)))
    grid on
end
suptitle('HR-series \beta s')

% Save
evector_savefile = 'vectors_HR-series';
save(evector_savefile,'vectorsPC');
betas_savefile = 'betas_HR-series';
save(betas_savefile,'HRbeta', 'maxbeta', 'minbeta');
reconstruction_savefile = 'recon_HR-series';
save(reconstruction_savefile,'HRreconstructedfullWF');

% Clean up work space
clearvars -except HRrawT HRrawWF HRT HRfullWF qvalue HRvectorsPC HRreconstructedfullWF HRbeta minbeta maxbeta
toc

% Next use HR_allmodes_match.m to calculate the match