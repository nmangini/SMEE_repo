% This script calculates the principle components of the Q-series catalogue
% and the associated beta values.

tic
% Variables
x = size(RO3fullWF);
L = x(1);
s = x(2);

% D = USV'
%   D is an mXn matrix containing the data
%   U and V are orthogonal matrices
%   S is a diagonal matrix
fprintf('Calculating principle components...')
D = RO3fullWF;

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

% Step #5: Calculate the beta values by projecting A onto U
fprintf('Calculating beta values...')
RO3beta = zeros(L,s);
minbeta = zeros(s,1);
maxbeta = zeros(s,1);
for i = 1:s
    for j = 1:s
        RO3beta(i,:) = A(:,i)' * vectorsPC;
        RO3reconstructedfullWF(:,i,j) = (RO3beta(i,1:j)*(vectorsPC(:,1:j)'))';
        minbeta(i,:) = min((RO3beta(i,:)));
        maxbeta(i,:) = max((RO3beta(i,:)));
    end
end 
fprintf('done\n')

% Plots
figure(1);
for i = 1:3
    subplot(3,1,i)
    plot(RO3T,vectorsPC(:,i))
    grid on
    xlabel('t (s)')
    ylabel('Strain')
end
suptitle('First 3 (normalized) Principal Components');

figure(2);
for i = 1:3
    subplot(3,1,i)
    plot(qvalue,RO3beta(i,:),'x-')
    xlabel('q')
    ylabel(strcat({'\beta '},num2str(i)))
    grid on
end
suptitle('RO3-series \beta s')

% Save
evector_savefile = 'vectors_RO3-series';
save(evector_savefile,'vectorsPC');
betas_savefile = 'betas_RO3-series';
save(betas_savefile,'RO3beta', 'maxbeta', 'minbeta');


% Clean up work space
clearvars -except RO3rawT RO3rawWF RO3T RO3fullWF qvalue RO3vectorsPC RO3reconstructedfullWF RO3beta minbeta maxbeta
toc

% Next use RO3_allmodes_match.m to calculate the match