function [x] = dQ_a2(PG,E,Nzeros,Nones,alpha1,alpha2, beta);

% this function computes the derivate of the Q function of EM
% according to alpha_2

% INPUT :
%      - PG : P_G^ij(Z'_ij=1|Z1,a_ij,o_ij,Alpha,beta,Theta)
%      - E : Eradication matrix
%      - Nzeros : matrix of same size as Z
%               Nzeros(i,j) is the number of zeros in the neighbours of
%               (i,j)
%      - Nones : matrix of same size as Z
%               Nones(i,j) is the number of ones in the neighbours of
%               (i,j)
%      - alpha, beta : parameters values of the hidden MRF
%        alpha1 if no eradication, alpha2 if eradication.
% OUTPUT :
%      - x = dQ/dalpha2



Alpha = alpha1*ones(size(E));
Alpha(find(E)) = alpha2;

% computation of dQ/dalpha1
M = exp(Alpha+beta*Nones)./(exp(beta*Nzeros)+exp(Alpha+beta*Nones));

x = sum(sum(PG(find(E))))-sum(sum(M(find(E))));

