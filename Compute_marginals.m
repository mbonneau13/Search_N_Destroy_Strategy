function [marginals]=Compute_marginals(A,O,E,alpha,beta,theta,nA,ZN,TrueObs)
% Computes the marginal probabilities of occupation of the conditional MRF
% from parameters alpha, beta, theta and from the data (A, O, E).

% INPUT :
%        - A : action matrix (a_ij) : a_ij is the number of observations of
%        cell (i,j) so far.
%        - O : observation matrix (o_ij)
%        - E : Eradication matrix (e_ij=1 if ij was eradicated last year)
%        - alpha : parameters values of the hidden MRF (singletons)
%          alpha(1) if no eradication, alpha(2) if eradication
%        - beta : correlation coefficient of the MRF
%        - theta : current observation probabilities
%          theta(1) if passive search, theta(2) if active search
%        - nA = number of observation steps so far
% OUTPUT :
%        - marginals(i,j) : Probability of occupation of cell (i,j) (x_ij=2).


nA = 1;



% Psi(i,j,x) : the value of psi_ij(x-1).
[Psi]=Compute_Psi(A,O,E,alpha,theta,nA,TrueObs);

P0 = squeeze(Psi(:,:,1));
P1 = squeeze(Psi(:,:,2));

P0(TrueObs == 0) = ones;
P0(TrueObs == 1) = zeros;

P1(TrueObs == 0) = zeros;
P1(TrueObs == 1) = ones;

P0(ZN == 0) = ones;
P1(ZN == 0) = zeros;

Psi(:,:,1) = P0;
Psi(:,:,2) = P1;

% messages(d,i,j,x) : fixed-point messages.
% epsilon : maximum absolute difference between to successive
% message computations for each time step.
precision=0.001;
Tmax=1000;


[messages,epsilon]=Compute_messages(Psi,beta,precision,Tmax,ZN);

P=zeros(size(Psi));
for i=1:size(Psi,1)
    for j=1:size(Psi,2)
        for xij=1:2
            P(i,j,xij) = Psi(i,j,xij)*prod(messages(:,i,j,xij));
        end
        P(i,j,:)=P(i,j,:)/sum(P(i,j,:));
    end
end

marginals=zeros(size(P,1),size(P,2));
marginals=P(:,:,2);