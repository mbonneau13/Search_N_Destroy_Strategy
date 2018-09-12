function [Psi]=Compute_Psi(A,O,E,alpha,theta,nA,TrueObs);
% Computes the singleton potentials Psi_ij(x_ij) of the HMRF.

% INPUT :
%        - A : action matrix (a_ij) : a_ij is the number of observations of
%          cell (i,j) so far.
%        - O : observation matrix (o_ij)
%        - E : Eradication matrix (e_ij=1 if ij was eradicated last year)
%        - alpha : parameters values of the hidden MRF (singletons)
%          alpha(1) if no eradication, alpha(2) if eradication
%        - beta : correlation coefficient of the MRF
%        - theta : current observation probabilities
%          theta(1) if passive search, theta(2) if active search
%        - nA = number of observation steps so far
% OUTPUT :
%        - Psi(i,j,x) : the value of psi_ij(x-1)
%          x=1 if no ants and x=2 if ants.


nA2 = A;                    %number of steps with active search
%nA1 = nA*ones(size(A)) - A; %number of steps without active search
nA1 = ones(size(A)); %passive search only at the first step

% Initialisation of Psi.
Psi=zeros(size(A,1),size(A,2),2);

% Computation of alpha_e_ij.
Alpha_e=E*alpha(2)+(1-E)*alpha(1);



%Computation of theta terms
Theta1 = theta(1)*ones(size(A)) ;
Theta2 = theta(2)*ones(size(A)) ;

Theta1(TrueObs ~= -1) = zeros;
Theta2(TrueObs ~= -1) = zeros;
% Theta3 = theta(3)*ones(size(A)) ;


%%x_ij = 0
% Places where ants have not been observed : Psi(i,j,1)=e if O(i,j)=0.
%Palces where ants have been observed : Psi(i,j,1)=0 if O(i,j)= 1.
% Psi(:,:,1)=(1-O);
O_1 = O;
O_0 = ~O;
Psi(:,:,1) = O_1.*Theta1 + O_0.*(1-Theta1); 



%%x_ij = 1
% Places where ants have been observed : Psi(i,j,2)=1 if O(i,j)= 1.
% Places where ants have not been observed :
% Psi(i,j,2)= exp(alpha_eij)*(1-theta(1))^nA1_ij*(1-theta(2))^nA2_ij  if O(i,j)= 0.
Psi(:,:,2) = O_1.*(exp(Alpha_e + log(1-Theta2))) + O_0.*(exp(Alpha_e + log(Theta2)));
% Psi(:,:,2)=ones(size(O)).*O + (1-O).*( exp(Alpha_e).*((1-Theta1).^nA1.*(1-Theta2) .^nA2)) ;

