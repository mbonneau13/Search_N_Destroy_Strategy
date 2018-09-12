function [messages,epsilon]=Compute_messages(Psi,beta,precision,Tmax,ZN);
% Computes the fixed-point messages of the sum-product algorithm

% INPUT :
%        - Psi(i,j,x) : the value of psi_ij(x-1).
%        - beta : correlation coefficient of the MRF, used to compute
%          psi_(ij,kl) = exp(beta*eq(x_ij,x_kl))
%        - precision : threshold on maximum absolute differences between
%          two successive time steps
%        - Tmax : maximum number of iterations
% OUTPUT :
%        - messages(d,i,j,x) : fixed-point messages (d is the direction of
%          message)
%        - epsilon : maximum absolute difference between to successive
%          message computations for each time step.

% initialization of all messages to 1/2.
messages=ones(4,size(Psi,1),size(Psi,2),size(Psi,3))/2;
%messages=ones(8,size(Psi,1),size(Psi,2),size(Psi,3))/2;

%% We fix the message to one out of the border
[l , c] = find(ZN == 0);
for i = 1:length(l)
    messages(1,l(i),c(i),:) = ones;
    messages(2,l(i),c(i),:) = ones;
    messages(3,l(i),c(i),:) = ones;
    messages(4,l(i),c(i),:) = ones;
end

epsilon=[ ];
t=0;
eps=precision+1;
while ((t<Tmax)&(eps>precision))
    t=t+1;
    [messages,eps]=Parallel_Update(messages,Psi,beta,4,ZN);
    %[messages,eps]=Parallel_Update(messages,Psi,beta,8);
    
    epsilon=[epsilon eps];
end

