function [Z1,PI1,PG1,W1,Nzeros,Nones] = Gibbs_Sampler(Z,A,O,E,alpha,beta,Theta,q,TrueObs,ZN);


% Gibbs Sampler : generate a realisation of the Markov Random Field
% conditional to the observations.
% The probabilistic distribution is defined by parameters alpha, beta, theta,
% and by action, observation,  eradication (A,O,E).

% use Gibbs_Sampler_3theta if model with 3 theta
% values : without active search in rural area, without active search in urban area, 
% and with active search



% INPUT 
%      - Z : initial value of (z_ij), the hidden variables matrix
%      - A : action matrix (a_ij) 
%      - O : observation matrix (o_ij)
%      - E : Eradication matrix (e_ij=1 if ij was eradicated last year)
%      - alpha, beta : parameters values of the hidden MRF
%        alpha(1) if no eradication, alpha(2) if eradication.
%      - theta : vector of observation probabilities : theta(1)=theta_0 and
%      theta(2)=theta_1
%      - q : number of Gibbs sampler passes 
% OUTPUT
%      - Z1 : new sample of Z after q passes
%      - PI1 : P_I^ij(Z'_ij=1|Z1,alpha,beta)
%      - PG1 : P_G^ij(Z'_ij=1|Z1,a_ij,o_ij,alpha,beta,Theta)
%      - W1 : local weights matrix used in PI1 computation
%      - Nzeros, Nones : number of neighbours
%        in state 1 and in state 0 of each cell (i,j) on Z1

%example
%A, O and E are provided by the field data or can be toy matrices
%
% Z = rand(size(A));
% Z = (Z<0.5);
% Z(find(O)) = 1;
%[Z1,PI1,PG1,W1,Nzeros,Nones] = Gibbs_Sampler(Z,A,O,E,[-1 0],0.5,[0.2 0.8],200);



Z1=Z;
PI1=zeros(size(Z1));
PG1=zeros(size(Z1));
W1=zeros(size(Z1));
% Z1(find(O))=1;
%% Verifier que l on respecte bien les vraies observations.
I = find(TrueObs ~= -1);
PI1(I) = TrueObs(I);
PG1(I) = TrueObs(I);
PI1(ZN == 0) = zeros;
PG1(ZN == 0) = zeros;

[Nzeros,Nones]=number_neighbours(Z1,4);
%[Nzeros,Nones]=number_neighbours(Z1,8);

% Compute matrices Alpha and Theta which will be used to compute PI1
% and PG1 without loops
Alpha = alpha(1)*ones(size(E));
Alpha(find(E)) = alpha(2);
% Theta = theta(1)*ones(size(A));
% Theta(find(A)) = theta(2);
    
%Simulation of a realisation  Z1 of PG(Z | O, E, A, alpha, beta, Theta)
% There will be q passes of updates of all Z1(i,j)
for t=1:q
    %clf;
    %spy(Z1);
 
    
    
    % Update each Z1(i,j) in turn
    for i=1:size(Z,1)
        for j=1:size(Z,2)
            if (ZN(i,j) == 1) && (TrueObs(i,j) == -1)
                             
                % Computation of the normalizing constant:
                W1(i,j)=exp(beta*Nzeros(i,j))+exp(Alpha(i,j)+beta*Nones(i,j));
                              
                %Computation ofP(Z_ij = 1 | Z1, alpha, beta, e_ij);
                PI1(i,j)=exp(Alpha(i,j)+beta*Nones(i,j))/W1(i,j);
            
                % Computation of P_G^ij(Z_ij = 1|Z1,e_ij, a_ij,o_ij,alpha,beta,Theta)
                if O(i,j) == 1
                    PG1(i,j)=(PI1(i,j)*(1-Theta(2)))/((1-Theta(2))*PI1(i,j) + Theta(1)*(1-PI1(i,j)));
                else
                    PG1(i,j)=(PI1(i,j)*Theta(2))/((1-Theta(1))*(1-PI1(i,j)) + Theta(2)*PI1(i,j));
                end
 
                % Zloc(i,j) sample
                p=rand();
                if (p<PG1(i,j))
                   Zloc(i,j)=1;
                else
                   Zloc(i,j)=0;
                end
                
                % Update of Nzeros and Nones and Z1
                if (Zloc(i,j)~=Z1(i,j))
                    k=Zloc(i,j)-Z1(i,j);
                    if (i>1)
                        Nones(i-1,j)=Nones(i-1,j)+k;
                        Nzeros(i-1,j)=Nzeros(i-1,j)-k;
                    end
                    if (i<size(Z,1))
                        Nones(i+1,j)=Nones(i+1,j)+k;
                        Nzeros(i+1,j)=Nzeros(i+1,j)-k;
                    end
                    if (j>1)
                        Nones(i,j-1)=Nones(i,j-1)+k;
                        Nzeros(i,j-1)=Nzeros(i,j-1)-k;
                    end
                    if (j<size(Z,2))
                        Nones(i,j+1)=Nones(i,j+1)+k;
                        Nzeros(i,j+1)=Nzeros(i,j+1)-k;
                    end
                    Z1(i,j)=Zloc(i,j);
                end
            end
        end
    end
end

%computation of P_G^ij(Z_ij = 1|Z1,e_ij, a_ij,o_ij,alpha,beta,Theta)
for i=1:size(Z,1)
        for j=1:size(Z,2)
            if (ZN(i,j) == 1) && (TrueObs(i,j) == -1)
                             
                % Computation of the normalizing constant:
                W1(i,j)=exp(beta*Nzeros(i,j))+exp(Alpha(i,j)+beta*Nones(i,j));
                              
                %Computation ofP(Z_ij = 1 | Z1, alpha, beta, e_ij);
                PI1(i,j)=exp(Alpha(i,j)+beta*Nones(i,j))/W1(i,j);
            
                % Computation of P_G^ij(Z_ij = 1|Z1,e_ij, a_ij,o_ij,alpha,beta,Theta)
                if O(i,j) == 1
                    PG1(i,j)=(PI1(i,j)*(1-Theta(2)))/((1-Theta(2))*PI1(i,j) + Theta(1)*(1-PI1(i,j)));
                else
                    PG1(i,j)=(PI1(i,j)*Theta(2))/((1-Theta(1))*(1-PI1(i,j)) + Theta(2)*PI1(i,j));
                end

            end;
        end;
end;
 


    
