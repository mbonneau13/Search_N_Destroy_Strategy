function Run_FunctionGUI(O,ZN,Search,Destroy,Budget,alpha1,beta1,Theta1,nTest,max_iter,Nit,ligne,colonne)
global K VV






%%% Simulate true Map
K = 1;
VV = [0 1];
% Nit = 20;

%%% Parameter estimation
E = zeros(ligne,colonne);
A = zeros(ligne,colonne);
TrueObs = -ones(ligne,colonne);

% alpha1 = [-0.5 -0.5];
% beta1 = 0.5;
% Theta1 = [0.5 0.5];
% nTest = 5;
MRF = zeros(nTest,2);
Noize = zeros(nTest,2);

for i = 1:nTest
    
    close all
    waitbar(i/10,['Parameters Estimation ' int2str(i) '/' int2str(nTest) ' ...'])    
    [alpha1,beta1,Theta1,Z1,PI1,PG1,Psi_vec,t,precision] = EM_test(A,O,E,alpha1,beta1,Theta1,max_iter,Nit,1,1,TrueObs,ZN);
    MRF(i,1) = alpha1(1);
    MRF(i,2) = beta1;
    Noize(i,:) = Theta1;
    disp(['Current MRF Parameters: \alpha  = ' num2str(MRF(i,1)) ', \beta = ' num2str(MRF(i,2))])
    disp(['Current noize Parameters: FP  = ' num2str(Noize(i,1)) ', FN = ' num2str(Noize(i,2))])
end


alpha = zeros(ligne,colonne,2);
%     alpha(:,:,2) = mean(MRF(:,1));
%     beta = mean(MRF(:,2));
%     Theta = mean(Noize,1);
alpha(:,:,2) = MRF(end,1);
beta = MRF(end,2);
Theta = Noize(end,:);

% alpha = zeros(ligne,colonne,2);
% alpha(:,:,2) =0.040198;
% beta = 0.85454;
% Theta = [0 0.34365];


disp('==== Estimated Parameters ====')
disp(['\alpha = ' num2str(MRF(end,1))])
disp(['\beta = ' num2str(MRF(end,2))])
disp(['\Theta_{FP} = ' num2str(Noize(end,1))])
disp(['\Theta_{FN} = ' num2str(Noize(end,2))])
[marginals] = Compute_marginals(A,O,E,[alpha alpha],beta,Theta,1,ZN,TrueObs);
close all
waitbar(0.8,'Computing the optimal strategy ...')

I_ZN = find(ZN == 1);

%%% LIP strategy
E_Cost = marginals.*(Search + Destroy) + (ones-marginals).*Search;
%     f = reshape(-marginals,1,ligne*colonne);
f = -1000.*marginals(I_ZN)';
nvars = length(f);
A = E_Cost(I_ZN)';
lb = zeros(nvars,1);
ub = ones(nvars,1);
intcon = 1:nvars;
options = optimoptions('intlinprog','Display','off');

b = Budget;
%     Aeq = double(reshape(~ZN,1,ligne*colonne));
%     beq = 0;
Aeq = [];
beq = [];
[x , ~ , ~ , ~] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);


x(x<0.5) = 0;
x(x>=0.5) = 1;
x_LIP = I_ZN(x == 1);

disp(['Estimated number of searched and occupied cells = ' num2str(sum(marginals(x_LIP)))])
close all

Border = [];
for y = 1:numel(ZN)
    
    if ZN(y) == 0
        x = [(y + 1) , (y - 1) , (y + size(ZN,1)) , (y - size(ZN,1))];
        x = x(x > 0);
        x = x(x <= numel(ZN));
        u = unique(ZN(x));
        if (numel(u) > 1)
            Border = [Border , y];
        end
    end
end
I = find(ZN == 0);
figure('Name','Observation')
S = O;
S(I) = 2;
imagesc(S)
axis image
axis off
colormap([1 1 1;1 1 1;0 0 0;.5 .5 .5])

figure('Name','Cost')

S = Search  + Destroy;
S(S < 0) = zeros;
imagesc(S)
axis image
axis off
colormap(jet)
colorbar

figure('Name','Marginals')
S = marginals;
S(I) = 0;
imagesc(S)
axis image
axis off
colormap(rot90(colormap(gray(200)),2))

%% Show LIP
S = ZN;
S(x_LIP) = 2;
figure('Name','Where to search')
imagesc(S)
axis image
axis off
colormap([.9 .9 .9;0.6 0.6 0.6;1 0 0])

load('Out.mat')
csvwrite([Out '.csv'],x_LIP)

