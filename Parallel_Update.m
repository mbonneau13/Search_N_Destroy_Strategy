function [mess_out,epsilon]=Parallel_Update(mess_in,Psi,beta,nei_type,ZN);
% Single step update of all messages m_ij(xj).

% INPUT :
%        - mess_in(d,i,j,x) : messages coming from direction d to cell (i,j).
%          with d = {1,2,3,4} = {top,right,bottom,left} if nei_type = 4
%          with d = {1,2,3,4, 5, 6, 7 ,8} = {top,right,bottom,left, NE, SE,
%          SO, NO} if nei_type = 8
%        - Psi(i,j,x) : the value of psi_ij(x-1).
%        - beta : correlation coefficient of the MRF, used to compute
%        psi_(ij,kl) = exp(beta*eq(x_ij,x_kl))
%        - nei_type = 4 or 8 (resp. 4 or 8 closest neighbours)
% OUTPUT :
%        - mess_out(d,i,j,x) : new messages.
%        - epsilon : maximum absolut difference between mess_in and mess_out

% Initialisations
epsilon=0;
n=size(Psi,1);
m=size(Psi,2);
mess_out = mess_in;


%%%%% 4 neighbours %%%%%%%%%%%%%%%
if (nei_type == 4)

    % Computation of the "fixed" border messages:
    mess_out(1,1,:,:)=1;
    mess_out(2,:,m,:)=1;
    mess_out(3,n,:,:)=1;
    mess_out(4,:,1,:)=1;


    % Message updates:

    % Messages from top
    for i=2:n % Messages from top to line 1 are not updated
    for j=1:m
        for xij=1:2
            if ZN(i,j) == 1
                mess_out(1,i,j,xij) = Psi(i-1,j,1)*exp(beta*(xij==1))*mess_in(2,i-1,j,1)*mess_in(1,i-1,j,1)*mess_in(4,i-1,j,1)+Psi(i-1,j,2)*exp(beta*(xij==2))*mess_in(2,i-1,j,2)*mess_in(1,i-1,j,2)*mess_in(4,i-1,j,2);
            else
                mess_out(1,i,j,xij) = 1;
            end
        end
        mess_out(1,i,j,:)=mess_out(1,i,j,:)/sum(mess_out(1,i,j,:));
    end
    end


    epsilon=max(max(abs(mess_in(1,:,:,1)-mess_out(1,:,:,1))));

    % Messages from right
    for i=1:n 
    for j=1:m-1 % Messages from right to column M are not updated
        for xij=1:2
            if ZN(i,j) == 1
                mess_out(2,i,j,xij) = Psi(i,j+1,1)*exp(beta*(xij==1))*mess_in(1,i,j+1,1)*mess_in(3,i,j+1,1)*mess_in(2,i,j+1,1)+Psi(i,j+1,2)*exp(beta*(xij==2))*mess_in(1,i,j+1,2)*mess_in(3,i,j+1,2)*mess_in(2,i,j+1,2);
            else
                mess_out(2,i,j,xij) = 1;
            end
        end
        mess_out(2,i,j,:)=mess_out(2,i,j,:)/sum(mess_out(2,i,j,:));
    end
    end


    epsilon=max(epsilon,max(max(abs(mess_in(2,:,:,1)-mess_out(2,:,:,1)))));

    % Messages from bottom
    for i=1:n-1 % Messages from bottom to line N are not updated
    for j=1:m 
        for xij=1:2
            if ZN(i,j) == 1
                mess_out(3,i,j,xij) = Psi(i+1,j,1)*exp(beta*(xij==1))*mess_in(3,i+1,j,1)*mess_in(2,i+1,j,1)*mess_in(4,i+1,j,1)+Psi(i+1,j,2)*exp(beta*(xij==2))*mess_in(3,i+1,j,2)*mess_in(2,i+1,j,2)*mess_in(4,i+1,j,2);
            else
                mess_out(3,i,j,xij) = 1;
            end
        end
        mess_out(3,i,j,:)=mess_out(3,i,j,:)/sum(mess_out(3,i,j,:));
    end
    end


    epsilon=max(epsilon,max(max(abs(mess_in(3,:,:,1)-mess_out(3,:,:,1)))));

    % Messages from left
    for i=1:n 
    for j=2:m % Messages from left to column 1 are not updated
        for xij=1:2
            if ZN(i,j) == 1
                mess_out(4,i,j,xij) = Psi(i,j-1,1)*exp(beta*(xij==1))*mess_in(1,i,j-1,1)*mess_in(4,i,j-1,1)*mess_in(3,i,j-1,1)+Psi(i,j-1,2)*exp(beta*(xij==2))*mess_in(1,i,j-1,2)*mess_in(4,i,j-1,2)*mess_in(3,i,j-1,2);
            else
                mess_out(4,i,j,xij) = 1;
            end
        end
        mess_out(4,i,j,:)=mess_out(4,i,j,:)/sum(mess_out(4,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(4,:,:,1)-mess_out(4,:,:,1)))));

%%%%% 8 neighbours %%%%%%%%%%%%%%%
else
    
    % Computation of the "fixed" border messages:
    mess_out([1 8 5],1,:,:)=1;
    mess_out([2 5 6],:,m,:)=1;
    mess_out([3 6 7],n,:,:)=1;
    mess_out([4 7 8],:,1,:)=1;


    % Message updates:

    % Messages from top
    for i=2:n % Messages from top to line 1 are not updated
    for j=1:m
        for xij=1:2
            mess_out(1,i,j,xij) = Psi(i-1,j,1)*exp(beta*(xij==1))*mess_in(2,i-1,j,1)*mess_in(1,i-1,j,1)*mess_in(4,i-1,j,1)*mess_in(5,i-1,j,1)*mess_in(6,i-1,j,1)*mess_in(7,i-1,j,1)*mess_in(8, i-1,j,1)+Psi(i-1,j,2)*exp(beta*(xij==2))*mess_in(2,i-1,j,2)*mess_in(1,i-1,j,2)*mess_in(4,i-1,j,2)*mess_in(5, i-1,j,2)*mess_in(6, i-1,j,2)*mess_in(7, i-1,j,2)*mess_in(8, i-1,j,2);
        end
        mess_out(1,i,j,:)=mess_out(1,i,j,:)/sum(mess_out(1,i,j,:));
    end
    end


    epsilon=max(max(abs(mess_in(1,:,:,1)-mess_out(1,:,:,1))));

    % Messages from right
    for i=1:n 
    for j=1:m-1 % Messages from right to column M are not updated
        for xij=1:2
            mess_out(2,i,j,xij) = Psi(i,j+1,1)*exp(beta*(xij==1))*mess_in(1,i,j+1,1)*mess_in(3,i,j+1,1)*mess_in(2,i,j+1,1)*mess_in(5,i,j+1,1)*mess_in(6,i,j+1,1)*mess_in(7,i,j+1,1)*mess_in(8,i,j+1,1) +Psi(i,j+1,2)*exp(beta*(xij==2))*mess_in(1,i,j+1,2)*mess_in(3,i,j+1,2)*mess_in(2,i,j+1,2)*mess_in(5,i,j+1,2)*mess_in(6,i,j+1,2)*mess_in(7,i,j+1,2)*mess_in(8,i,j+1,2);
        end
        mess_out(2,i,j,:)=mess_out(2,i,j,:)/sum(mess_out(2,i,j,:));
    end
    end


    epsilon=max(epsilon,max(max(abs(mess_in(2,:,:,1)-mess_out(2,:,:,1)))));

    % Messages from bottom
    for i=1:n-1 % Messages from bottom to line N are not updated
    for j=1:m 
        for xij=1:2
            mess_out(3,i,j,xij) = Psi(i+1,j,1)*exp(beta*(xij==1))*mess_in(3,i+1,j,1)*mess_in(2,i+1,j,1)*mess_in(4,i+1,j,1)*mess_in(5,i+1,j,1)*mess_in(6,i+1,j,1)*mess_in(7,i+1,j,1)*mess_in(8,i+1,j,1) +Psi(i+1,j,2)*exp(beta*(xij==2))*mess_in(3,i+1,j,2)*mess_in(2,i+1,j,2)*mess_in(4,i+1,j,2)*mess_in(5,i+1,j,2)*mess_in(6,i+1,j,2)*mess_in(7,i+1,j,2)*mess_in(8,i+1,j,2);
        end
        mess_out(3,i,j,:)=mess_out(3,i,j,:)/sum(mess_out(3,i,j,:));
    end
    end


    epsilon=max(epsilon,max(max(abs(mess_in(3,:,:,1)-mess_out(3,:,:,1)))));

    % Messages from left
    for i=1:n 
    for j=2:m % Messages from left to column 1 are not updated
        for xij=1:2
            mess_out(4,i,j,xij) = Psi(i,j-1,1)*exp(beta*(xij==1))*mess_in(1,i,j-1,1)*mess_in(4,i,j-1,1)*mess_in(3,i,j-1,1)*mess_in(5,i,j-1,1)*mess_in(6,i,j-1,1)*mess_in(7,i,j-1,1)*mess_in(8,i,j-1,1)+Psi(i,j-1,2)*exp(beta*(xij==2))*mess_in(1,i,j-1,2)*mess_in(4,i,j-1,2)*mess_in(3,i,j-1,2)*mess_in(5,i,j-1,2)*mess_in(6,i,j-1,2)*mess_in(7,i,j-1,2)*mess_in(8,i,j-1,2);
        end
        mess_out(4,i,j,:)=mess_out(4,i,j,:)/sum(mess_out(4,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(4,:,:,1)-mess_out(4,:,:,1)))));

    % Messages from NE
    for i=2:n % Messages from NE to line 1 are not updated
    for j=1:m - 1 % Messages from NE to column M are not updated
        for xij=1:2
            mess_out(5,i,j,xij) = Psi(i-1,j+1,1)*exp(beta*(xij==1))*mess_in(1,i-1,j+1,1)*mess_in(2,i-1,j+1,1)*mess_in(3,i-1,j+1,1)*mess_in(4,i-1,j+1,1)*mess_in(5,i-1,j+1,1)*mess_in(6,i-1,j+1,1)*mess_in(8,i-1,j+1,1)  + Psi(i-1,j+1,2)*exp(beta*(xij==2))*mess_in(1,i-1,j+1,2)*mess_in(2,i-1,j+1,2)*mess_in(3,i-1,j+1,2)*mess_in(4,i-1,j+1,2)*mess_in(5,i-1,j+1,2)*mess_in(6,i-1,j+1,2)*mess_in(8,i-1,j+1,2);
        end
        mess_out(5,i,j,:)=mess_out(5,i,j,:)/sum(mess_out(5,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(5,:,:,1)-mess_out(5,:,:,1)))));
    
   % Messages from SE
    for i=1:n - 1 % Messages from SE to line N are not updated
    for j=1:m - 1 % Messages from SE to column M are not updated
        for xij=1:2
            mess_out(6,i,j,xij) = Psi(i+1,j+1,1)*exp(beta*(xij==1))*mess_in(1,i+1,j+1,1)*mess_in(2,i+1,j+1,1)*mess_in(3,i+1,j+1,1)*mess_in(4,i+1,j+1,1)*mess_in(5,i+1,j+1,1)*mess_in(6,i+1,j+1,1)*mess_in(7,i+1,j+1,1)+ Psi(i+1,j+1,2)*exp(beta*(xij==2))*mess_in(1,i+1,j+1,2)*mess_in(2,i+1,j+1,2)*mess_in(3,i+1,j+1,2)*mess_in(4,i+1,j+1,2)*mess_in(5,i+1,j+1,2)*mess_in(6,i+1,j+1,2)*mess_in(7,i+1,j+1,2);
        end
        mess_out(6,i,j,:)=mess_out(6,i,j,:)/sum(mess_out(6,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(6,:,:,1)-mess_out(6,:,:,1)))));
     
    
      % Messages from SO
    for i=1:n - 1 % Messages from SO to line N are not updated
    for j=2:m  % Messages from SO to column 1 are not updated
        for xij=1:2
            mess_out(7,i,j,xij) = Psi(i+1,j-1,1)*exp(beta*(xij==1))*mess_in(1,i+1,j-1,1)*mess_in(2,i+1,j-1,1)*mess_in(3,i+1,j-1,1)*mess_in(4,i+1,j-1,1)*mess_in(6,i+1,j-1,1)*mess_in(7,i+1,j-1,1)*mess_in(8,i+1,j-1,1)+ Psi(i+1,j-1,2)*exp(beta*(xij==2))*mess_in(1,i+1,j-1,2)*mess_in(2,i+1,j-1,2)*mess_in(3,i+1,j-1,2)*mess_in(4,i+1,j-1,2)*mess_in(6,i+1,j-1,2)*mess_in(7,i+1,j-1,2)*mess_in(8,i+1,j-1,2) ;
        end
        mess_out(7,i,j,:)=mess_out(7,i,j,:)/sum(mess_out(7,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(7,:,:,1)-mess_out(7,:,:,1))))); 
    
          % Messages from NO
    for i=2:n  % Messages from NO to line 1 are not updated
    for j=2:m  % Messages from NO to column 1 are not updated
        for xij=1:2
            mess_out(8,i,j,xij) = Psi(i-1,j-1,1)*exp(beta*(xij==1))*mess_in(1,i-1,j-1,1)*mess_in(2,i-1,j-1,1)*mess_in(3,i-1,j-1,1)*mess_in(4,i-1,j-1,1)*mess_in(5,i-1,j-1,1)*mess_in(7,i-1,j-1,1)*mess_in(8,i-1,j-1,1) + Psi(i-1,j-1,2)*exp(beta*(xij==2))*mess_in(1,i-1,j-1,2)*mess_in(2,i-1,j-1,2)*mess_in(3,i-1,j-1,2)*mess_in(4,i-1,j-1,2)*mess_in(5,i-1,j-1,2)*mess_in(7,i-1,j-1,2)*mess_in(8,i-1,j-1,2) ;
        end
        mess_out(8,i,j,:)=mess_out(8,i,j,:)/sum(mess_out(8,i,j,:));
    end
    end
 

    epsilon=max(epsilon,max(max(abs(mess_in(8,:,:,1)-mess_out(8,:,:,1)))));
end;




 

