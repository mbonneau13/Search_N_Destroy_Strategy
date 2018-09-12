function [Nzeros,Nones]=number_neighbours(Z,nei_type);

%for a given 0/1 matrix Z, this function computes the number of neighbours
%in state 1 and in state 0 of each cell (i,j)
%neighbours of a cell (i,j) are cells (i,j+1), (i,j-1), (i+1,j), (i-1,j)
%except on boundaries

%INPUT 
%     - hidden variable Z
%     - nei_type = 4 or 8 (resp. 4 or 8 closest neighbours)
% OUTPUT
%     - Nzeros : matrix of same size as Z
%               Nzeros(i,j) is the number of zeros in the neighbours of
%               (i,j)
%     - Nones : matrix of same size as Z
%               Nones(i,j) is the number of ones in the neighbours of
%               (i,j)



if (nei_type == 4) 
    % Z is shifted up, down, right, left
    Zup = [Z(2:size(Z,1),:) ; zeros(1,size(Z,2))];
    Zdown = [zeros(1,size(Z,2)) ; Z(1:(size(Z,1)-1),:)];
    Zright = [zeros(size(Z,1),1) Z(:,1:(size(Z,2)-1))];
    Zleft = [Z(:,2:size(Z,2)) zeros(size(Z,1),1)];

    Nones = Zup+Zdown+Zright+Zleft;

    Nzeros = 4-Nones;
    Nzeros(1,:) = Nzeros(1,:)-1;
    Nzeros(size(Z,1),:) = Nzeros(size(Z,1),:)-1;
    Nzeros(:,1) = Nzeros(:,1)-1;
    Nzeros(:,size(Z,2)) = Nzeros(:,size(Z,2))-1;

else
    % Z is shifted up, down, right, left
    Zup = [Z(2:size(Z,1),:) ; zeros(1,size(Z,2))];
    Zdown = [zeros(1,size(Z,2)) ; Z(1:(size(Z,1)-1),:)];
    Zright = [zeros(size(Z,1),1) Z(:,1:(size(Z,2)-1))];
    Zleft = [Z(:,2:size(Z,2)) zeros(size(Z,1),1)];
    Zne = zeros(size(Z));
    Zne(2:size(Z,1), 1:size(Z,2) -1) = Z(1:size(Z,1) -1, 2:size(Z,2));
    Zse = zeros(size(Z));
    Zse(1:size(Z,1)-1, 1:size(Z,2) -1) = Z(2:size(Z,1) , 2:size(Z,2));
    Zso = zeros(size(Z));
    Zso(1:size(Z,1)-1, 2:size(Z,2) ) = Z(2:size(Z,1) , 1:size(Z,2)-1);
    Zno = zeros(size(Z));
    Zno(2:size(Z,1), 2:size(Z,2) ) = Z(1:size(Z,1) -1, 1:size(Z,2)-1);

    Nones = Zup+Zdown+Zright+Zleft + Zno + Zso + Zse + Zne;

    Nzeros = 8-Nones;
    Nzeros(1,:) = Nzeros(1,:)-3;
    Nzeros(size(Z,1),:) = Nzeros(size(Z,1),:)-3;
    Nzeros(:,1) = Nzeros(:,1)-3;
    Nzeros(:,size(Z,2)) = Nzeros(:,size(Z,2))-3;
    
    Nzeros(1,1) = Nzeros(1,1) +1;
    Nzeros(1,size(Z,2)) = Nzeros(1,size(Z,2)) +1;
    Nzeros(size(Z,1),1) = Nzeros(size(Z,1),1) +1;
    Nzeros(size(Z,1),size(Z,2)) = Nzeros(size(Z,1),size(Z,2)) +1;
    
    
    
end;



































