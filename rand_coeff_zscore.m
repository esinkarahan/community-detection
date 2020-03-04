function randcoeff = rand_coeff_zscore(A,B)
% A, B: partitions/graphs (n x n)
% calculate the Rand z-score to compare community structure of graphs
% community: mesoscopic groups of nodes with more internal connections than
% external connections
% Ref: Traud 2011, SIAM
%
% Esin Karahan, Jan, 2020, 

% number of nodes
n = length(A);
% number of classes in each partition
na = length(unique(A)); 
nb = length(unique(B));

% Contingency table
% % %
% The element CT_ij of the contingency table indicates the number of nodes 
% that are classified into the ith group of the first partition and 
% the jth group of the second partition.
CT = zeros(na,nb);

for i=1:na
    for j=1:nb
        %for a in cluster i, for b in cluster j
        CT(i,j) = sum((A==i).*(B==j));
    end
end

CTi = sum(CT,2); %row sum
CTj = sum(CT,1); %column sum

%total number of pairs
M  = n*(n-1)/2;  
M1 = comb(CTi);
M2 = comb(CTj);
% pair counting index
wab = 0;
for i=1:na
    for j = 1:nb
        if CT(i,j)>0
            wab = wab + (CT(i,j)*(CT(i,j)-1))/2;
        end
    end
end
C1 = n*(n^2 - 3*n - 2) - 8*(n+1)*M1 + 4*sum(CTi.^3);
C2 = n*(n^2 - 3*n - 2) - 8*(n+1)*M2 + 4*sum(CTj.^3);

sigmawab = M/16 - ((4*M1-2*M)^2*(4*M2-2*M)^2/(256*M^2)) + C1*C2/(16*n*(n-1)*(n-2)) ...
    +((4*M1 - 2*M)^2 - 4*C1 - 4*M)*((4*M2 - 2*M)^2 - 4*C2 - 4*M)/(64*n*(n-1)*(n-2)*(n-3));

randcoeff = 1/sqrt(sigmawab)*(wab - M1*M2/M);
end

function z = comb(n1)
% sum_i /n_i\ 
%       \ 2 /
z = 0;
for i=1:length(n1)
    if n1>0
        z = z+(n1(i)*(n1(i)-1))/2;
    end
end
end