function T = nodal_assoc_matrix(C)
% C: npartition x nnodes
% Tij: Nodal association matrix, 
% how many times nodes i & j have been assigned to the same community 
%
% Esin Karahan, Jan, 2020, 

% number of nodes
[~,nn] = size(C);

T = zeros(nn);

for i=1:nn
    for j=1:nn
        %for node i in cluster a, for j in cluster a
        T(i,j) = sum(C(:,i)==C(:,j));
    end
end
