% Author: Ernst Roos
% Last Edited: 24-07-2020

function [ D, d ] = RandomPolyhedronMittal( L, M )
%RANDOMPOLYHEDRONMITTAL Constructs a random polyhedron for testing purposes
%as first explained by Mittal and Hanasusanto (2018) and later used in Zhen
%et al. (2020)

D = eye(L);
d = ones(L,1);

c = 0.5 * ones(L,1);

s = zeros(M, L);
r = zeros(M);

for j = 1 : M
    s(j, :) = rand(1, L);
    s(j, :) = s(j,:)/norm(s(j,:), 2);
    
    r(j) = norm(s(j,:), 1) * rand(1) - 0.5 * norm(s(j,:), 1);
    
    if r(j) > 0
        D = [D; s(j,:)];
        d = [d; r(j) + s(j,:) * c];
    else
        D = [D; -s(j,:)];
        d = [d; -r(j) - s(j,:) * c];
    end
end


end

