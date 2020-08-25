% Author: Ernst Roos
% Last Edited: 24-07-2020

function [ vertices ] = find_vertices( D ,d )
%FIND_VERTICES Finds the vertices of the polyhedron defined by 
% {x >= 0, D * x <= d}

[q, L] = size(D);

full_D = [D; -eye(L)]; full_d = [d; zeros(L,1)];
n_constr = length(full_d);

vertices = [];

combinations = VChooseK(1:n_constr, L);
n_comb = size(combinations, 1);

for i = 1 : n_comb

   D_p = full_D(combinations(i,:),:);

   vertex = D_p \ full_d(combinations(i,:));
   
   if full_D * vertex <= full_d + 1e-5
       vertices = [vertices; vertex'];
   end
end


end

