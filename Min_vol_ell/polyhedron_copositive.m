% Author : Areesh Mittal

function [Q,c,obj,out] = polyhedron_copositive(D, d)

% Inputs S, t represent polytope P = {x : Sx <= t}
% Outputs A,b represent approximation to minimum volume ellipsoid 
% E = {x : ||Ax+b|| <= 1} of the polytope P
% Output obj = 1 / det(A) which is proportional to volume of E


options = sdpsettings('verbose', 0, 'solver', 'mosek');

D = [D; -eye(size(D, 2))];
d = [d; zeros(size(D, 2), 1)];

[M,N] = size(D);
Sc = [-D, d];

Q = sdpvar(N);
c = sdpvar(N,1);
NN = sdpvar(M);
obj = -logdet(Q);

Ac = [Q, c];

H = -Sc'*NN*Sc;
H(end,end) = H(end,end) + 1; 

cons = [NN(:) >= 0; [H,Ac';Ac,eye(N)] >= 0];

out = optimize(cons, obj, options);

Q = double(Q); c = double(c);
obj = exp(double(obj));
end