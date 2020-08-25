% Author: Ernst Roos
% Last Edited: 24-07-2020

function [ Q_val, c_val ] = SOCP_MVE_exact( vertices )

[n, L] = size(vertices);

%% Variables
Q = sdpvar(L);
c = sdpvar(L, 1);

%% Objective
obj = - log(det(Q));

%% Constraints
constraints = [];
for i = 1 : n
   constraints = [constraints;
       norm(Q * vertices(i,:)' + c, 2) <= 1];
end

%% Solve
approx_settings = sdpsettings('solver', 'mosek', 'verbose', 0);
optimize(constraints, obj, approx_settings);

Q_val = double(Q);
c_val = double(c);

end

