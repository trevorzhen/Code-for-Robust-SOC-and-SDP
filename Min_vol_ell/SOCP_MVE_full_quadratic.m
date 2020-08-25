% Author: Ernst Roos
% Last Edited: 24-07-2020

function [ Q_val, c_val ] = SOCP_MVE_full_quadratic( D, d, n_FME )
%SOCP_MVE_QUADRATIC Decision rule with all quadratic & interaction terms

[q, L] = size(D);

%% FME

q = q - n_FME;

[lambda_coef, other_coef, b_FME] = fme_mve(D, d, n_FME);

n_c = size(lambda_coef, 1);

 
%% Variables
Q = sdpvar(L, L);
c = sdpvar(L, 1);

other_terms = [c, Q];

%% Decision Rule
u = sdpvar(q, 1, 'full');
V = sdpvar(q, L, 'full');
R = cell(q, 1);
for i = 1 : q
    R{i} = sdpvar(L, L);
end

%% Auxiliary
y = sdpvar(2, n_c, 'full');

%% Useful Expressions
A1 = [1 zeros(1, L+1); zeros(L+1, L+2)];
A2 = [zeros(1, L+2); zeros(L+1, 1), eye(L+1)];

R_hat = cell(n_c, 1);

for i = 1 : n_c
    R_hat{i} = zeros(L);
    for j = 1 : q
        R_hat{i} = R_hat{i} + lambda_coef(i,j) * R{j};
    end
end

%% Objective
obj = - log(det(Q));
%obj = 1;

%% Constraints

constraints = [];

for i = 1 : n_c
    constraints = [constraints;
        -sum(y(:,i)) <= b_FME(i);
        y(1,i) * A1 + y(2,i) * A2 <= - [lambda_coef(i,:) * u, (V'*lambda_coef(i,:)' + other_terms * other_coef(i,:)')'/2, 0;
        (V'*lambda_coef(i,:)' + other_terms * other_coef(i,:)')/2, R_hat{i}, zeros(L,1);
        zeros(1, L + 2)]];
end


%% Solve
approx_settings = sdpsettings('solver', 'mosek', 'verbose', 0);
optimize(constraints, obj, approx_settings);

Q_val = double(Q);
c_val = double(c);

yalmip('clear')
end