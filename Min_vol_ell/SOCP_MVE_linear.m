% Author: Ernst Roos
% Last Edited: 24-07-2020

function [ Q_val, c_val, Q_lb, c_lb, time_lb] = SOCP_MVE_linear( D, d, n_FME, to_lb)

[q, L] = size(D);

%% FME

q = q - n_FME;

[lambda_coef, other_coef, b_FME] = fme_mve(D, d, n_FME);

n_c = size(lambda_coef, 1);


%% Variables
Q = sdpvar(L);
c = sdpvar(L, 1);

other_terms = [c, Q];

%% Auxiliary Variables
u = sdpvar(q, 1, 'full');
V = sdpvar(q, L, 'full');

%% Objective
obj = - log(det(Q));

%% Constraints
constraints = [];

for i = 1 : n_c
    constraints = [constraints;
        lambda_coef(i,:)*u + norm(V'*lambda_coef(i,:)' + other_terms * other_coef(i,:)', 2) <= b_FME(i)];
end

%% Solve
approx_settings = sdpsettings('solver', 'mosek', 'verbose', 0);
optimize(constraints, obj, approx_settings);

Q_val = double(Q);
c_val = double(c);



%% Lower Bound Computation
if to_lb
    vol = @ (Q) 1 / det(Q);
    
    tic
    
    V_val = double(V);
    other_val = double(other_terms);
    
    %% Obtain Auxiliary Scenarios
    w_scen = ((V_val'*lambda_coef' + other_val * other_coef') ./ vecnorm(V_val'*lambda_coef' + other_val * other_coef'))';
    
    w_scen = w_scen(~any(isnan(w_scen), 2), :);
    w_scen = unique(round(w_scen, 4), 'rows');
    
    k_w = size(w_scen, 1);
 
        
    % Container
    x_scen = zeros(k_w, L);
    
    
    % Setup Gurobi model
    model.A = sparse(D);
    model.sense = '<';
    model.rhs = d;
    model.lb = zeros(L,1);
    
    params.outputflag = 0;
    params.method = 1;
    
    % Solve model for each scenario w
    for i = 1 : k_w
        model.obj = w_scen(i,:) * double(Q);
        result = gurobi(model, params);
        x_scen(i,:) = result.x;
    end
    
    % Filter unique scenarios    
    x_scen = unique(round(x_scen, 4), 'rows');
    k_x = size(x_scen, 1);
    
    % Use exact method based on only scenarios found
    [Q_lb, c_lb] = SOCP_MVE_exact(x_scen);
    
    
    time_lb = toc;
else
    Q_lb = zeros(L);
    c_lb = zeros(L,1);
    time_lb = 0;
end

yalmip('clear')
end

