% Author: Ernst Roos
% Last Edited: 24-07-2020

%% Initialization
warning('off', 'MATLAB:singularMatrix')
clearvars
clc
clf

% Setting that toggles whether redundant constraints are removed
to_remove_redundant = true;

% Setting that toggles whether a random instance is generated or one is loaded from memory
generate_random = true;

% Setting that toggles whether the exact solution is determined (can only
% be used in extremely small dimensions!)
todo_exact = true;

% Setting that toggles whether a lower bound is computed
todo_lb = true;

% TODO - comment
lb_once = true;

% Parameter in which the instance numbers to be loaded are specified
instance_numbers = [1];
N_trials = length(instance_numbers);

% Number of adjustable variables to be eliminated with Fourier-Motzkin
% Elimination
n_FME = 0;

% Dimensions of the problem (L is number of variables, M is number of
% polyhedral constraints)
L = 2; M = 2;

% Containers for results
errors = zeros(N_trials, 5);
time = zeros(N_trials, 5);

%% Main Loop over instances
for iter = instance_numbers
    fprintf('Processing Instance %d\n', iter)
    
    %%  Load Instance
    
    if generate_random
        [D, d] = RandomPolyhedronMittal(L, M);
        q = length(d);
    else
        load(strcat('Instances/', string(L), '_', string(M), '_', string(iter), '.mat'))
    end
    
    % Redundant constraints are removed if toggled
    if to_remove_redundant
        remove_indices = [];
        aux_opts = sdpsettings('verbose', 0);
        for i = 1 : q
            % Setup auxiliary LP to determine whether constraint is
            % redundant
            x = sdpvar(L, 1);
            
            aux_obj = D(i,:)*x;
            aux_constr = [D([1:i-1, i+1:q],:)*x <= d([1:i-1, i+1:q]); x >= 0];
            
            
            status = optimize(aux_constr, -aux_obj, aux_opts);
            
            if value(aux_obj) <= d(i) && status.problem == 0
                remove_indices = [remove_indices; i];
            end
        end
        
        % Remove redundant constraints
        D(remove_indices,:) = [];
        d(remove_indices) = [];
        q = length(d);
    end
    
    
    %% Find exact solution
    tic
    if todo_exact
        if exist('Q_ex', 'var') == 0
            [Q_ex, c_ex] = SOCP_MVE_exact(find_vertices(D, d));
        end
    else
        Q_ex = eye(L); c_ex = zeros(L, 1);
    end
    time(iter, 2) = toc;
    
    %% Find LDR solution
    tic
    if todo_lb % Toggled when lower bound is to be computed
        [Q_lin, c_lin, Q_lb, c_lb, time_lb] = SOCP_MVE_linear(D, d, n_FME, todo_lb);
    else
        [Q_lin, c_lin] = SOCP_MVE_linear(D, d, n_FME, false);
        time_lb = 0; Q_lb = eye(L); c_lb = zeros(L, 1);
    end
    time(iter, 1) = time_lb;
    time(iter, 3) = toc - time_lb;
    
    %% Find Full Quadratic Decision Rule solution
    tic
    [Q_full_quad, c_full_quad] = SOCP_MVE_full_quadratic(D, d, n_FME);
    time(iter, 4) = toc;
    
    %% Find Copositive approach solution
    tic
    [Q_CP, c_CP] = polyhedron_copositive(D, d);
    time(iter, 5) = toc;
    
    %% Calculate error
    
    vol = @ (Q) 1 / det(Q);
    
    % Calculate errors based on proper value (lower bound/exact)
    if todo_exact || lb_once
        errors(iter, :) = [100 * ( (vol(Q_lb)/vol(Q_ex))^(1/L) - 1), ...
            0, ...
            100 * ( (vol(Q_lin)/vol(Q_ex))^(1/L) - 1), ...
            100 * ( (vol(Q_full_quad)/vol(Q_ex))^(1/L) - 1), ...
            100 * ( (vol(Q_CP)/vol(Q_ex))^(1/L) - 1)];
    else
        errors(iter, :) = [100 * ( (vol(Q_lb)/vol(Q_lb))^(1/L) - 1), ...
            0, ...
            100 * ( (vol(Q_lin)/vol(Q_lb))^(1/L) - 1), ...
            100 * ( (vol(Q_full_quad)/vol(Q_lb))^(1/L) - 1), ...
            100 * ( (vol(Q_CP)/vol(Q_lb))^(1/L) - 1)];
    end
    
    % Ensure proper clearing of saved exact/lower bound solutions
    if iter ~= instance_numbers(end)
        clear Q_ex c_ex Q_lb c_lb
    end
end

%% Stats

indices = [3; 4; 5];
to_print = [mean(errors(:, indices), 1); mean(time(:, indices), 1)];

fprintf('& %2.1f\\%% & %2.1f\\%% & %2.1f\\%%\n', to_print(1,1), to_print(1,2), to_print(1,3))
fprintf('& %1.3f & %1.3f & %1.3f\n', to_print(2,1), to_print(2,2), to_print(2,3))

%% 2D Illustration
if L == 2
    figure(1)
    clf
    x = sdpvar(L, 1);
    hold on
    plot(norm(Q_lin*x + c_lin, 2) <= 1, [],'b')
    plot(norm(Q_full_quad*x + c_full_quad, 2) <= 1, [], 'm')
    plot(norm(Q_ex*x + c_ex, 2) <= 1, [], 'g')
    plot([D * x <= d; x >= 0], [], 'r')
    hold off
    axis([-0.5 1.5 -0.5 1.5])
end
