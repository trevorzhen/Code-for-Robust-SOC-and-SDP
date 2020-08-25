% Author: Ernst Roos
% Last Edited: 24-07-2020

function [lambda_coef, other_coef, b_FME] = fme_mve( D, d, n_FME )
%FME_MVE Constructs the appropriate constraint matrices to apply
%Fourier-Motzkin Elimination for the adjustable formulation of the Minimum
%Volume Ellipsoid Problem

[q, L] = size(D);

if n_FME > 0
    A_FME = [[d'; -D'; -eye(q)], [eye(L+1); zeros(q,L+1)]];
    b_FME = [1; zeros(L+q, 1)];
    
    [A_FME, b_FME] = core_fme(A_FME, b_FME, (q + L + 1 - n_FME));
    
    q = q - n_FME;
    
    lambda_coef = A_FME(:, 1:q);
    other_coef = A_FME(:, q+1:end);
else    
    lambda_coef = [d'; -D'; -eye(q)];
    other_coef = [eye(L+1); zeros(q,L+1)];
    
    b_FME = [1; zeros(L+q, 1)];
end


end

