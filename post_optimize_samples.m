function [R_opt, D_opt] = post_optimize_samples(R, D, a0, mu_ln_R, sigma_R, mu_ln_D, sigma_D, varargin)
% POST_OPTIMIZE_SAMPLES  Iterative rearrangement under geometric constraints
%   to minimize the GF-discrepancy.
%
%   [R_opt, D_opt] = post_optimize_samples(R, D, a0, mu_ln_R, sigma_R, mu_ln_D, sigma_D)
%   performs quantile mapping while preserving non-overlap constraints.
%
%   Optional parameters:
%       'maxIter', 20  - maximum number of iterations
%       'tol', 1e-4    - convergence tolerance on GF-discrepancy change

p = inputParser;
addParameter(p, 'maxIter', 20);
addParameter(p, 'tol', 1e-4);
parse(p, varargin{:});
maxIter = p.Results.maxIter;
tol     = p.Results.tol;

N = size(R,1);
N_side = size(R,2)/2;
R_opt = R; D_opt = D;

GF_old = compute_GF_discrepancy(R_opt, D_opt, mu_ln_R, sigma_R, mu_ln_D, sigma_D);
fprintf('  Post-optimization start, initial GF = %.4f\n', GF_old);

for iter = 1:maxIter
    % Quantile mapping to ideal values
    R_ideal = zeros(size(R_opt));
    for j = 1:2*N_side
        [~, idx] = sort(R_opt(:,j));
        R_ideal(idx, j) = logninv(((1:N)'-0.5)/N, mu_ln_R, sigma_R);
    end
    D_ideal = zeros(size(D_opt));
    for j = 1:2*N_side
        [~, idx] = sort(D_opt(:,j));
        D_ideal(idx, j) = logninv(((1:N)'-0.5)/N, mu_ln_D, sigma_D);
    end

    % Check constraints pointwise
    R_new = R_opt; D_new = D_opt;
    accept = false(N,1);
    for i = 1:N
        R_left_ideal  = R_ideal(i, 1:N_side);
        R_right_ideal = R_ideal(i, N_side+1:end);
        D_left_ideal  = D_ideal(i, 1:N_side);
        D_right_ideal = D_ideal(i, N_side+1:end);

        [valid, ~, ~] = check_overlap(a0, fliplr(R_left_ideal), R_right_ideal, ...
                                      fliplr(D_left_ideal), D_right_ideal);
        if valid
            R_new(i,:) = R_ideal(i,:);
            D_new(i,:) = D_ideal(i,:);
            accept(i) = true;
        end
    end

    R_opt = R_new; D_opt = D_new;
    GF_new = compute_GF_discrepancy(R_opt, D_opt, mu_ln_R, sigma_R, mu_ln_D, sigma_D);
    fprintf('  Iteration %d: GF = %.6f, acceptance rate = %.1f%%\n', iter, GF_new, mean(accept)*100);

    if abs(GF_old - GF_new) < tol
        fprintf('  GF-discrepancy converged, stopping.\n');
        break;
    end
    GF_old = GF_new;
end
end