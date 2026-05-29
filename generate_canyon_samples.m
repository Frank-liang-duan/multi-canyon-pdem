function [R_samples, D_samples] = generate_canyon_samples(mu_R, cov_R, mu_D, cov_D, N_side, N_sel, a0_fixed)
% GENERATE_VALLEY_SAMPLES  Generate initial representative points with geometry constraints.
%   [R_samples, D_samples] = generate_valley_samples(...)
%   uses a Sobol sequence and inverse CDF transform to draw candidate radii
%   and spacings, then filters out overlapping configurations.
%
%   Inputs:
%       mu_R, cov_R - mean and CoV of canyon radius (log-normal)
%       mu_D, cov_D - mean and CoV of canyon spacing (log-normal)
%       N_side      - number of valleys on each side of the central one
%       N_sel       - required number of valid samples
%       a0_fixed    - radius of the central valley (default 1)
%
%   Outputs:
%       R_samples - N_sel x 2*N_side matrix, [R_right, R_left]
%       D_samples - N_sel x 2*N_side matrix, spacing in the same order

if nargin < 7, a0_fixed = 1.0; end

sigma_ln_R = sqrt(log(1 + cov_R^2));
mu_ln_R    = log(mu_R) - 0.5 * sigma_ln_R^2;
sigma_ln_D = sqrt(log(1 + cov_D^2));
mu_ln_D    = log(mu_D) - 0.5 * sigma_ln_D^2;

total_vars = 4 * N_side;
try
    sob = sobolset(total_vars, 'Skip', 1e3, 'Leap', 1e2);
    u_all = net(sob, N_sel * 10);  % larger candidate pool to reduce filtering bias
catch
    error('Statistics and Machine Learning Toolbox (sobolset) is required.');
end

R_samples = zeros(N_sel, 2*N_side);
D_samples = zeros(N_sel, 2*N_side);
count_valid = 0;
idx_cand = 0;

while count_valid < N_sel && idx_cand < size(u_all,1)
    idx_cand = idx_cand + 1;
    u = u_all(idx_cand, :);
    half = 2 * N_side;
    u_R = u(1:half);
    u_D = u(half+1:end);

    R_cand = logninv(u_R, mu_ln_R, sigma_ln_R);
    D_cand = logninv(u_D, mu_ln_D, sigma_ln_D);

    R_left  = R_cand(1:N_side);
    R_right = R_cand(N_side+1:end);
    D_left  = D_cand(1:N_side);
    D_right = D_cand(N_side+1:end);

    [valid, R_ord, D_ord] = check_overlap(a0_fixed, R_left, R_right, D_left, D_right);
    if valid
        count_valid = count_valid + 1;
        R_samples(count_valid, :) = R_ord;
        D_samples(count_valid, :) = D_ord;
    end
end

if count_valid < N_sel
    warning('Only %d valid samples generated (requested %d).', count_valid, N_sel);
    R_samples = R_samples(1:count_valid, :);
    D_samples = D_samples(1:count_valid, :);
end
end