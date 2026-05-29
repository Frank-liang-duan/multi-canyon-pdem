function [A1data, aa1PDF, aa2CDF] = func_pdem(alph, eta, mu_R, cov_R, mu_D, cov_D, N_side, a0, xx_range)
% FUNC_PDEM  Core PDEM analysis for one parameter combination.
%   [A1data, aa1PDF, aa2CDF] = func_pdem(alph, eta, mu_R, cov_R, mu_D, cov_D, N_side, a0, xx_range)
%   performs representative point selection, deterministic solving, KDE,
%   and reliability analysis.
%
%   Inputs:
%       alph    - incident angle (degrees)
%       eta     - dimensionless frequency
%       mu_R, cov_R - mean and CoV of radius
%       mu_D, cov_D - mean and CoV of spacing
%       N_side  - number of valleys on each side (total = 2*N_side+1)
%       a0      - reference half-width (default 1)
%       xx_range- [xmin, dx, xmax] for observation points (default [-8,0.02,8])
%
%   Outputs:
%       A1data  - matrix [xx, mean, 95% quantile, mean-std, mean+std, CoV, Pf]
%       aa1PDF  - matrix [z, PDF at 3 key points]
%       aa2CDF  - matrix [z, CDF at 3 key points]

if nargin < 8, a0 = 1; end
if nargin < 9, xx_range = [-1, 0.02, 1]; end

Ja = -N_side:N_side;
NJ = length(Ja);
alph_rad = alph * pi/180;
cs = 1;
w = pi * eta * cs / a0;
k = w / cs;

% Truncation order
if eta < 1
    NN = 12;
elseif eta < 2
    NN = 16;
elseif eta < 3
    NN = 20;
else
    NN = 24;
end

if NJ > 3 
    xx = (xx_range(1) : xx_range(2) : xx_range(3))';
else
    xx = (-8 : 0.02 : 8)';
end

Numx = length(xx);

% Log-normal parameters
sigma_R = sqrt(log(1 + cov_R^2));
mu_ln_R = log(mu_R) - 0.5 * sigma_R^2;
sigma_D = sqrt(log(1 + cov_D^2));
mu_ln_D = log(mu_D) - 0.5 * sigma_D^2;

% Generate optimized point set
N_sel = 2500;
[R0, D0] = generate_canyon_samples(mu_R, cov_R, mu_D, cov_D, N_side, N_sel, a0);
[R_samples, D_samples] = post_optimize_samples(R0, D0, a0, mu_ln_R, sigma_R, mu_ln_D, sigma_D);

N_sel = size(R_samples, 1);
P_q = ones(N_sel,1) / N_sel;
fprintf('Effective representative points: %d\n', N_sel);

% Deterministic solutions
Y_all = zeros(Numx, N_sel);
for q = 1:N_sel
    R_left  = R_samples(q, 1:N_side);
    R_right = R_samples(q, N_side+1:end);
    D_left  = D_samples(q, 1:N_side);
    D_right = D_samples(q, N_side+1:end);
    Raj = [fliplr(R_left), a0, R_right];
    d_left_cum = fliplr(cumsum(fliplr(D_left)));
    d_right_cum = cumsum(D_right);
    d_a = [d_left_cum, 0, d_right_cum];
    
    if NJ > 3
        Y_all(:, q) = compute_displacement_vectorN(Ja, Raj, d_a, alph_rad, k, NN);
    else
        Y_all(:, q) = compute_displacement_vector3(Ja, Raj, d_a, alph_rad, k, NN);
    end
        
end

% KDE post-processing
u_max = 1.2 * max(Y_all(:));
n_z = 500;
z = linspace(0, u_max, n_z)';
PDF_all = zeros(Numx, n_z);
CDF_all = zeros(Numx, n_z);
for ix = 1:Numx
    [pdf, ~] = ksdensity(Y_all(ix,:)', z, 'Weights', P_q);
    pdf = pdf / trapz(z, pdf);
    PDF_all(ix,:) = pdf;
    CDF_all(ix,:) = cumtrapz(z, pdf);
end

% Reliability analysis
u_limit = 3.0;
[~, idx_limit] = min(abs(z - u_limit));
Pf = 1 - CDF_all(:, idx_limit);

U_mean = zeros(Numx,1);
U_std  = zeros(Numx,1);
for ix = 1:Numx
    U_mean(ix) = trapz(z, z .* PDF_all(ix,:)');
    U_std(ix)  = sqrt(trapz(z, (z - U_mean(ix)).^2 .* PDF_all(ix,:)'));
end
CoV_u = U_std ./ U_mean;

U_95 = zeros(Numx,1);
for ix = 1:Numx
    cdf = CDF_all(ix,:);
    idx = find(cdf >= 0.95, 1);
    if isempty(idx)
        U_95(ix) = z(end);
    elseif idx == 1
        U_95(ix) = z(1);
    else
        x1 = cdf(idx-1); x2 = cdf(idx);
        y1 = z(idx-1);   y2 = z(idx);
        U_95(ix) = y1 + (0.95 - x1) * (y2 - y1) / (x2 - x1);
    end
end

% Key point indices (center and ±0.3a)
idx_center = round(Numx/2);
idx_left   = round((Numx/2) - 0.3/xx_range(2));
idx_right  = round((Numx/2) + 0.3/xx_range(2));
A1data = [xx, U_mean, U_95, U_mean-U_std, U_mean+U_std, CoV_u, Pf];
aa1PDF = [z, PDF_all([idx_left, idx_center, idx_right], :)'];
aa2CDF = [z, CDF_all([idx_left, idx_center, idx_right], :)'];
end