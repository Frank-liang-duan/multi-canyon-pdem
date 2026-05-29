function GF = compute_GF_discrepancy(R, D, mu_ln_R, sigma_R, mu_ln_D, sigma_D)
% COMPUTE_GF_DISCREPANCY  Compute the generalized F-discrepancy.
%   GF = compute_GF_discrepancy(R, D, mu_ln_R, sigma_R, mu_ln_D, sigma_D)
%   returns the maximum marginal F-discrepancy across all dimensions.
%
%   Inputs:
%       R, D     - sample matrices (N x 2*N_side)
%       mu_ln_*, sigma_* - log-normal distribution parameters
%
%   Output:
%       GF       - GF-discrepancy value

N_side = size(R,2)/2;
DVals = zeros(1, 4*N_side);
for i = 1:2*N_side
    [FR, xR] = ecdf(R(:,i));
    FtheoR = logncdf(xR, mu_ln_R, sigma_R);
    DVals(2*i-1) = max(abs(FR - FtheoR));

    [FD, xD] = ecdf(D(:,i));
    FtheoD = logncdf(xD, mu_ln_D, sigma_D);
    DVals(2*i) = max(abs(FD - FtheoD));
end
GF = max(DVals);
end