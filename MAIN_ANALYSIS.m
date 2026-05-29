% MAIN_ANALYSIS  Full parametric study for PDEM valley-group analysis.
%   Performs four sequential analyses:
%     1. Incident angle
%     2. Frequency
%     3. Canyon number (N_side)
%     4. Random parameters (mean & CoV of radius/spacing)
%   Results are saved in subdirectories under the current folder.

clear; clc; close all;

% Add core functions to path (adjust relative path if necessary)
if exist('core', 'dir')
    addpath('core');
end

%% ==================== Common baseline parameters ====================
base_alph = 45;            % incident angle [deg]
base_eta  = 1;             % dimensionless frequency
base_mu_R = 1.0;           % mean radius
base_cov_R = 0.25;         % CoV of radius
base_mu_D = 4.0;           % mean spacing
base_cov_D = 0.25;         % CoV of spacing
base_N_side = 1;           % number of valleys on each side (3 canyons total)
a0 = 1;                    % reference half‑width
xx_range = [-8, 0.02, 8];  % observation window [min, step, max]

%% ================================================================
%  1. INCIDENT ANGLE ANALYSIS
%  ================================================================
disp('========== 1. Incident Angle Analysis ==========');
Lalph = [5, 30, 60, 90];
save_dir = 'ParamAnalysis_Angle';
if ~isfolder(save_dir), mkdir(save_dir); end

tt = tic;
for i_alph = 1:length(Lalph)
    alph = Lalph(i_alph);
    fprintf('  Angle = %d deg\n', alph);
    [A1data, aa1PDF, aa2CDF] = func_pdem(alph, base_eta, base_mu_R, base_cov_R, ...
                                         base_mu_D, base_cov_D, base_N_side, a0, xx_range);
    fname = sprintf('alpha%02d_eta%d.mat', alph, base_eta);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'alph', 'base_eta');
end
fprintf('  Angle analysis finished: %.1f sec\n\n', toc(tt));

%% ================================================================
%  2. FREQUENCY ANALYSIS
%  ================================================================
disp('========== 2. Frequency Analysis ==========');
Leta = [0.5, 1, 2, 4];
save_dir = 'ParamAnalysis_Freq';
if ~isfolder(save_dir), mkdir(save_dir); end

tt = tic;
for i_eta = 1:length(Leta)
    eta = Leta(i_eta);
    fprintf('  eta = %.1f\n', eta);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, eta, base_mu_R, base_cov_R, ...
                                         base_mu_D, base_cov_D, base_N_side, a0, xx_range);
    fname = sprintf('alpha%d_eta%.1f.mat', base_alph, eta);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'base_alph', 'eta');
end
fprintf('  Frequency analysis finished: %.1f sec\n\n', toc(tt));

%% ================================================================
%  3. CANYON NUMBER ANALYSIS (N_side)
%  ================================================================
disp('========== 3. Canyon Number Analysis ==========');
N_side_list = [1, 2, 4, 7, 11];  % corresponds to 3,5,7,9,11 canyons
save_dir = 'ParamAnalysis_Number';
if ~isfolder(save_dir), mkdir(save_dir); end

% For canyon number analysis, restrict observation to the central valley
xx_range_center = [-1, 0.02, 1];   % narrower window around center

tt = tic;
for i_n = 1:length(N_side_list)
    N_side = N_side_list(i_n);
    fprintf('  N_side = %d (total %d canyons)\n', N_side, 2*N_side+1);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, base_eta, base_mu_R, base_cov_R, ...
                                         base_mu_D, base_cov_D, N_side, a0, xx_range_center);
    fname = sprintf('Nside%d.mat', N_side);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'N_side');
end
fprintf('  Canyon number analysis finished: %.1f sec\n\n', toc(tt));

%% ================================================================
%  4. RANDOM PARAMETER ANALYSIS
%     (mean & CoV of radius and spacing)
%  ================================================================
disp('========== 4. Random Parameter Analysis ==========');
save_dir = 'ParamAnalysis_Random';
if ~isfolder(save_dir), mkdir(save_dir); end

% Parameter grids
mu_R_list  = [0.5, 1.0, 1.5, 2.0];
cov_R_list = [0.1, 0.2, 0.3, 0.4];
mu_D_list  = [3.0, 4.0, 5.0, 6.0];
cov_D_list = [0.1, 0.2, 0.3, 0.4];

tt = tic;

% --- 4a. Mean radius mu_R ---
disp('  ---- 4a. mu_R ----');
for mu_R = mu_R_list
    fprintf('    mu_R = %.1f\n', mu_R);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, base_eta, mu_R, base_cov_R, ...
                                         base_mu_D, base_cov_D, base_N_side, a0, xx_range);
    fname = sprintf('muR%.1f.mat', mu_R);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'mu_R');
end

% --- 4b. CoV of radius cov_R ---
disp('  ---- 4b. cov_R ----');
for cov_R = cov_R_list
    fprintf('    cov_R = %.2f\n', cov_R);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, base_eta, base_mu_R, cov_R, ...
                                         base_mu_D, base_cov_D, base_N_side, a0, xx_range);
    fname = sprintf('covR%.2f.mat', cov_R);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'cov_R');
end

% --- 4c. Mean spacing mu_D ---
disp('  ---- 4c. mu_D ----');
for mu_D = mu_D_list
    fprintf('    mu_D = %.1f\n', mu_D);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, base_eta, base_mu_R, base_cov_R, ...
                                         mu_D, base_cov_D, base_N_side, a0, xx_range);
    fname = sprintf('muD%.1f.mat', mu_D);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'mu_D');
end

% --- 4d. CoV of spacing cov_D ---
disp('  ---- 4d. cov_D ----');
for cov_D = cov_D_list
    fprintf('    cov_D = %.2f\n', cov_D);
    [A1data, aa1PDF, aa2CDF] = func_pdem(base_alph, base_eta, base_mu_R, base_cov_R, ...
                                         base_mu_D, cov_D, base_N_side, a0, xx_range);
    fname = sprintf('covD%.2f.mat', cov_D);
    save(fullfile(save_dir, fname), 'A1data', 'aa1PDF', 'aa2CDF', 'cov_D');
end

fprintf('  Random parameter analysis finished: %.1f sec\n\n', toc(tt));

disp('========== All analyses completed. ==========');