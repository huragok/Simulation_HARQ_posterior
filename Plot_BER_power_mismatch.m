clear all;
close all;
clc;

addpath('./functions');

%% 0. Load the data file that contains the test result
% Prior channel
load('Test_201611120047265.mat') % 64QAM, retransmission 2 and 3 0 : 2 : 16, designed with beta = 8 and K = 1
%load('Test_201611118477452.mat')% 64QAM, retransmission 2 - 5 -4 : 2 : 12, designed with beta = 8 and K = 1

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;

nu = test_cases(1).param_origin.nu; % Pathloss factor

theta = test_cases(1).param_origin.theta;
M = test_cases(1).param_origin.M; % Total number of transmissions

n_sigma2 = length(test_cases);
dB_inv_sigma2 = zeros(1, n_sigma2);
for i_sigma2 = 1 : n_sigma2
    dB_inv_sigma2(i_sigma2) = test_cases(i_sigma2).param_origin.dB_inv_sigma2;
end

test_cases_mismatched = test_cases;
d_mismatched = test_cases(1).param_origin.d; % Distance between S and D
K_mismatched = test_cases(1).param_origin.K;

% Posterior channel
load('Test_2016111202530315.mat') % 64QAM, retransmission 2 and 3 0 : 2 : 16, designed with beta = 8 and K = 1
%load('Test_2016111155518076.mat')% 64QAM, retransmission 2 - 5 -4 : 2 : 12, designed with beta = 6.5 and K = 0.8
d = test_cases(1).param_origin.d; % Distance between S and D
K = test_cases(1).param_origin.K;
%% 1. Simulation settings
N_batch = 5; % Number of batches,
N_per_batch = 1e6; % Number of monte-carlo run per batch, restricted by memory size
seed = 8;

%% 2. Initialization
Q = 2 ^ Nbps;

constellation = get_constellation(Nbps, type_mod, 1);

sigma2 = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes

beta = d ^ -nu;

map_noncore = get_map_noncore(Q, M);
map_core = get_map_core(Q, M);

map_QAP = cell(n_sigma2, 1);
for i_sigma2 = 1 : n_sigma2
    map_QAP{i_sigma2} = [1 : Q; test_cases(i_sigma2).map];
end

map_QAP_mismatched = cell(n_sigma2, 1);
for i_sigma2 = 1 : n_sigma2
    map_QAP_mismatched{i_sigma2} = [1 : Q; test_cases_mismatched(i_sigma2).map];
end

%% 3. Not let us test the bit error rate
BER_analytical_noncore = zeros(M, n_sigma2); % The BER upperbound of non-core tested on posterior channel settings
BER_analytical_core = zeros(M, n_sigma2); % The BER upperbound of core tested on posterior channel settings
BER_analytical_QAP = zeros(M, n_sigma2); % The BER upperbound of CoRe designed with posterior channel settings tested on posterior channel settings
BER_analytical_QAP_mismatched = zeros(M, n_sigma2); % The BER upperbound of CoRe designed with prior channel settings tested on posterior channel settings
BER_analytical_QAP_optimistic = zeros(M, n_sigma2); % The BER upperbound of CoRe designed with prior channel settings tested on prior channel settings

BER_MC_noncore = zeros(M, n_sigma2); % The Monte-Carlo BER of non-core tested on posterior channel settings
BER_MC_core = zeros(M, n_sigma2); % The Monte-Carlo BER of core tested on posterior channel settings
BER_MC_QAP = zeros(M, n_sigma2); % The Monte-Carlo BER of CoRe designed with posterior channel settings tested on posterior channel settings
BER_MC_QAP_mismatched = zeros(M, n_sigma2); % The Monte-Carlo BER of CoRe designed with prior channel settings tested on posterior channel settings
BER_MC_QAP_optimistic = zeros(M, n_sigma2); % The Monte-Carlo BER of CoRe designed with prior channel settings tested on prior channel settings

for i_sigma2 = 1 : n_sigma2
    tic

    % Compute the bit error rate using our analytical upper bound
    BER_analytical_noncore(:, i_sigma2) = get_BER_upper_bound(constellation, map_noncore, beta, K, sigma2(i_sigma2));
    BER_analytical_core(:, i_sigma2) = get_BER_upper_bound(constellation, map_core, beta, K, sigma2(i_sigma2));
    BER_analytical_QAP(:, i_sigma2) = get_BER_upper_bound(constellation, map_QAP{i_sigma2}, beta, K, sigma2(i_sigma2));
    BER_analytical_QAP_mismatched(:, i_sigma2) = get_BER_upper_bound(constellation, map_QAP_mismatched{i_sigma2}, beta, K, sigma2(i_sigma2));
    %BER_analytical_QAP_optimistic(:, i_sigma2) = get_BER_upper_bound(constellation, map_QAP_mismatched{i_sigma2}, beta_sr_mismatched, beta_rd_mismatched, g_mismatched(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    
    % Compute the bit error rate using Monte-Carlo simulation
    BER_MC_noncore(:, i_sigma2) = get_BER(constellation, map_noncore, beta, K, theta, sigma2(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC_core(:, i_sigma2) = get_BER(constellation, map_core, beta, K, theta, sigma2(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC_QAP(:, i_sigma2) = get_BER(constellation, map_QAP{i_sigma2}, beta, K, theta, sigma2(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC_QAP_mismatched(:, i_sigma2) = get_BER(constellation, map_QAP_mismatched{i_sigma2}, beta, K, theta, sigma2(i_sigma2), N_per_batch, N_batch, seed);
    %BER_MC_QAP_optimistic(:, i_sigma2) = get_BER(constellation, map_QAP_mismatched{i_sigma2}, beta_sr_mismatched, beta_rd_mismatched, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    
    toc;
    disp(['BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - BER upper bounds, non-CoRe actual: ', num2str(BER_analytical_noncore(:, i_sigma2)')])
    disp([' - BER upper bounds, CoRe actual: ', num2str(BER_analytical_core(:, i_sigma2)')])
    disp([' - BER upper bounds, QAP actual: ', num2str(BER_analytical_QAP(:, i_sigma2)')])
    disp([' - BER upper bounds, QAP mismatched: ', num2str(BER_analytical_QAP_mismatched(:, i_sigma2)')])
    %disp([' - BER upper bounds, QAP optimistic: ', num2str(BER_analytical_QAP_optimistic(:, i_sigma2)')])
    
    disp([' - BER empirical, non-CoRe actual: ', num2str(BER_MC_noncore(:, i_sigma2)')])
    disp([' - BER empirical, CoRe actual: ', num2str(BER_MC_core(:, i_sigma2)')])
    disp([' - BER empirical, QAP actual: ', num2str(BER_MC_QAP(:, i_sigma2)')])
    disp([' - BER empirical, QAP mismatched: ', num2str(BER_MC_QAP_mismatched(:, i_sigma2)')])
    %disp([' - BER empirical, QAP optimistic: ', num2str(BER_MC_QAP_optimistic(:, i_sigma2)')])
end

save(['BER_noise_power_mismatch_', num2str(Q), 'QAM.mat'], 'dB_inv_sigma2', 'BER_analytical_QAP', 'BER_analytical_QAP_mismatched', 'BER_analytical_QAP_optimistic', 'BER_MC_QAP', 'BER_MC_QAP_mismatched', 'BER_MC_QAP_optimistic');

%% Visualization
% The BER upperbound
cmap = [0, 0, 0; 0, 0, 1; 1, 0, 0; 0, 1, 0; 1, 1, 0];
legend_item = cell(4 * M - 3, 1);
h = figure;
semilogy(dB_inv_sigma2, BER_analytical_noncore(1, :), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    semilogy(dB_inv_sigma2, BER_analytical_noncore(m, :), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    semilogy(dB_inv_sigma2, BER_analytical_core(m, :), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical_QAP(m, :), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical_QAP_mismatched(m, :), 's:', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    
    legend_item{4 * m - 6} = ['NM', num2str(m-1)];
    legend_item{4 * m - 5} = ['CR', num2str(m-1)];
    legend_item{4 * m - 4} = ['QAP', num2str(m-1)];
    legend_item{4 * m - 3} = ['QAPM', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
%ylim([1e-6, 10]), xlim([2, 28])
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_noise_power_mismatch_upperbound_', num2str(Q), 'QAM_new.fig']);

%% The empirical BER
h = figure;
semilogy(dB_inv_sigma2, BER_MC_noncore(1, :), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    semilogy(dB_inv_sigma2, BER_MC_noncore(m, :), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    semilogy(dB_inv_sigma2, BER_MC_core(m, :), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_MC_QAP(m, :), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_MC_QAP_mismatched(m, :), 's:', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{4 * m - 6} = ['NM', num2str(m-1)];
    legend_item{4 * m - 5} = ['CR', num2str(m-1)];
    legend_item{4 * m - 4} = ['QAP', num2str(m-1)];
    legend_item{4 * m - 3} = ['QAPM', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
%ylim([1e-6, 10]), xlim([2, 28])
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_noise_power_mismatch_MonteCarlo_', num2str(Q), 'QAM_new.fig']);