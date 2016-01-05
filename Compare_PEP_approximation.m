clear all;
close all;
clc;

addpath('./functions');


%% 1. Simulation settings
% Constellation specification
Nbps = 6; % 4, 5, 6
type_mod = 'QAM';


dB_inv_sigma2 = [-2 : 1.5 : 10]; % 1/sigma2 in dB, 64 QAM
d = 0.5; % Distance between S and R, R and D

nu = 3; % Pathloss factor
K = 1;
theta = 0;
% We set M = 2, since for lar Number of retransmission

N_batch = 5; % Number of batches,
N_per_batch = 1e6; % Number of monte-carlo run per batch, restricted by memory size
seed = 8;

% 64QAM
p0 = 0; q0 = 6;
p1 = 0; q1 = 48;

% 16QAM
% p0 = 0; q0 = 3;
% p1 = 0; q1 = 15;

%% 2. Initialization: generate and save all test cases
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, 1);

sigma2 = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes are the same n_sigma2-by-1 vector
n_sigma2 = length(dB_inv_sigma2);


beta = d.^ -nu; % n_d-by-1 vector
h_mean = sqrt(beta * K / (K + 1)) * exp(theta * 1j);

dist_sqr0 = abs(constellation(p0 + 1) - constellation(q0 + 1)) ^ 2;
dist_sqr1 = abs(constellation(p1 + 1) - constellation(q1 + 1)) ^ 2;

pep_num = zeros(n_sigma2, 1); % PEP evaluated according to its definition in (7) numerically
pep_chernoff_approx = zeros(n_sigma2, 1);  % PEP evaluated using Chernoff bound and the approximation by Proposition 1
pep_doubleexp_approx = zeros(n_sigma2, 1); % PEP evaluated using more accurate sum of two exponential approximation and the approximation by Proposition 1
pep_chernoff_num = zeros(n_sigma2, 1); % PEP evaluated using Chernoff bound with E evaluated numerically
pep_doubleexp_num = zeros(n_sigma2, 1); % PEP evaluated using the more accurate sum of two exponential approximation with E evaluated numerically
rng(seed);

for i_batch = 1 : N_batch
    % generate the random channels for the first transmission
    delta1_0 = abs(h_mean + sqrt(beta / 2 / (K + 1)) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    
    % generate the random channels for the first retransmission
    delta1_1 = abs(h_mean + sqrt(beta / 2 / (K + 1)) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
   
    for i_sigma2 = 1 : n_sigma2
        tmp0 = delta1_0 ./ sigma2(i_sigma2);
        tmp1 = delta1_1 ./ sigma2(i_sigma2);
        
        pep_num(i_sigma2) = pep_num(i_sigma2) + mean(qfunc(sqrt(dist_sqr0 * tmp0 / 2 + dist_sqr1 * tmp1 / 2)));
        pep_chernoff_num(i_sigma2) = pep_chernoff_num(i_sigma2) + 0.5 * mean(exp(-dist_sqr0 * tmp0 / 4)) * mean(exp(-dist_sqr1 * tmp1 / 4));
        pep_doubleexp_num(i_sigma2) = pep_doubleexp_num(i_sigma2) + mean(exp(-dist_sqr0 * tmp0 / 4)) * mean(exp(-dist_sqr1 * tmp1 / 4)) / 12 + mean(exp(-dist_sqr0 * tmp0 / 3)) * mean(exp(-dist_sqr1 * tmp1 / 3)) / 4;
        
        pep_chernoff_approx(i_sigma2) =...
        pep_chernoff_approx(i_sigma2) +...
        0.5 * get_factor_PEP_update_param(dist_sqr0, beta, K, sigma2(i_sigma2), 4)...
            * get_factor_PEP_update_param(dist_sqr1, beta, K, sigma2(i_sigma2), 4);
        
        pep_doubleexp_approx(i_sigma2) =...
        pep_doubleexp_approx(i_sigma2) +...
        1 / 12 * get_factor_PEP_update_param(dist_sqr0, beta, K, sigma2(i_sigma2), 4)...
               * get_factor_PEP_update_param(dist_sqr1, beta, K, sigma2(i_sigma2), 4)...
        +1 / 4 * get_factor_PEP_update_param(dist_sqr0, beta, K, sigma2(i_sigma2), 3)...
               * get_factor_PEP_update_param(dist_sqr1, beta, K, sigma2(i_sigma2), 3);
    end
end

pep_num = pep_num / N_batch;
pep_chernoff_approx = pep_chernoff_approx / N_batch;
pep_doubleexp_approx = pep_doubleexp_approx / N_batch;
pep_chernoff_num = pep_chernoff_num / N_batch;
pep_doubleexp_num = pep_doubleexp_num / N_batch;

figure;
semilogy(dB_inv_sigma2, pep_num, 'k+:', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2, pep_doubleexp_num, 'rs-.', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_chernoff_num, 'ro-.', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_doubleexp_approx, 'bs-', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_chernoff_approx, 'bo-', 'linewidth', 2);
legend('Qfunc+Num', 'ExpSum+Num', 'Chernoff+Num', 'ExpSum+Approx', 'Chernoff+Approx');
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('PEP');