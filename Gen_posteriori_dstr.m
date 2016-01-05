clear all;
close all;
clc;

addpath('functions');

%% 1. Simulation settings
% 64QAM
load('./Test_201614185534739.mat');
% Non-block fading
% n_df = 2400 / 6;
% dB_inv_sigma2 = -4.7; % 9.4, 1.6, -2.4, -4.7
% N = 4; % 1, 2, 3, 4
% max_frame = 50;

% pure block fading
n_df = 1;
dB_inv_sigma2 = -5.2; % 8.3, 0.6, -3.2, -5.2
N = 4; % 1, 2, 3, 4
max_frame = 20000;

% 100 RB
% n_df = 100;
% dB_inv_sigma2 = -4.7; % 9.4, 1.6, -2.5, -4.7
% N = 4; % 1, 2, 3, 4
% max_frame = 200;

% 10 RB
% n_df = 10;
% dB_inv_sigma2 = -4.7; % 9.2, 1.5, -2.5, -4.7
% N = 4; % 1, 2, 3, 4
% max_frame = 2000;

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
d = test_cases(1).param_origin.d; % Distance between S and R, R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions


iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, 1);

beta = d ^ -nu;

map  = [1 : Q; test_cases(1).map];

sigma2 = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes

%% 3. Generate the channel samples corresponding to the successful transmission and the failed transmissions
[h_success, v_success, h_failure, v_failure] = get_channel_noise_samples(N, constellation, map, beta, sigma2, max_frame, iter_max, coding_rate, nldpc, seed, n_df);

%% 4. Visualization, plot the mean and covaraince matrix
points_success_channel = [real(h_success), imag(h_success)];
points_failure_channel = [real(h_failure), imag(h_failure)];

points_success_noise = [real(v_success), imag(v_success)];
points_failure_noise = [real(v_failure), imag(v_failure)];

% The covariance plot of the channel samples
colormap('hot');

cov_success_channel = cov(points_success_channel);
cov_failure_channel = cov(points_failure_channel);

% figure;
% imagesc(cov_success_channel);
% axis equal;
% xlabel('[real(h(1 : N)), imag(h(1 : N))]');
% ylabel('[real(h(1 : N)), imag(h(1 : N))]');
% xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
% set(gca, 'XTick', 1 : 2 * N);
% set(gca, 'YTick', 1 : 2 * N);
% set(gca, 'Fontsize', 18);
% colorbar;
% caxis([0, max(max(max(cov_success_channel)), max(max(cov_failure_channel)))]);
% 
% figure;
% imagesc(cov_failure_channel);
% axis equal;
% xlabel('[real(h(1 : N)), imag(h(1 : N))]');
% ylabel('[real(h(1 : N)), imag(h(1 : N))]');
% xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
% set(gca, 'XTick', 1 : 2 * N);
% set(gca, 'YTick', 1 : 2 * N);
% set(gca, 'Fontsize', 18);
% colorbar;
% caxis([0, max(max(max(cov_failure_channel)), max(max(cov_failure_channel)))]);
% 
% % The covariance plot of the noise samples
% cov_success_noise = cov(points_success_noise);
% cov_failure_noise = cov(points_failure_noise);
% 
% figure;
% imagesc(cov_success_noise);
% axis equal;
% xlabel('[real(n(1 : N)), imag(n(1 : N))]');
% ylabel('[real(n(1 : N)), imag(n(1 : N))]');
% xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
% set(gca, 'XTick', 1 : 2 * N);
% set(gca, 'YTick', 1 : 2 * N);
% set(gca, 'Fontsize', 18);
% colorbar;
% caxis([0, max(max(max(cov_success_noise)), max(max(cov_failure_noise)))]);
% 
% figure;
% imagesc(cov_failure_noise);
% axis equal;
% xlabel('[real(n_R(1 : N)), imag(n_R(1 : N))]');
% ylabel('[real(n_R(1 : N)), imag(n_R(1 : N))]');
% xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
% set(gca, 'XTick', 1 : 2 * N);
% set(gca, 'YTick', 1 : 2 * N);
% set(gca, 'Fontsize', 18);
% colorbar;
% caxis([0, max(max(max(cov_failure_noise)), max(max(cov_failure_noise)))]);

%% 5. Output the result to a xml file used for normality test, each row corresponds to 1 observation in the order of real(h(1)), imag(h(1)), ...real(h(N)), imag(h(N)), real(g(1)), imag(g(1)), ...real(g(N)), imag(g(N)) separated by comma
doc = com.mathworks.xml.XMLUtils.createDocument('ChannelSamples');

channel_samples = doc.getDocumentElement;
channel_samples.setAttribute('nSuccess', num2str(size(h_success, 1)));
channel_samples.setAttribute('nFailure', num2str(size(h_failure, 1)));
channel_samples.setAttribute('beta', num2str(beta));
channel_samples.setAttribute('N', num2str(N));
channel_samples.setAttribute('Nprb', num2str(n_df));
channel_samples.setAttribute('sigma2', num2str(sigma2));

% The success part
success = doc.createElement('Success');
for i_success = 1 : size(h_success, 1)
    entry = doc.createElement('Entry');
    entry.appendChild(doc.createTextNode(num2str([real(h_success(i_success, :)),...
                                                  imag(h_success(i_success, :)),...
                                                  real(v_success(i_success, :)),...
                                                  imag(v_success(i_success, :))], '%-f ')));
    success.appendChild(entry);
end
channel_samples.appendChild(success);

% The failure part
failure = doc.createElement('Failure');
for i_failure = 1 : size(h_failure, 1)
    entry = doc.createElement('Entry');
    entry.appendChild(doc.createTextNode(num2str([real(h_failure(i_failure, :)),...
                                                  imag(h_failure(i_failure, :)),...
                                                  real(v_failure(i_failure, :)),...
                                                  imag(v_failure(i_failure, :))], '%-f ')));
    failure.appendChild(entry);
end
channel_samples.appendChild(failure);

xmlwrite(['samples_', num2str(N), '_', num2str(n_df), '.xml'], doc);
