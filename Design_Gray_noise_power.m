clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% Constellation specification
Nbps = 6; % 4, 5, 6
type_mod = 'QAM';

% Now we are considering a one hop link (S-D)

dB_inv_sigma2 = [0]; % 1/sigma2 in dB
d = 0.5; % Distance between S and R, R and D

K = 1; % The Rician coefficients
theta = 0; % The phase of the LOS component

nu = 3; % Pathloss factor
M = 5; % Number of retransmission

epsilon = 0.01; % Tolerance to control the error of scaling the 2 cost matrices to integer
n_itr = 1000000; % Number of iterations for the tabu QAP solver

%% 2. Initialization: generate and save all test cases
test_cases = construct_test_cases(Nbps, type_mod, dB_inv_sigma2, d, nu, K, theta, M, true);
n_case = length(test_cases);
time_step = regexprep(num2str(clock),'[^\w'']',''); % The time step used to label all saved files as a suffix
Q = 2 ^ Nbps;

%% 3. Start the computation
% Compute the expected number of pairwise error bit based on the cost 
% matrix for the previous transmissions, or initialize this value to be 1 /
% 2 of the hamming distance matrix

for i_case = 1 : n_case
    disp(['Test case ', num2str(i_case), '/', num2str(n_case)]);
    
    E = test_cases(i_case).param_derived.E;
    xpcd_PBER = get_hamming_dist(test_cases(i_case).param_origin.Nbps) / 2 / test_cases(i_case).param_derived.Q / test_cases(i_case).param_origin.Nbps; % Initialize the expected pairwise BER before any transmission
    xpcd_PBER = xpcd_PBER .* E; % The expected pairwise BER after the first transmission (Gray mapping)
    
    test_cases(i_case).map = zeros(test_cases(i_case).param_origin.M - 1, test_cases(i_case).param_derived.Q);
    for m = 2 : test_cases(i_case).param_origin.M
  
        map = 1 : Q;
        test_cases(i_case).map(m - 1, :) = map;
        
        % Update the expected PBER
        xpcd_PBER = get_xpcd_PBER(xpcd_PBER, E, map);

        % toc;
        disp([' - Transmission ', num2str(m), ': BER = ', num2str(sum(sum(xpcd_PBER)))]);
    end
    
    % Save the test results regularly in case of a crash
    save(['Test_', time_step, '.mat'], 'test_cases');
    disp(['Test case ', num2str(i_case), '/', num2str(n_case), ' saved.']);
end
