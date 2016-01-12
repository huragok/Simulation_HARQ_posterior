function BER = get_BER(constellation, map, beta, K, theta, sigma2, N_per_batch, N_batch, seed)
%   BER = get_BER(constellation, map, beta, K, theta, sigma2, N_per_batch, N_batch, seed)
%   Get the BER by actually running Monte-Carlo simulation on the ML
%   demodulators
% _________________________________________________________________________
%	Inputs:
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta:
%       K:
%       theta:
%       sigma2:
%       N_per_batch:    Scalar, number of Monte-Carlo run per batch (size 
%                       of vectorization)
%       N_batch:        Scalar, number of batches (for-loop size)
%       seed:           Scalar, seed for the random number generator
%	Outputs:
%		BER:			M-by-1 vector, the BER after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/11/2015
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map);
Nbps = round(log2(Q)); % Number of bit per symbol
symbols_mapped = constellation(map); % The mapped symbols at all transmissions

BER = zeros(M, N_batch);
B = get_hamming_dist(Nbps); % The hamming distance matrix

h_mean = sqrt(beta * K / (K + 1)) * exp(theta * 1j);

rng(seed);
for i_batch = 1 : N_batch
    p = randi(Q, 1, N_per_batch); % Generate the random transmitted index
    d = zeros(Q, N_per_batch); % The distance measurement used for the ML demodulator

    for m = 1 : M
        h =  h_mean + sqrt(beta / 2 / (K + 1)) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));
        symbol = symbols_mapped(m, p); % Generate the random symbol

        y = h .* symbol + sqrt(sigma2 / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));

        d = d + abs(repmat(y, Q, 1) - symbols_mapped(m, :).' * h) .^ 2 / sigma2; % Compute the ML measurement: the weighted distance square between the received signal and the Q symbols for each of N realization
        [~, p_demod] = min(d, [], 1);
        BER(m, i_batch) = mean(B((p_demod - 1) * Q + p)) / Nbps; % BER for the m-th transmission
    end
end

BER = mean(BER, 2);