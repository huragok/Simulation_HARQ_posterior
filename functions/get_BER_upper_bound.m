function BER = get_BER_upper_bound(constellation, map, beta, K, sigma2)
%   BER = get_BER_upper_bound(constellation, map, beta, K, theta, sigma2)
%   Get the BER upper bounds after each retransmission based on our
%   approximation
% _________________________________________________________________________
%	Inputs:
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta:           Scalar, the variance of the Rayleigh channel from
%                       source to relay
%       K:
%       sigma2:         Scalar, the variance of AWGN noise at the
%                       destination
%	Outputs:
%		BER:            M-by-1 vector, the BER after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/11/2016
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map); % Number of transmissions and size of the constellation
Nbps = round(log2(Q)); % Number of bit per symbol

xpcd_PBER = get_hamming_dist(Nbps) / 2 / Q / Nbps; % Initialization
BER = zeros(M, 1);
for m = 1 : M
    dist_sqr = abs(repmat(constellation(map(m, :)), 1, Q) - repmat(constellation(map(m, :)).', Q, 1)) .^ 2; % A Q-by-Q matrix containing the distance square mesurements
    E = get_factor_PEP_update(dist_sqr, beta, K, sigma2); % Get this thing fully vectorized
    xpcd_PBER = xpcd_PBER .* E; % Update the expected number of pairwise error bit
    BER(m) = sum(sum(xpcd_PBER));
end

