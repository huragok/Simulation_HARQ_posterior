function E = get_factor_PEP_update(dist_sqr, beta, sigma_sqr)
%   E = get_factor_PEP_update(dist_sqr, beta, sigma_sqr)
%   Get the factor to be multiplied to PEP between indices p and q for the
%   first (M-1) retransmissions to compute the PEP for the first M 
%   retransmissions
% _________________________________________________________________________
%	Inputs:
% 		dist_sqr:       Vector, the square norm of the 2 comstellation 
%                       points to which p and q are mapped in the M-th
%                       retransmission
%       beta:           Scalar, the variance of the Rayleigh channel
%       sigma_sqr:      Scalar, the variance of AWGN noise at the
%                       destination
%	Outputs:
%		E:              Scalar, the factor used to update PEP for the first
%                       (M-1) retransmissions into PEP for the first M
%                       retransmissions between p and q
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/04/2016
% Codename: Dunkirk
% _________________________________________________________________________

E = 1 ./ (1 + dist_sqr * beta / 4 / sigma_sqr); % Note the difference between the exponential integral function definition in the reference and Matlab

