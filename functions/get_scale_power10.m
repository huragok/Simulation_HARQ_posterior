function pwr_a = get_scale_power10(a, epsilon)
%   pwr_a = get_scale_power10(a, epsilon)
%   Estimate the power of 10 needed to scale both input matrix to integer
% _________________________________________________________________________
%	Inputs:
% 		a:              Q-by-Q postitive vector, the row-major order 
%                       matrix to be scaled
%       epsilon:        scalar, the relative tolerance of the rounding
%                       error
%	Outputs:
%		pwr_a:          Integer scalar, the minimun power of 10 needed to 
%                       achieve the tolerance
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/08/2015
% Codename: Dunkirk
% _________________________________________________________________________

min_a = min(a(a > 0)); % Minimum non-zero elements in pwr_a
pwr_a = ceil(log10(1 / epsilon / min_a));
