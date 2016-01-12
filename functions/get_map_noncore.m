function map = get_map_noncore(Q, M)
%   map = get_map_noncore(Q, M)
%   Get the non-CoRe mapping for a Q-QAM constellation and a total number 
%   of M transmissions, i.e. use Gray mapping for all transmissions
% _________________________________________________________________________
%	Inputs:
% 		Q:      scalar, size of the constellation
%       M:      scalar, total number of transmissions
%	Outputs:
%		map:    M-by-Q matrix, the mapping at each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/14/2015
% Codename: Dunkirk
% _________________________________________________________________________

map = repmat(1 : Q, M, 1);

