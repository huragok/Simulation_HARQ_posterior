function [map, cost] = solve_QAP_tabu(f, d, pow10_f, pow10_d, n_itr)
%   map = solve_QAP_tabu(f, d, n_itr)
%   A QAP solver function based on tabu search. The QAP problem is 
%   specified in the Koopmans-Beckmann form, i.e. we would like to solve
%   \begin{align}
%       \min_x\sum_{p}\sum_{i}\sum_{q}\sum_{k}a_{pq}b_{ik}x_{pi}x_{qk} \\
%       \mbox{s.t.}\sum_{p}x_{pi} = 1,\sum_{i}x_{pi} = 1, x_{pi}\in\{0, 1\}
%   \end{align}
%   A wrapper function that calls the tabu search method written in c++ by
%   E. Taillard.
% _________________________________________________________________________
%	Inputs:
% 		f:          Q-by-Q matrix, the 2D "Flow" matrix f_{pq}
%       d:          Q-by-Q matrix, the 2D "Distance" matrix d_{ik}
%       pow10_f:    Integer scalar, the power of 10 neede to scale and
%                   round the "Flow" matrix
%       pow10_d:    Integer scalar, the power of 10 neede to scale and
%                   round the "Distance" matrix
%       n_itr:      Integer scalar, number of iterations
%	Outputs:
%		map:    1-by-Q vector, a permutaion of 1 : Q indicating how the Q
%               indices are mapped to constellation points
%       cost:   Scalar, the cost of the resulting mapping scheme
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/08/2015
% Codename: Dunkirk
% _________________________________________________________________________

MAX_LONG = 2^63 - 1; % Maximum number represented by 
f_int = int64(round(f * 10 ^ pow10_f)); % Convert the flow matrix to integer
d_int = int64(round(d * 10 ^ pow10_d)); % Convert the distanve matrix to integer
n_itr = int64(n_itr);%

if (max(max(f_int)) > MAX_LONG) || (max(max(d_int)) > MAX_LONG) % Check whether the scaling is out of range
    error('Integer out of c++ range!');
end

% Call the cmex worker function. All the input must be of type int64
[map, cost] = tabou_QAP_mex(f_int, d_int, n_itr);