function hamming_dist = get_hamming_dist(Nbps)
%   hamming_dist = get_hamming_dist(Nbps)
%   Get the hamming distance between integers in [0, 2 ^ Nbps - 1]
%   Reused from function get_n_diff_bits() from Q3AP simulation
% _________________________________________________________________________
%	Inputs:
% 		Nbps:           scalar, number of bits per symbol
%	Outputs:
%		hamming_dist:   Q-by-Q vector, the row major order vector
%                       of the hamming distance matrix
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = 2 ^ Nbps;
hamming_dist = zeros(Q, Q);
for i = 0 : Q - 1
    for j = 0 : Q - 1
        hamming_dist(i + 1, j + 1) = sum(de2bi(bitxor(i, j)));
    end
end
