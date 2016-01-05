function idxs = bit2idx(bits, Nbps)
%   idxs = bit2idx(bits)
%   Convert a sequence of bits to a sequence of indices of 1,...,Q. If
%   mod(length(bits), Nbps) ~= 0, then bits is is zero padded at the end
%   and then converted to idxs
% _________________________________________________________________________
%	Inputs:
%       bits:       1-by-nldpc vector, the encoded bits
%       Nbps:       Scalar, number of bits per symbol
%	Outputs:
%		idxs:		1-by-n_idx vector, the index symbols corresponding to
%                   the zero padded coded bits
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/25/2015
% Codename: Dunkirk
% _________________________________________________________________________

nldpc = length(bits);
n_idx = floor(nldpc / Nbps);
bits = [bits, zeros(1, n_idx * Nbps - nldpc)];
idxs = 2 .^ (Nbps - 1 : -1 : 0) * reshape(bits, Nbps, n_idx) + 1;

