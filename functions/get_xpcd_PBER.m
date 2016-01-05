function xpcd_PBER_new = get_xpcd_PBER(xpcd_PBER, E, map)
%   xpcd_PBER = get_xpcd_PBER(xpcd_PBER, E, map)
%   Update the expected pairwise bit error rate given the update another
%   round of retransmission
% _________________________________________________________________________
%	Inputs:
%       xpcd_PBER:          Q-by-Q matrix, the original expected pairwise
%                           bit error rate in row major order
%       E:                  Q-by-Q matrix, the updating matrix in row
%                           major order.
%       map:                Q-by-Q vector, how does the Q indices are 
%                           mapped to constellation points
%	Outputs:
%       xpcd_PBER_new:      Q-by-Q matrix, the updated expected pairwise
%                           bit error rate in row major order
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = length(map);
xpcd_PBER_new = zeros(Q, Q);

for p = 1 : Q
    for q = 1 : Q
        if p ~= q
            xpcd_PBER_new(p, q) = xpcd_PBER(p, q) * E(map(p), map(q));
        end
    end
end