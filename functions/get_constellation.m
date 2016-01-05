function X = get_constellation(Nbps, type_mod, pwr)
%   X = get_constellation(Nbps, type_mod)
%   Generate all the constellation points using Matlab comm toolbox, using
%   the default Gray mapping.
% _____________________________________________________________________________
%	Inputs:
% 		Nbps:       scalar, number of bits per symbol
%       type_mod:   string, 'PSK' for phase shift keying and QAM for quadrature 
%                   amplitude modulation
%       pwr:        scalar, average power of the modulated symbols
%	Outputs:
%		X:			2 ^ Nbps-by-1 vector, the modulated constellations
% _____________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/18/2014
% Codename: Dunkirk
% _____________________________________________________________________________

Q = 2 ^ Nbps; % Number of constellation points

if strcmpi(type_mod, 'PSK')
    hmod = comm.PSKModulator(Q);
elseif strcmpi(type_mod, 'QAM')
    if Q == 64
        CustomSymbolMapping = [47, 46, 42, 43, 59, 58, 62, 63,...
                               45, 44, 40, 41, 57, 56, 60, 61,...
                               37, 36, 32, 33, 49, 48, 52, 53,...
                               39, 38, 34, 35, 51, 50, 54, 55,...
                               7, 6, 2, 3, 19, 18, 22, 23,...
                               5, 4, 0, 1, 17, 16, 20, 21,...
                               13, 12, 8, 9, 25, 24, 28, 29,...
                               15, 14, 10, 11, 27, 26, 30, 31];
    elseif Q == 16
        CustomSymbolMapping = [11, 10, 14, 15, 9, 8, 12, 13, 1, 0, 4, 5, 3, 2, 6, 7];
    end
    hmod = comm.RectangularQAMModulator(Q, 'SymbolMapping', 'Custom', 'CustomSymbolMapping', CustomSymbolMapping);
else
    error('type_mode must either be PSK or QAM')
end

X = step(hmod, (0 : Q - 1)');
pwr_tmp = X' * X / Q; % Average power of the constellation points
X = X * sqrt(pwr / pwr_tmp); % Normalization

end

