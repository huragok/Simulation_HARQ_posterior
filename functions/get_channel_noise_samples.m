function [h_success, v_success, h_failure, v_failure] = get_channel_noise_samples(N, constellation, map, beta, K, theta, sigma_sqr, max_frame, iter_max, coding_rate, nldpc, seed, n_df)
%   [h_success, v_success, h_failure, v_failure] = get_channel_noise_samples(N, constellation, map, beta, sigma_sqr, max_frame, iter_max, coding_rate, nldpc, seed, n_df)
%   Generate the samples of channels and noises corresponding to the
%   successful packet transmissions and the failed packet transmissions.
% _________________________________________________________________________
%	Inputs:
%       N:              scalar, the number of transmissions when we
%                       decide whether a packet has been transmitted 
%                       successfully or not, must be no larger than M and
%                       larget than 0
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta:           Scalar, the variance of the Rayleigh channel
%       K:
%       theta:
%       sigma_sqr:      Scalar, the variance of AWGN noise at the
%                       destination
%       max_frame:      Scalar, number of LDPC frames in simulation.
%       iter_max:       Sclar, maximum iteration time within the iterative 
%                       receiver.
%       coding_rate:    coding rate of LDPC, {1/2,2/3,3/4,5/6}.
%       nldpc:          Scalar, bit length after channel coding, 
%                       mod(nldpc,24)=0.
%       seed:           Scalar, seed for the random number generator
%       n_df:           Scalar, degree of freedom, i.e. the number of
%                       independent fading channel values corresponding to
%                       each transport block. Must be able to divide nldpc
%                       / Q
%	Outputs:
%		h_success:      n_success-by-N matrix, the realization of channel
%                       corresponding to the successful transmissions
%       v_success:      n_success-by-N matrix, the realization of noise
%                       corresponding to the successful transmissions
%       h_failure:      n_failure-by-(N matrix, the realization of channel
%                       corresponding to the failed transmissions
%       v_failure:     n_success-by-N matrix, the realization of noise
%                       corresponding to the failed transmissions
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/04/2015
% Codename: Dunkirk
% _________________________________________________________________________
% The design of this function is based on MIMO_BER_sim function by  
% Weiliang Zeng, 03/04/10.
% _________________________________________________________________________

[M, Q] = size(map);
Nbps = round(log2(Q)); % Number of bit per symbol
rng(seed);

% Config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;

[H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
bits_per_frame = length(H_cols) - length(P_matrix);

%generate bit vectors and mapped symbols used in MAP detector
bit_mat = (dec2bin(0 : Q - 1) > '0') + 0;
bit_mat_anti = 1 - 2 * bit_mat; % antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
sym_mod_mat = constellation(map(1 : N, :)).';
if (N == 1)
    sym_mod_mat = sym_mod_mat(:);
end

% Start the transmission frame by frame
h_success = zeros(0, N);
v_success = zeros(0, N);
h_failure = zeros(0, N);
v_failure = zeros(0, N);

h_mean = sqrt(beta * K / (K + 1)) * exp(theta * 1j);
for i = 1 : max_frame
	if mod(i, 5)==0 % Print a 'x' for every 5 frames
		fprintf('x');
	end
	if mod(i, 400)==0 % Print enter for every 400 frames
		fprintf('\n');
    end

	% source, coding, random interlever and modulation
    data = round(rand(1, bits_per_frame)); % Randomly generated information-bearing bit, 1-by-bits_per_frame
    codewordTemp = LdpcEncode(data, H_rows, P_matrix ); % LDPC encoded codeword, 1-by-nldpc
    codeword_invTemp = randintrlv(codewordTemp, 0);  % Interleaved codeword, 1-by-nldpc
    coded_index = bit2idx(codeword_invTemp, Nbps);

    transmit_Mod = zeros(N, ceil(nldpc / Nbps)); % Modulated symbol for each (re)transmission
    for m = 1 : N % 
        transmit_Mod(m, :) = constellation(map(m, coded_index)); % Modulation and normalization
    end
	
	% Lets start simulating the transmission
	numSymbol = size(transmit_Mod, 2); % Number of transmitted symbols per stream per frame
	y = zeros(N, numSymbol); % The received signal at each (re)transmission
    % Generate the channels. Assume channel to be independently fading
    % across symbols
    h_unique = h_mean + sqrt(beta / 2 / (K + 1)) * (randn(N, n_df) + 1i * randn(N, n_df)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
    v = sqrt(sigma_sqr / 2) * (randn(N, numSymbol) + 1i * randn(N, numSymbol)); % Noise samples at the relay
    
    n_rep = numSymbol / n_df;
    h = repmat(h_unique, 1, n_rep);
    h = complex(real(h), imag(h));
    for m = 1 : N 	
		% Received signal and equivalent channel generation
        y(m, :) = h(m, :) .* transmit_Mod(m, :) + v(m, :);
        y(m, :) = complex(real(y(m, :)), imag(y(m, :))); % make sure y is complex, used for MAP
    end
    
     %iterative receiver
    LextC = zeros(1, nldpc);

    error_perFrame_perMS = zeros(iter_max, 1);
    for iter = 1 : iter_max

        LextDemodulation = MAP_demod(y, h, bit_mat_anti, LextC, sym_mod_mat, sigma_sqr);

        LextDemo_deinv = randdeintrlv(LextDemodulation, 0); %de-interleave
        [LLR_output_tmp, errors] = MpDecode(LextDemo_deinv, H_rows, H_cols, max_iterations, decoder_type, 1, 1, data); %ldpc decoder
        error_perFrame_perMS(iter) = errors(end); %count errors in each iteration
        LLR_output = LLR_output_tmp(end, :);
        LextDecoder = LLR_output - LextDemo_deinv; %calculate extrinic information from channel decoder
        LextC = randintrlv(LextDecoder, 0); %interleave
    end
    error_all = error_perFrame_perMS(end);

    % Detection of the current frame is successful, compute the
    % effective thourhput
    if error_all == 0
        h_success = [h_success; h_unique.'];

        
        % Randomly samples the noise samples (must correspond to the value)
        idx = randi(n_rep, N, n_df) - 1;
        v_success_frame = zeros(n_df, N);
        for m = 1 : N
            for i_df = 1 : n_df
                v_success_frame(i_df, m) = v(m, i_df + n_df * idx(m, i_df));
            end
        end

        v_success = [v_success; v_success_frame];
    else
        h_failure = [h_failure; h_unique.'];

        idx = randi(n_rep, N, n_df) - 1;
        v_failure_frame = zeros(n_df, N);
        for m = 1 : N
            for i_df = 1 : n_df
                v_failure_frame(i_df, m) = v(m, i_df + n_df * idx(m, i_df));
            end
        end

        v_failure = [v_failure; v_failure_frame];
    end
    
end