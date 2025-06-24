clear;
clc;

% Simulation Parameters
SNRs        = 0:1:25;

constellation       = struct();
constellation.BPSK  = [-1, 1];
constellation.QPSK  = [-1+1i, 1+1i, -1-1i, 1-1i]./sqrt(2);
constellation.QAM16 = [-3+3i, -3+1i, -3-3i, -3-1i, -1+3i, -1+1i, -1-3i, -1-1i, 3+3i, 3+1i, 3-3i, 3-1i, 1+3i, 1+1i, 1-3i, 1-1i]./sqrt(10);
names       = fieldnames(constellation);

N_sc        = 20;
num_bits    = 4;
num_frames  = 3;
L           = 5;
syms        = constellation.(names{log2(num_bits) + 1});

% Baseline - Known Channel, Single UE, No diversity
ber_base    = baseline(N_sc, L, num_frames, 1, constellation.BPSK, SNRs);
ber_qp      = baseline(N_sc, L, num_frames, 2, constellation.QPSK, SNRs);
ber_q16     = baseline(N_sc, L, num_frames, num_bits, syms, SNRs);

% Final Design - 16QAM, Unknown Channel, Multiple UE, Repetition Coding
ber_est_multi_diversity = main_ofdm_sys_est_rep_multi_user(N_sc, L, num_frames, num_bits, syms, SNRs);


% Plot baseline system with different modulations and final system
figure;
semilogy(SNRs, ber_base);
xlabel("SNR");
ylabel("BER");
title(sprintf("BER vs SNR"));
grid on;
hold on;
semilogy(SNRs, ber_qp);
semilogy(SNRs, ber_q16);
semilogy(SNRs, ber_est_multi_diversity);
legend("Baseline", "Baseline with QPSK", ...
    "Baseline with 16QAM", "16QAM, Unknown Channel, 2UE, Diversity");
hold off;


%% Number of Pilots vs BER
% Study trade-off between number of pilots and BER
% Channel Estimation, 16QAM modulation,
% 1 UE, 10dB SNR, Pilots from 1 to N_sc

ber_trade_off   = trade_off_system(N_sc, L, num_frames, num_bits, syms);

figure;
semilogy(1:N_sc, ber_trade_off);
title("Number of Pilots vs BER")
xlabel("Pilots");
ylabel("BER");
grid on;

%% Functions
function y = convolution(h, X)
    y = zeros(size(X));    
    for n = 1:length(y)
        for p = 1:length(h)
            if (n - p + 1 > 0) && (n - p + 1 <= length(X))
                y(n) = y(n) + X(n - p + 1) * h(p);
            end
        end
    end  
end

function [y_noisy, h] = simulate_channel(x, N_sc, L, SNR_dB)
% SIMULATE_CHANNEL
% Description: Simulates a multipath fading channel with additive white 
% Gaussian noise (AWGN). Computes the channel impulse response, generates 
% noisy channel output, and calculates noise power.
% Input:
%   - x: Input signal vector.
%   - N_sc: Number of subcarriers in each frame.
%   - L: Length of the channel impulse response (number of taps).
%   - SNR_dB: Signal-to-Noise Ratio in decibels.
% Output:
%   - y_noisy: Output signal after passing through channel and adding noise.
%   - h: Channel impulse response.    
    
    % Simulate channel with random taps
    h = (randn(L, 1) + 1i .* randn(L, 1)) .* sqrt(1/(2*L));        
    y = convolution(h, x);

    % Get Expected Symbol Energy (approx. 1)
    Es = 0;
    for frame = 1:3
        start_idx = (N_sc + L - 1)*(frame - 1) + (L - 1);
        end_idx = frame *(N_sc + L - 1);
        Es = Es + sum(abs(x(start_idx:end_idx)) .^ 2) / N_sc; % Symbol energy
    end
    Es = Es / 3;

    % Noise power
    SNR_linear = 10^(SNR_dB/10);
    N0 = Es / SNR_linear;

    % Generate Gaussian noise
    noise = sqrt(N0/2) .* (randn(size(y)) + 1i .* randn(size(y)));
    
    % Add noise to simulated channel output
    y_noisy = y + noise;
end

function decoded_bits = decode_symbols(symbols, num_bits)
% Description: Decodes complex symbols into binary bits for different 
% modulation schemes based on the decoding regions for each scheme
% Input:
%   - symbols: Received complex symbol array.
%   - num_bits: Number of bits per symbol (1, 2, or 4).
% Output:
%   - decoded_bits: Binary bit array corresponding to input symbols.
    
    switch num_bits
        case 1
            decoded_bits = real(symbols) > 0;
        case 2
            re_bits = real(symbols) > 0;
            im_bits = imag(symbols) < 0;
            decoded_bits = [im_bits, re_bits];
        case 4
            % [bit1 bit2 bit3 bit4]
            bit1 = real(symbols) > 0; % 0xxx or 1xxx
            bit2 = abs(real(symbols)) < 2/sqrt(10); % x0xx or x1xx
            bit3 = imag(symbols) < 0; % xx0x or xx1x
            bit4 = abs(imag(symbols)) < 2/sqrt(10); % xxx0 or xxx1

            decoded_bits = [bit1, bit2, bit3, bit4];
        otherwise
            disp("Only 1, 2, or 4 bits accepted");
    end

end

function [X_tx, bitstring] = gen_baseline(N_sc, cp_length, num_frames, syms, IDFT_matrix)
% GEN_BASELINE
% Generates an OFDM signal for a single user.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   cp_length   - Length of Cyclic Prefix
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx        - Serialized OFDM signal with all frames
%   bitstring   - Random bitstring transmitted
%

    % Initialize parameters and output signal
    num_bits = log2(length(syms)); % Number of bits per symbol
    X_tx = zeros((N_sc + cp_length) * num_frames, 1); % OFDM signal
    
    % Generate random bit sequence
    total_bits = N_sc * num_frames * num_bits;
    bitstring = randi([0, 1], total_bits, 1);

    % Process each frame
    for frame_idx = 1:num_frames
        % Extract and reshape bits for the current frame
        start_idx = (frame_idx - 1) * (N_sc * num_bits) + 1;
        end_idx = frame_idx * (N_sc * num_bits);
        frame_bits = bitstring(start_idx:end_idx);
        frame_bits = reshape(frame_bits, [N_sc, num_bits]);

        % Map bits to modulation symbols
        X_freq = syms(bit2int(frame_bits.', num_bits) + 1).';

        % Perform IDFT to convert to time domain
        X_time = IDFT_matrix * X_freq;

        % Add cyclic prefix
        X_cp = [X_time(end - cp_length + 1:end); X_time];

        % Append the cyclic-prefixed signal to the output
        start_output_idx = (frame_idx - 1) * (N_sc + cp_length) + 1;
        end_output_idx = frame_idx * (N_sc + cp_length);
        X_tx(start_output_idx:end_output_idx) = X_cp;
    end
end

function [BER, Total_capacity] = baseline(N_sc, L, num_frames, num_bits, syms, SNRs)
% BASELINE
% This function simulates a basic ofdm system with no diversity, known 
% channel, and a single user. The simulation runs over a range of SNRs with
% BER returned at the end of the simulation.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   num_bits    - Number of bits in one symbol
%   syms        - Modulation symbol mapping
%   SNRs        - SNR range to test the simulation
%
% OUTPUT:
%   BER         - Bit Error Rate over all SNR range
%
    cp_length = L - 1;
    BER = zeros(size(SNRs));
    Total_capacity = zeros(size(SNRs));
    
    % Generate IDFT and DFT matrices
    IDFT_matrix = exp(1i * 2 * pi * (0:N_sc-1)' * (0:N_sc-1) / N_sc) / sqrt(N_sc);
    DFT_matrix = conj(IDFT_matrix);

    num_iters = 1e5;
    
    for k = 1:length(SNRs)
        SNR_dB = SNRs(k);
        num_errors = 0;
       
        for m = 1:num_iters
            % Generate OFDM signal
            [x, bitstring] = gen_baseline(N_sc, cp_length, num_frames, syms, IDFT_matrix);
            
            % Simulate channel and add noise
            [y_noisy, h] = simulate_channel(x, N_sc, L, SNR_dB);
    
            % Channel Frequency Response
            H_freq = sqrt(N_sc) .* DFT_matrix(1:N_sc, 1:L) * h;
            
            % Preallocate reconstructed bitstring
            hatbitstring = zeros(N_sc*num_frames*num_bits, 1);
            
            % Retrieve data and compute BER
            for j = 1:num_frames
                % Remove CP
                y_no_cp = y_noisy((N_sc + cp_length)*(j-1)+cp_length+1:(N_sc + cp_length)*(j));
        
                % DFT at receiver
                Y_freq = DFT_matrix * y_no_cp;   
                
                % Channel equalization
                R = Y_freq ./ H_freq;
               
                % Decode symbols
                decoded_bits = decode_symbols(R, num_bits);
                
                % Reconstruct bitstring
                hatbitstring(N_sc*num_bits*(j-1) + 1: N_sc*num_bits*j) = reshape(decoded_bits, [N_sc*num_bits, 1]);
            end
            
            % Error Accumulation
            num_errors = num_errors + sum(hatbitstring ~= bitstring, "all");
        end

        % BER calculation
        BER(k) = num_errors / N_sc / num_frames / num_bits / num_iters;
        
    end
end

function indices = get_user_indices(N_sc, num_frames, pilot, repition)
    indices = cell(num_frames, 2);
    idx1 = 1:N_sc/2;
    idx2 = N_sc/2 + 1:N_sc;
    if ~pilot && ~repition
        indices{1, 1} = idx1;
        indices{1, 2} = idx2;
        indices{2, 1} = idx2;
        indices{2, 2} = idx1;
        indices{3, 1} = idx1;
        indices{3, 2} = idx2;
    elseif ~repition
        indices{1, 1} = [2, 3, 4, 6, 7, 8, 10];
        indices{1, 2} = [11, 12, 14, 15, 16, 18, 19];
        indices{2, 1} = idx2;
        indices{2, 2} = idx1;
        indices{3, 1} = idx1;
        indices{3, 2} = idx2;
    else
        indices{1, 1} = [2, 3, 4, 10, 11, 12, 18, 19, 20];
        indices{1, 2} = [6, 7, 8, 14, 15, 16];

        indices{2, 1} = [1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20];
        indices{2, 2} = setdiff(1:20, indices{2, 1});

        indices{3, 1} = [5, 6, 7, 8, 13, 14, 15, 16];
        indices{3, 2} = setdiff(1:20, indices{3, 1});
    end
end

function rep_data = get_repetition_data()
% [
%       repeated sym in prev frame, 
%       num repeated symbols in this frame,
%       num unique symbols
%       repetition_offset
% ]
    rep_data = cell(3, 2);
    rep_data{1, 1} = [0, 3, 3, 0];
    rep_data{2, 1} = [0, 2, 8, 4];
    rep_data{3, 1} = [1, 1, 8, 4];

    rep_data{1, 2} = [0, 2, 3, 0];
    rep_data{2, 2} = [0, 2, 4, 0];
    rep_data{3, 2} = [0, 2, 8, 4];
end

function [X_tx, bitstring1, bitstring2] = gen_multi_user_repetition(N_sc, L, num_frames, syms, user_indices, pilot_indices, IDFT_matrix)
% GEN_MULTI_USER_REPETITION
% Generates an OFDM signal for multiple users with pilot symbols 
% inserted in the first frame. Uses repetition coding to increase
% diversity.
%
% INPUTS:
%   N_sc            - Number of subcarriers in an OFDM frame
%   L               - Number of channel taps
%   num_frames      - Number of OFDM frames
%   syms            - Modulation symbol mapping
%   user_indices    - Cell array containing subcarrier indices allocated to
%                       each user in each frame
%   pilot_indices   - Indices of the inserted pilots
%   IDFT_matrix     - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx            - Serialized OFDM signal with all frames
%   bitstring1      - Random bitstring transmitted for user 1
%   bitstring2      - Random bitstring transmitted for user 2
%
    num_bits = log2(length(syms));

    cp_length = L - 1; % Cyclic prefix length
    pilot = 1; % Pilot symbol (usually average energy of constellation)

    % Generate random bit sequence for all users
    data_bits = num_bits * 15;
    bitstring1 = randi([0, 1], data_bits, 1);
    bitstring2 = randi([0, 1], data_bits, 1);

    rep_bitstring1 = [repmat(bitstring1(1:3*num_bits), 3, 1); repmat(bitstring1(3*num_bits + 1 : 7*num_bits), 2, 1); repmat(bitstring1(7*num_bits + 1:11*num_bits), 2, 1); bitstring1(11*num_bits + 1:end)];
    rep_bitstring2 = [repmat(bitstring2(1:3*num_bits), 2, 1); repmat(bitstring2(3*num_bits + 1 : 7*num_bits), 2, 1); repmat(bitstring2(7*num_bits + 1:11*num_bits), 2, 1); bitstring2(11*num_bits + 1:end)];
    
    % Initialize output signal
    X_tx = zeros((N_sc + cp_length) * num_frames, 1);
    
    bit_offset = [0, 0];
    for j = 1:num_frames
        X_freq = zeros(N_sc, 1);

        if j == 1
            X_freq(pilot_indices) = pilot;
        end

        [user1_idx, user2_idx] = user_indices{j, :};

        user1_bits = rep_bitstring1(bit_offset(1) + 1: bit_offset(1) + length(user1_idx) * num_bits);
        user1_bits = reshape(user1_bits, [num_bits, length(user1_idx)]);

        user2_bits = rep_bitstring2(bit_offset(2) + 1: bit_offset(2) + length(user2_idx) * num_bits);
        user2_bits = reshape(user2_bits, [num_bits, length(user2_idx)]);

        u1_syms = syms(bit2int(user1_bits, num_bits) + 1);
        u2_syms = syms(bit2int(user2_bits, num_bits) + 1);

        X_freq(user1_idx) = u1_syms;
        X_freq(user2_idx) = u2_syms;

        bit_offset(1) = bit_offset(1) + length(user1_idx)*num_bits;
        bit_offset(2) = bit_offset(2) + length(user2_idx)*num_bits;

        X_time = IDFT_matrix * X_freq;

        X_cp = [X_time(end - cp_length + 1:end); X_time];
        X_tx((N_sc + cp_length) * (j - 1) + 1 : (N_sc + cp_length) * j) = X_cp;

    end
end

function BER_avg = main_ofdm_sys_est_rep_multi_user(N_sc, L, num_frames, num_bits, syms, SNRs)
%MAIN_OFDM_SYS_EST_REP_MULTI_USER
% Simulates an OFDM system with multi-user, diversity (repetition coding), 
% and unknown channel. Implements channel estimation using Maximum 
% Likelihood (ML).
% Decodes symbols using MRC.
% 
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   num_bits    - Number of bits in one symbol
%   syms        - Modulation symbol mapping
%   SNRs        - SNR range to test the simulation
%
% OUTPUT:
%   BER_avg     - Avg Bit Error Rate for a range of SNRs for 2 UEs
%
    
    % Cyclic prefix length
    cp_length = L - 1;
    
    % BER initialization
    BER1 = zeros(size(SNRs));
    BER2 = zeros(size(SNRs));
    
    % Define pilot and subcarrier allocation
    pilot_indices = 1:4:N_sc; % Pilot subcarriers
    user_indices = get_user_indices(N_sc, num_frames, 1, 1);
    pilot = 1; % Pilot value
    
    % Get Repetition information
    rep_data = get_repetition_data();

    % Define DFT/IDFT matrices
    IDFT_matrix = exp(1i * 2 * pi * (0:N_sc-1)' * (0:N_sc-1) / N_sc) ./ sqrt(N_sc);
    DFT_matrix = sqrt(N_sc) * conj(IDFT_matrix);
    
    % Simulation parameters
    num_iters = 1e5;
    
    % Loop over SNR values
    for k = 1:length(SNRs)
        SNR_dB = SNRs(k);
        num_errors1 = 0;
        num_errors2 = 0;
        
        % Iterate over multiple simulation runs for averaging
        for m = 1:num_iters
            % Generate transmitted signal with repetition coding
            [x_tx, bitstring1, bitstring2] = gen_multi_user_repetition(N_sc, L, num_frames, syms, user_indices, pilot_indices, IDFT_matrix);
            
            % Channel simulation for both users
            y_noisy = zeros(length(x_tx), 2);
            for user = 1:2
                y_noisy(:, user) = simulate_channel(x_tx, N_sc, L, SNR_dB);
            end
            
            % Separate output by user and rehape by frames
            y_reshape = reshape(y_noisy, [N_sc + cp_length, num_frames, 2]);
            
            % Remove CP and get main subcarriers
            y_nocp = y_reshape(cp_length + 1: end, :, :);
            
            % Convert to freq domain
            Y_freq = zeros(size(y_nocp));
            for user = 1:2
                Y_freq(:, :, user) = DFT_matrix * y_nocp(:, :, user) ./ sqrt(N_sc);
            end
    
            % Estimate channel
            Z = DFT_matrix(pilot_indices, 1:L);
            Y_pilot = squeeze(Y_freq(pilot_indices, 1, :)); % Pilot observations
            hat_h = inv(Z' * Z) * Z' * Y_pilot ./ pilot; % Estimated channel taps         
            H_freq = (DFT_matrix(:, 1:L) * hat_h); % Full channel response
            
            % Decoding Block
            hatbitstring = zeros(length(bitstring1), 2);

            % Find repetitions of the symbols and use MRC for decoding
            total_sc = [0, 0];
            bit_offset = [0, 0];
            
            for frame = 1:num_frames
                for user = 1:2

                    % Get Indices for User
                    user_idx = user_indices{frame, user};
                    total_sc(user) = total_sc(user) + length(user_idx);

                    % Get frame repetition info for user
                    frame_rep = rep_data{frame, user};

                    % Get Estimated Channel and Received Symbols
                    H = H_freq(user_idx, user);
                    Y = Y_freq(user_idx, frame, user);
                    
                    % If repeated symbols are spaced apart in different
                    % frames, need to check previous frame last 4 symbols.
                    if frame_rep(1) && frame > 1

                        % Get previous frame info
                        prev_frame_idx = user_indices{frame - 1, user};
                        prev_frame_rep = rep_data{frame - 1, user};

                        % Channel and Symbols based on prev frame info
                        H_prev = H_freq(prev_frame_idx, user);
                        Y_prev = Y_freq(prev_frame_idx, frame  - 1, user);

                        % Combine last symbols of prev and starting symbols
                        % of current frame to extract repeating symbols
                        H_new = [H_prev(end - prev_frame_rep(4) + 1:end); H(1:end-frame_rep(4))];
                        Y_new = [Y_prev(end - prev_frame_rep(4) + 1:end); Y(1:end-frame_rep(4))];
                        
                        % Rearrange to matrix shape as defined.
                        num_reps = frame_rep(2) + frame_rep(1);
                        H_k = reshape(H_new, [length(H_new) / num_reps, num_reps]);
                        Y_k = reshape(Y_new, [length(Y_new) / num_reps, num_reps]);
                    else
                        % Extract repeating symbols and rearrange to matrix
                        % for ease. Symbols are reshaped to:
                        % |     S1      |       S1      |       S1      |
                        % |     S2      |       S2      |       S2      |
                        % |     S3      |       S3      |       S3      |
                        H_k = reshape(H(1:end-frame_rep(4)), [(length(user_idx) - frame_rep(4)) / frame_rep(2), frame_rep(2)]);
                        Y_k = reshape(Y(1:end-frame_rep(4)), [(length(user_idx) - frame_rep(4)) / frame_rep(2), frame_rep(2)]);
                    end
                    
                    % Equalize Channel (MRC) and Decode Symbols
                    R = sum((conj(H_k) .* Y_k), 2) ./ (vecnorm(H_k, 2, 2).^2);
                    decoded_bits = decode_symbols(R, num_bits);
                    
                    % Reconstruct bitstring
                    start_idx = bit_offset(user) + 1;
                    end_idx = bit_offset(user) + (frame_rep(3) - frame_rep(4)) * num_bits;
                    hatbitstring(start_idx:end_idx, user) = reshape(decoded_bits.', [1, (frame_rep(3) - frame_rep(4)) * num_bits]);
                    bit_offset(user) = end_idx;

                    % Final 4 symbols don't repeat => independent
                    % equalization
                    if (frame == 3)

                        % Get last four symbols (non-repeating)
                        H_k = H(end - frame_rep(4) + 1: end);
                        Y_k = Y(end - frame_rep(4) + 1: end);

                        % Equalize and decode
                        R = Y_k ./ H_k;
                        decoded_bits = decode_symbols(R, num_bits);
                        
                        % Reconstruct bitstring
                        start_idx = bit_offset(user) + 1;
                        end_idx = bit_offset(user) + (frame_rep(3) - frame_rep(4)) * num_bits;
                        hatbitstring(start_idx:end_idx, user) = reshape(decoded_bits.', [1, (frame_rep(3) - frame_rep(4)) * num_bits]);
                        bit_offset(user) = end_idx;
                    end
                end
            end

            % Separate reconstructed bitstrings
            hatbitstring1 = hatbitstring(:, 1);
            hatbitstring2 = hatbitstring(:, 2);
        

            % Accumulate iteration errors
            num_errors1 = num_errors1 + sum(hatbitstring1 ~= bitstring1, "all");
            num_errors2 = num_errors2 + sum(hatbitstring2 ~= bitstring2, "all");
        end

        % Individual BER of each user
        BER1(k) = num_errors1 / total_sc(1) / num_bits / num_iters;
        BER2(k) = num_errors2 / total_sc(2) / num_bits / num_iters;

    end

    % Avg BER of system
    BER_avg = (BER1 + BER2) ./ 2;
end

function [X_tx, bitstring, pilot_indices] = gen_trade_off(N_sc, L, num_frames, syms, num_pilots, IDFT_matrix)
% GEN_TRADE_OFF
% Generates an OFDM signal with pilot symbols inserted for channel estimation.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   num_pilots  - Number of Pilots to be inserted
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx            - Serialized OFDM signal with all frames
%   bitstring       - Random bitstring transmitted
%   pilot_indices   - Indices of the inserted pilots
%

    % Initialize constants and parameters
    num_bits = log2(length(syms)); % Number of bits per symbol
    cp_length = L - 1; % Cyclic prefix length

    % Generate random bit sequence
    num_data_bits = num_bits * (N_sc * num_frames - num_pilots);
    bitstring = randi([0, 1], num_data_bits, 1);
    pilot = 1; % Pilot symbol - average energy of constellation

    % Determine subcarrier indices for pilots and data
    pilot_indices = round(1:N_sc/num_pilots:N_sc);
    data_indices = setdiff(1:N_sc, pilot_indices);

    % Initialize output signal
    X_tx = zeros((N_sc + cp_length) * num_frames, 1);
    prev_frame_bits = 0; % Track bits processed across frames

    % Generate OFDM signal for each frame
    for frame_idx = 1:num_frames
        % Determine the number of data subcarriers for the current frame
        if frame_idx == 1
            num_data_subcarriers = length(data_indices); % Reserve space for pilots
        else
            num_data_subcarriers = N_sc; % All subcarriers used for data
        end

        % Extract bits for the current frame
        if num_data_subcarriers > 0
            start_idx = prev_frame_bits + 1;
            end_idx = prev_frame_bits + (num_data_subcarriers * num_bits);
            frame_bits = bitstring(start_idx:end_idx);
            frame_bits = reshape(frame_bits, [num_data_subcarriers, num_bits]);
           
            % Map bits to modulation symbols
            data_symbols = syms(bit2int(frame_bits.', num_bits) + 1).';
        end

        % Construct frequency-domain OFDM frame
        X_freq = zeros(N_sc, 1);
        if frame_idx == 1
            % Insert pilots and data in the first frame
            X_freq(pilot_indices) = pilot; % Pilot symbols
            if num_data_subcarriers > 0
                X_freq(data_indices) = data_symbols; % Data symbols
            end
        else
            % Use all subcarriers for data in subsequent frames
            X_freq = data_symbols;
        end

        % Transform to time domain using IDFT
        X_time = IDFT_matrix * X_freq;

        % Add cyclic prefix
        X_cp = [X_time(end - cp_length + 1:end); X_time];
        frame_start_idx = (frame_idx - 1) * (N_sc + cp_length) + 1;
        frame_end_idx = frame_idx * (N_sc + cp_length);
        X_tx(frame_start_idx:frame_end_idx) = X_cp;

        % Update processed bits tracker
        prev_frame_bits = prev_frame_bits + (num_data_subcarriers * num_bits);
    end
end

function BER = trade_off_system(N_sc, L, num_frames, num_bits, syms)
% TRADE_OFF_SYSTEM
% Simulates an OFDM system with single-user, no diversity, and unknown 
% channel. Implements channel estimation using Maximum Likelihood (ML).
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   num_bits    - Number of bits in one symbol
%   syms        - Modulation symbol mapping
%
% OUTPUT:
%   BER         - Bit Error Rate for 10dB SNR for all pilots
%

    % System Parameters
    cp_length = L - 1;          % Length of Cyclic Prefix
    pilot_range = 1:1:N_sc;     % Range of pilots from L to all subcarriers
    BER = zeros(size(pilot_range));     % BER results for pilots

    % Generate IDFT and DFT matrices
    IDFT_matrix = exp(1i * 2 * pi * (0:N_sc-1)' * (0:N_sc-1) / N_sc) / sqrt(N_sc);
    DFT_matrix = conj(IDFT_matrix);

    % Simulation Parameters
    SNR_dB = 10;
    num_iters = 1e5;
    pilot = 1;
    
    % Loop over each SNR value
    for num_pilots = pilot_range
            
        num_errors = 0;

        for m = 1:num_iters
            % Generate OFDM signal with pilot symbols
            [x, bitstring, pilot_indices] = gen_trade_off(N_sc, L, num_frames, syms, num_pilots, IDFT_matrix);
        
            % Simulate channel with random taps
            y_noisy = simulate_channel(x, N_sc, L, SNR_dB);
                
            % Initialize variables for decoding
            hatbitstring = zeros(length(bitstring), 1);
            H_freq = zeros(N_sc, L);
            bit_offset = 0;
    
            % Process each frame
            for j = 1:num_frames
                % Extract frame and remove CP
                frame_start = (N_sc + cp_length) * (j-1) + cp_length + 1;
                frame_end = (N_sc + cp_length) * j;
                y_no_cp = y_noisy(frame_start:frame_end);
    
                % Apply DFT to get frequency domain representation
                Y_freq = DFT_matrix * y_no_cp;   
                
                % Estimate channel response using ML
                if j == 1
                    % Separate pilot and data subcarriers
                    data_indices = setdiff(1:N_sc, pilot_indices);
                    
                    % Pilot matrix
                    Z = sqrt(N_sc) .* DFT_matrix(pilot_indices, 1:L); 
                    
                    % Pilot observations
                    Y_pilot = Y_freq(pilot_indices); 
                    
                    % Estimated channel taps
                    hat_h = inv(Z' * Z) * Z' * Y_pilot ./ pilot;
                    
                    % Full channel response
                    H_freq = sqrt(N_sc) .* (DFT_matrix(1:N_sc, 1:L) * hat_h);
                else
                    % All subcarriers contain data
                    data_indices = 1:N_sc;
                end
    
                % Channel Equalization
                R = Y_freq(data_indices) ./ H_freq(data_indices);
               
                % Decode symbols
                decoded_bits = decode_symbols(R, num_bits);
                
                % Reconstruct bitstring
                frame_bit_count = length(data_indices) * num_bits;
                hatbitstring(bit_offset + 1:bit_offset + frame_bit_count) = ...
                    reshape(decoded_bits, [frame_bit_count, 1]);
                bit_offset = bit_offset + frame_bit_count;
            end
            
            % Compute BER for the current SNR
            num_errors = num_errors + sum(hatbitstring ~= bitstring, "all");
        end
        BER(num_pilots) = num_errors / (N_sc * num_frames - num_pilots) / num_bits / num_iters;
    end
end
