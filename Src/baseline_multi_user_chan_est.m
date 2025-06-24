function BER_avg = baseline_multi_user_chan_est(N_sc, L, num_frames, num_bits, syms, SNRs)
% BASELINE_MULTI_USER_CHAN_EST
% Simulates an OFDM system with multi-user, no diversity, and unknown 
% channel. Implements channel estimation using Maximum Likelihood (ML).
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
%   BER_avg     - Bit Error Rate for a range of SNRs
%

    % System Parameters
    cp_length = L - 1;           % Length of Cyclic Prefix
    BER1 = zeros(size(SNRs));
    BER2 = zeros(size(SNRs));
    
    % Subcarrier indices
    user_indices = get_user_indices(N_sc, num_frames, 1, 0); 
    
    % Pilot subcarrier indices
    pilot_indices = [1, 5, 9, 13, 17, 20];
    pilot = 1;

    % Generate IDFT and DFT matrices
    IDFT_matrix = exp(1i * 2 * pi * (0:N_sc-1)' * (0:N_sc-1) / N_sc) / sqrt(N_sc);
    DFT_matrix = conj(IDFT_matrix);
    num_iters = 1e5;

    % Loop over each SNR value
    for k = 1:length(SNRs)
        SNR_dB = SNRs(k);        
        num_errors1 = 0;
        num_errors2 = 0;
        for m = 1:num_iters
            % Generate OFDM signal with pilot symbols
            [x, bitstring1, bitstring2] = gen_multi_user_chan_est(N_sc, L, num_frames, syms, user_indices, pilot_indices, IDFT_matrix);
        
            % Simulate channel for each user
            y_noisy = zeros(length(x), 2);
            for l = 1:2
                y = simulate_channel(x, N_sc, L, SNR_dB);
                y_noisy(:, l) = y;
            end
            
            % Initialize variables for decoding
            hatbitstring1 = zeros(size(bitstring1));
            hatbitstring2 = zeros(size(bitstring2));
            bit_offset = 0;
            total_sc = [0, 0];

            % Process each frame
            for j = 1:num_frames
                % Extract frame and remove CP
                frame_start = (N_sc + cp_length) * (j-1) + cp_length + 1;
                frame_end = (N_sc + cp_length) * j;
                y_no_cp = y_noisy(frame_start:frame_end, :);
    
                % Apply DFT to get frequency domain representation
                Y1_freq = DFT_matrix * y_no_cp(:, 1);
                Y2_freq = DFT_matrix * y_no_cp(:, 2);
                
                % Estimate channel response using ML
                if j == 1
                    % Pilot matrix
                    Z = sqrt(N_sc) .* DFT_matrix(pilot_indices, 1:L);
    
                    % Pilot observations
                    Y1_pilot = Y1_freq(pilot_indices, :);
                    Y2_pilot = Y2_freq(pilot_indices, :);
    
                    % Estimated channel taps
                    hat_h1 = inv(Z' * Z) * Z' * Y1_pilot ./ pilot;
                    hat_h2 = inv(Z' * Z) * Z' * Y2_pilot ./ pilot;
    
                    % Full channel response
                    H1_freq = sqrt(N_sc) .* (DFT_matrix(1:N_sc, 1:L) * hat_h1);
                    H2_freq = sqrt(N_sc) .* (DFT_matrix(1:N_sc, 1:L) * hat_h2); 
                end
                
                % Get data indices for frame and user
                [user1_idx, user2_idx] = user_indices{j, :};
                total_sc(1) = total_sc(1) + length(user1_idx);
                total_sc(2) = total_sc(2) + length(user2_idx);
    
                % Channel equalization for data subcarriers
                R1 = Y1_freq(user1_idx) ./ H1_freq(user1_idx);
                R2 = Y2_freq(user2_idx) ./ H2_freq(user2_idx);
               
                % Decode symbols from equalized data
                decoded_bits1 = decode_symbols(R1, num_bits);
                decoded_bits2 = decode_symbols(R2, num_bits);
                
                % Reconstruct bitstring for the current frame
                frame_bit_count = length(user1_idx) * num_bits;
                hatbitstring1(bit_offset + 1:bit_offset + frame_bit_count) = ...
                    reshape(decoded_bits1, [frame_bit_count, 1]);
    
                hatbitstring2(bit_offset + 1:bit_offset + frame_bit_count) = ...
                    reshape(decoded_bits2, [frame_bit_count, 1]);
                bit_offset = bit_offset + frame_bit_count;
            end
    
            % Compute BER for the current SNR
            num_errors1 = num_errors1 + sum(hatbitstring1 ~= bitstring1, "all");
            num_errors2 = num_errors2 + sum(hatbitstring2 ~= bitstring2, "all");
        end
        BER1(k) = num_errors1 / total_sc(1) / num_bits / num_iters;
        BER2(k) = num_errors2 / total_sc(2) / num_bits / num_iters;
    end

    BER_avg = (BER1 + BER2) ./ 2;
end
