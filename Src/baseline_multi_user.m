function BER_avg = baseline_multi_user(N_sc, L, num_frames, num_bits, syms, SNRs)
% BASELINE_MULTI_USER
% This function simulates a basic ofdm system with no diversity, known 
% channel, and a single user. RANDOM ALLOCATION
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
%   BER_avg     - Bit Error Rate over all SNR range
%

    % System Parameters
    cp_length = L - 1;           % Length of Cyclic Prefix
    BER1 = zeros(size(SNRs));
    BER2 = zeros(size(SNRs));
    user_indices = get_user_indices(N_sc, num_frames, 0, 0); % Subcarrier indices

    % Generate IDFT and DFT matrices
    IDFT_matrix = exp(1i * 2 * pi * (0:N_sc-1)' * (0:N_sc-1) / N_sc) / sqrt(N_sc);
    DFT_matrix = conj(IDFT_matrix);
    num_iters = 1e5;

    % Loop over each SNR value
    for k = 1:length(SNRs)
        SNR_dB = SNRs(k);        
        
        num_error1 = 0;
        num_error2 = 0;

        for m = 1:num_iters
            % Generate OFDM signal with indices of each user's data
            [x, bitstring1, bitstring2] = gen_multi_user(N_sc, cp_length, num_frames, syms, user_indices, IDFT_matrix);
            
            % Simulate channel for each user
            y_noisy = zeros(length(x), 2);
            H_freq = zeros(N_sc, 2);
            for l = 1:2
                [y, h] = simulate_channel(x, N_sc, L, SNR_dB);
                H_freq(:, l) = sqrt(N_sc) .* DFT_matrix(1:N_sc, 1:L) * h; 
                y_noisy(:, l) = y;
            end
        
            % Preallocate reconstructed bitstring
            hatbitstring1 = zeros(size(bitstring1));
            hatbitstring2 = zeros(size(bitstring2));
            total_sc = [0, 0];
            
            % Retrieve data and compute BER
            for j = 1:num_frames
                % Remove CP
                y1_no_cp = y_noisy((N_sc + cp_length)*(j-1)+cp_length+1:(N_sc + cp_length)*(j), 1);
                y2_no_cp = y_noisy((N_sc + cp_length)*(j-1)+cp_length+1:(N_sc + cp_length)*(j), 2);
        
                % DFT at receiver
                Y1_freq = DFT_matrix * y1_no_cp;
                Y2_freq = DFT_matrix * y2_no_cp;  
                
                % Channel equalization
                R1 = Y1_freq ./ H_freq(:, 1);
                R2 = Y2_freq ./ H_freq(:, 2);
               
                % Decode symbols
                decoded_bits1 = decode_symbols(R1(user_indices{j, 1}), num_bits);
                decoded_bits2 = decode_symbols(R2(user_indices{j, 2}), num_bits);

                % Total subcarriers per user
                total_sc(1) = total_sc(1) + length(user_indices{j, 1});
                total_sc(2) = total_sc(2) + length(user_indices{j, 2});
                
                % Reconstruct bitstring per by user
                hatbitstring1(N_sc*num_bits/2*(j-1) + 1: N_sc*num_bits/2*j) = reshape(decoded_bits1, [N_sc*num_bits/2, 1]);
                hatbitstring2(N_sc*num_bits/2*(j-1) + 1: N_sc*num_bits/2*j) = reshape(decoded_bits2, [N_sc*num_bits/2, 1]);
            end
    
            % BER calculation
            num_error1 = num_error1 + sum(hatbitstring1 ~= bitstring1, "all");
            num_error2 = num_error2 + sum(hatbitstring2 ~= bitstring2, "all");
        end
        BER1(k) = num_error1 / total_sc(1) / num_bits / num_iters;
        BER2(k) = num_error2 / total_sc(2) / num_bits / num_iters;

    end

    BER_avg = (BER1 + BER2) ./ 2;
end

