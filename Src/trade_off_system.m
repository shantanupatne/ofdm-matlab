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
