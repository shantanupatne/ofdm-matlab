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

