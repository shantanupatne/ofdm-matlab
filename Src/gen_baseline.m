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
