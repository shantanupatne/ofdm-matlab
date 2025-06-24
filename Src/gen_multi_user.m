function [X_tx, bitstring1, bitstring2] = gen_multi_user(N_sc, cp_length, num_frames, syms, user_indices, IDFT_matrix)
% GEN_MULTI_USER
% Generates an OFDM signal for multiple users.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   cp_length   - Length of Cyclic Prefix
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   num_users   - Number of users
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx         - Serialized OFDM signal with all frames
%   user_indices - Cell array of user-specific subcarrier indices per frame
%   bitstring    - Bitstring transmitted per user (columns for each user)
%

    % Initialize constants and outputs
    num_bits = log2(length(syms)); % Bits per symbol
    X_tx = zeros((N_sc + cp_length) * num_frames, 1); % OFDM signal
    
    % Generate random bitstrings for all users
    user_bits = N_sc * num_frames * num_bits / 2;
    bitstring1 = randi([0, 1], user_bits, 1);
    bitstring2 = randi([0, 1], user_bits, 1);

    % Generate OFDM signal frame by frame
    for frame_idx = 1:num_frames
        X_freq = zeros(N_sc, 1); % Frequency-domain OFDM frame
        
        % Extract and reshape bitstring for current user and frame
        start_idx = (frame_idx - 1) * (N_sc * num_bits / 2) + 1;
        end_idx = frame_idx * (N_sc * num_bits / 2);
        
        % Get bits for UE1
        frame_bits1 = bitstring1(start_idx:end_idx, 1);
        frame_bits1 = reshape(frame_bits1, [N_sc / 2, num_bits]);

        % Get bits for UE2
        frame_bits2 = bitstring2(start_idx:end_idx, 1);
        frame_bits2 = reshape(frame_bits2, [N_sc / 2, num_bits]);

        % Map bits to modulation symbols
        mapped_symbols1 = syms(bit2int(frame_bits1.', num_bits) + 1).';
        mapped_symbols2 = syms(bit2int(frame_bits2.', num_bits) + 1).';
        
        % Assign symbols to allocated subcarriers
        X_freq(user_indices{frame_idx, 1}) = mapped_symbols1;
        X_freq(user_indices{frame_idx, 2}) = mapped_symbols2;
        
        % Perform IDFT and add cyclic prefix
        time_domain_signal = IDFT_matrix * X_freq;
        cyclic_prefixed_signal = [time_domain_signal(end - cp_length + 1:end); time_domain_signal];
        
        % Append the cyclic-prefixed signal to the serialized output
        start_idx = (frame_idx - 1) * (N_sc + cp_length) + 1;
        end_idx = frame_idx * (N_sc + cp_length);
        X_tx(start_idx:end_idx) = cyclic_prefixed_signal;
    end
end
