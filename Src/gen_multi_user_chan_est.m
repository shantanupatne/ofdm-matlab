function [X_tx, bitstring1, bitstring2] = gen_multi_user_chan_est(N_sc, L, num_frames, syms, user_indices, pilot_indices, IDFT_matrix)
% GEN_MULTI_USER_CHAN_EST
% Generates an OFDM signal for multiple users with pilot symbols inserted in the first frame.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   num_users   - Number of users
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx         - Serialized OFDM signal with all frames
%   bitstring    - Random bitstring transmitted for each user
%   pilot        - Pilot symbol used
%   user_indices - Cell array containing subcarrier indices allocated to
%                  each user in each frame
%

    % Initialize constants
    num_bits = log2(length(syms)); % Number of bits per symbol
    cp_length = L - 1; % Cyclic prefix length
    pilot = 1; % Pilot symbol (usually average energy of constellation)

    % Generate random bit sequence for all users
    data_bits = num_bits * (N_sc - length(pilot_indices) + (num_frames - 1) * N_sc)/2;
    bitstring1 = randi([0, 1], data_bits, 1);
    bitstring2 = randi([0, 1], data_bits, 1);

    % Initialize output signal
    X_tx = zeros((N_sc + cp_length) * num_frames, 1);

    % Generate OFDM signal frame by frame
    bit_offset = [0, 0];
    for frame_idx = 1:num_frames
        X_freq = zeros(N_sc, 1); % Frequency-domain representation of the OFDM frame

        if frame_idx == 1
            % Insert pilot symbol
            X_freq(pilot_indices) = pilot;
        end

        % Get User Allocations
        [user_idx1, user_idx2] = user_indices{frame_idx, :};

        % Extract bits for the current frame
        % UE1
        start_idx1 = bit_offset(1) + 1;
        end_idx1 = bit_offset(1) + (length(user_idx1) * num_bits);
        frame_bits1 = bitstring1(start_idx1:end_idx1);
        frame_bits1 = reshape(frame_bits1, [length(user_idx1), num_bits]);

        % UE2
        start_idx2 = bit_offset(2) + 1;
        end_idx2 = bit_offset(2) + (length(user_idx2) * num_bits);
        frame_bits2 = bitstring2(start_idx2:end_idx2);
        frame_bits2 = reshape(frame_bits2, [length(user_idx2), num_bits]);

        % Map bits to modulation symbols
        data_symbols1 = syms(bit2int(frame_bits1.', num_bits) + 1).';
        data_symbols2 = syms(bit2int(frame_bits2.', num_bits) + 1).';

        % Insert user data symbols into allocated subcarriers
        X_freq(user_idx1) = data_symbols1;
        X_freq(user_idx2) = data_symbols2;

        % Update bit offset
        bit_offset(1) = end_idx1;
        bit_offset(2) = end_idx2;

        % Perform IDFT to convert to time domain
        time_domain_signal = IDFT_matrix * X_freq;

        % Add cyclic prefix
        cyclic_prefixed_signal = [time_domain_signal(end - cp_length + 1:end); time_domain_signal];
        frame_start_idx = (frame_idx - 1) * (N_sc + cp_length) + 1;
        frame_end_idx = frame_idx * (N_sc + cp_length);
        X_tx(frame_start_idx:frame_end_idx) = cyclic_prefixed_signal;
    end
end
