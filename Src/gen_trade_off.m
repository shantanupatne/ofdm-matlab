function [X_tx, bitstring, pilot_indices] = gen_trade_off(N_sc, L, num_frames, syms, num_pilots, IDFT_matrix)
% GEN_BASELINE_CHAN_EST
% Generates an OFDM signal with pilot symbols inserted for channel estimation.
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx        - Serialized OFDM signal with all frames
%   bitstring   - Random bitstring transmitted
%   pilot       - Pilot symbol used in the first frame
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
