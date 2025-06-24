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

